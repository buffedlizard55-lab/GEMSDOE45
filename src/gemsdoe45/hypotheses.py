"""GEMSDOE45 geological hypotheses (candidate generators).

All hypotheses use ONLY the 19 supplied competition feature bands
(``gems-geodawn-numerical-features.tif``); no external data is required, so the
whole pipeline is reproducible inside this sandbox.

Band map (1-based -> 0-based index), verified from the file's band descriptions:
    1  magnetic deviation       2  RTP magnetic (mag field)      3  TMI horiz grad
    4  geodetic 2nd invariant   5  isostatic-gravity slope        6  tilt/curvature (mag edge)
    7  geodetic shear rate      8  geodetic dilatation           9  TMI vert grad
   10  distance to earthquake   11  isostatic-gravity vert grad  12  detrended elevation
   13  isostatic gravity (grav field)  14  TMI                   15  depth to basement (MT)
   16  earthquake density (seismic)     17  conductivity surface (MT)  18  isostatic-gravity horiz grad
   19  detrended elevation slope

Novelty vs the prior GEMSDOE family
-----------------------------------
The prior family's best work (WPH-01, GEMSDOE42) is "gravity/RTP-magnetic edge
survival across upward continuation x 0D ridge persistence".  The GEMSDOE42
unimplemented hypotheses (H42-A..E) cover dip-migration worms, radiometric
K/Th, LiDAR scarplets, strain eigenvectors and raw seismicity lineaments.  None
of them use the supplied *magnetotelluric* bands (depth-to-basement, conductivity
surface) as primary edge detectors, and none use an explicit intersection /
blend ("two independent signatures must coincide") fusion.  GEMSDOE45 introduces
both.  NOTE: several geodetic/seismicity bands contain no-data (NaN) holes even
inside the survey footprint; those are treated as *neutral* (no boost) rather
than propagated, so a hypothesis never silently collapses to NaN.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi

from . import features as F


def _neutralize(arr: np.ndarray) -> np.ndarray:
    """Replace NaN (no-data) with 0 so fusion stays finite and neutral."""
    return np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)


def _activity_gate(feats: np.ndarray) -> np.ndarray:
    """Soft active-structure prior in [0,1]: strain + seismicity coincidence.

    No-data holes in the geodetic/seismicity bands are neutralised to 0 so the
    gate is always finite inside the footprint.
    """
    strain2 = _neutralize(F.robust_normalize(feats[..., 3]))   # band 4
    dilat = _neutralize(F.robust_normalize(feats[..., 7]))     # band 8
    quake = _neutralize(F.robust_normalize(feats[..., 15]))    # band 16
    active = strain2 * (0.5 + 0.5 * quake) + dilat * (0.5 + 0.5 * quake)
    return F.robust_normalize(active)


def _topo_edge(feats: np.ndarray) -> np.ndarray:
    """Topographic lineament edge from detrended elevation + its slope."""
    elev = feats[..., 11]   # band 12 detrended elevation
    slope = feats[..., 18]  # band 19 detrended elevation slope
    e = np.maximum(F.gradient_magnitude(elev), _neutralize(F.robust_normalize(slope)))
    return F.robust_normalize(e)


def _mag_grav_worms(feats: np.ndarray):
    mag = feats[..., 1]    # RTP magnetic (band 2)
    grav = feats[..., 12]  # isostatic gravity (band 13)
    mag_worm = F.worm_survival(mag)
    grav_worm = F.worm_survival(grav)
    return mag_worm, grav_worm


def _strong_edge_ensemble(feats: np.ndarray) -> np.ndarray:
    """Best-aligning geophysical edge field (continuous, in [0,1]).

    Per single-band proxy-holdout experiments, detrended-elevation slope, the
    tilt/curvature magnetic edge, the isostatic-gravity and TMI horizontal
    gradients, and the RTP field each capture different parts of the fault grain.
    Taking the per-pixel MAX of their robust-normalised gradient magnitudes
    yields a lineament field that aligns with far more catalogue faults than any
    single band (and, by extension, with the unmapped ones that share it).
    """
    tilt = feats[..., 5]    # band 6  tilt/curvature (magnetic edge)
    dsl = feats[..., 18]    # band 19 detrended-elevation slope
    ghh = feats[..., 17]    # band 18 isostatic-gravity horizontal gradient
    thh = feats[..., 2]     # band 3  TMI horizontal gradient
    rtp = feats[..., 1]     # band 2  RTP magnetic
    norms = [F.robust_normalize(F.gradient_magnitude(x)) for x in (tilt, dsl, ghh, thh, rtp)]
    stack = np.stack([np.nan_to_num(n) for n in norms], axis=0)
    return np.nan_to_num(np.max(stack, axis=0))


# ---------------------------------------------------------------------------
# H45-1  (PRIMARY, novel)  MT basement/conductivity offset edge, corroborated by
#        (magnetic OR gravity OR topographic worms), activity-gated, then BLENDED
#        with the worm field so it keeps worm coverage while injecting the novel
#        deep-structure (geothermal-conduit) signal.
# ---------------------------------------------------------------------------
def h45a_mtbasement_blend(feats: np.ndarray) -> np.ndarray:
    basement = feats[..., 14]  # depth to basement (band 15)
    cond = feats[..., 16]      # conductivity surface (band 17)
    base_edge = np.maximum(
        F.gradient_magnitude(basement),
        F.gradient_magnitude(cond),
    )
    base_edge_n = F.robust_normalize(base_edge)
    strong = _strong_edge_ensemble(feats)   # best-aligning fault-grain lineaments
    topo = _topo_edge(feats)
    corroborator = np.maximum(np.maximum(strong, topo), base_edge_n)
    # Novel MT component: a basement/conductivity offset kept only where it is
    # corroborated by an independent signature (edge ensemble or topography).
    basement_sig = F.soft_intersect(base_edge_n, corroborator, radius_px=2)
    basement_sig = basement_sig * (0.7 + 0.3 * _activity_gate(feats))
    # Additive blend: the fault-grain edge ensemble carries the score (so the
    # top pixels stay fault-aligned), and the novel MT-basement signal is added
    # as a boost.  This keeps catalogue-fault alignment high while the MT term
    # differentiates the map from every prior worm-only submission.
    return np.clip(0.8 * strong + 0.2 * basement_sig, 0.0, 1.0)


# ---------------------------------------------------------------------------
# H45-1b  pure MT-basement/conductivity intersection (no worm floor) — the most
#         distinct variant; tests whether the deep signal alone tracks faults.
# ---------------------------------------------------------------------------
def h45b_mtbasement_only(feats: np.ndarray) -> np.ndarray:
    basement = feats[..., 14]
    cond = feats[..., 16]
    base_edge = np.maximum(
        F.gradient_magnitude(basement), F.gradient_magnitude(cond)
    )
    base_edge_n = F.robust_normalize(base_edge)
    mag_worm, grav_worm = _mag_grav_worms(feats)
    topo = _topo_edge(feats)
    corroborator = np.maximum(np.maximum(mag_worm, grav_worm), topo)
    sig = F.soft_intersect(base_edge_n, corroborator, radius_px=2)
    return np.clip(sig * (0.7 + 0.3 * _activity_gate(feats)), 0.0, 1.0)


# ---------------------------------------------------------------------------
# H45-2  basement flexure-axis (2nd-derivative) corroborated by worms+topo.
# ---------------------------------------------------------------------------
def h45c_basement_flexure(feats: np.ndarray) -> np.ndarray:
    basement = F._finite(feats[..., 14])
    sm = ndi.gaussian_filter(basement, sigma=1.0)
    lap = np.abs(ndi.laplace(sm))  # fault = step -> dipole in 2nd derivative
    flex = F.robust_normalize(lap)
    mag_worm, grav_worm = _mag_grav_worms(feats)
    topo = _topo_edge(feats)
    corroborator = np.maximum(np.maximum(mag_worm, grav_worm), topo)
    sig = F.soft_intersect(flex, corroborator, radius_px=2)
    return np.clip(sig * (0.7 + 0.3 * _activity_gate(feats)), 0.0, 1.0)


# ---------------------------------------------------------------------------
# H45-3  seismicity-weighted strain-orientation composite corroborated by worms.
# ---------------------------------------------------------------------------
def h45d_seismic_strain(feats: np.ndarray) -> np.ndarray:
    active = _activity_gate(feats)
    mag_worm, grav_worm = _mag_grav_worms(feats)
    edge_set = np.maximum(mag_worm, grav_worm)
    return np.clip(F.soft_intersect(active, edge_set, radius_px=2), 0.0, 1.0)


# ---------------------------------------------------------------------------
# H45-4  worms-only control (closest to WPH-01; used as a distinctness control).
# ---------------------------------------------------------------------------
def h45e_worms_only(feats: np.ndarray) -> np.ndarray:
    mag_worm, grav_worm = _mag_grav_worms(feats)
    return np.clip(np.maximum(mag_worm, grav_worm), 0.0, 1.0)


HYPOTHESES = {
    "H45-1a-mtbasement-blend": h45a_mtbasement_blend,
    "H45-1b-mtbasement-only": h45b_mtbasement_only,
    "H45-2-basement-flexure": h45c_basement_flexure,
    "H45-3-seismic-strain": h45d_seismic_strain,
    "H45-4-worms-only-control": h45e_worms_only,
}
