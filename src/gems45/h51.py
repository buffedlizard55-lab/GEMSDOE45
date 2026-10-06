"""H51 fields: multi-physics corroboration, Kaplan-Meier tip hazard, catalogue novelty.

THE SCIENTIFIC HYPOTHESIS (frozen in research/hypotheses.md before evaluation)
-----------------------------------------------------------------------------
An expert adds a fault that a previous mapper missed for one of three reasons, and each has a
different physical signature:

(a) **The trace exists but was not drawn** -- it is a scarp that is real, sharp and *linear* in
    detrended-elevation slope.  Signature: a ridge of ``|grad sigma|_det_elev_slope`` at 1-3 px.
(b) **The trace is buried** -- no topographic scarp, but a magnetic or gravity gradient edge lies
    on the *same line*, at the same strike.  Signature: an edge that two or more independent
    potential fields agree on, plus the topographic edge.
(c) **The trace is a continuation or a relay step-over of a mapped one** -- the mapped sections
    stop short.  Signature: an along-strike distance from a catalogue tip that is short relative
    to the catalogue's own relay-gap survival curve.

The prior family used (a) and (c) with a *guessed* extension radius and a *geometric* decimation
(Poisson-disk at d = 2.8 px).  This module replaces both: (c) becomes the Kaplan-Meier
product-limit curve of the catalogue's own tip-to-tip relay gaps, and the decimation becomes an
explicit expected-coverage maximisation (``coverage.py``).
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage

from .grid import BAND_NAMES

# Channels chosen for PHYSICS, not for measured catalogue enrichment.  Each is a documented
# fault/contact indicator: a derivative of a potential field or of detrended topography.
TOPO_CHANNELS = ["det_elev_slope"]
GRAV_CHANNELS = ["iso_grav_anom_hg", "iso_grav_anom_vg", "iso_grav_anom"]
MAG_CHANNELS = ["tmi_hg", "tmi_vg", "tc", "rtp"]
ALL_CHANNELS = TOPO_CHANNELS + GRAV_CHANNELS + MAG_CHANNELS
ELONG = 0.5          # Sato-style blob-suppression constant
SCALES = (1.2, 2.5)  # pixels; 120 m and 250 m -- the kernel scale and one octave up


def line_response(z: np.ndarray, sigma: float, elong: float = ELONG) -> np.ndarray:
    """Sato/Frangi-style ridge (line) response of a 2-D field: |l2| * exp(-l1^2 / (2 c^2 l2^2)).

    ``l1`` is the Hessian eigenvalue of larger magnitude (across-strike curvature) and ``l2`` the
    smaller one (along-strike).  A ridge has |l2| large and |l1| small, so the exponential
    suppresses blobs -- exactly the discrimination a fault trace needs against a lithologic
    contact or a hill.
    """
    s = ndimage.gaussian_filter(z.astype(np.float32), sigma)
    hyy = ndimage.gaussian_filter(z, sigma, order=(2, 0))
    hxx = ndimage.gaussian_filter(z, sigma, order=(0, 2))
    hxy = ndimage.gaussian_filter(z, sigma, order=(1, 1))
    trace = hxx + hyy
    det = hxx * hyy - hxy * hxy
    disc = np.sqrt(np.maximum(trace * trace / 4.0 - det, 0.0))
    l1 = trace / 2.0 + disc
    l2 = trace / 2.0 - disc
    return np.abs(l2) * np.exp(-(l1 * l1) / (2.0 * elong * elong * (l2 * l2 + 1e-30)))


def normalized_response(stack: np.ndarray, footprint: np.ndarray, name: str,
                        sigma: float, pct: float = 99.5) -> np.ndarray:
    """Per-channel ridge response scaled by its own high percentile inside the footprint."""
    i = BAND_NAMES.index(name)
    b = stack[i].astype(np.float32)
    v = b[footprint]
    z = (b - float(v.mean())) / (float(v.std()) + 1e-12)
    r = line_response(z, sigma)
    ref = float(np.percentile(r[footprint], pct)) if footprint.any() else 1.0
    return np.clip(r / (ref + 1e-30), 0.0, 1.0)


def corroboration(stack: np.ndarray, footprint: np.ndarray,
                  scales: tuple[float, ...] = SCALES) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Mean multi-physics line corroboration, plus per-physics-family means.

    Returns ``(mean_all, mean_topo, mean_potential)``.  The arithmetic mean is deliberate: the
    sibling H47 experiment used a *geometric* mean over 16 terms and measured only 1.2x
    enrichment, and a product over many terms is destroyed by any single near-zero factor.
    """
    acc = np.zeros(stack.shape[1:], dtype=np.float32)
    topo = np.zeros_like(acc)
    pot = np.zeros_like(acc)
    n = 0
    for sigma in scales:
        for name in ALL_CHANNELS:
            r = normalized_response(stack, footprint, name, sigma).astype(np.float32)
            acc += r
            if name in TOPO_CHANNELS:
                topo += r
            else:
                pot += r
            n += 1
    acc /= n
    topo /= (len(scales) * len(TOPO_CHANNELS))
    pot /= (len(scales) * len(GRAV_CHANNELS + MAG_CHANNELS))
    return acc, topo, pot


def novelty_profile(d_cat: np.ndarray, r_zero: float, tau: float) -> np.ndarray:
    """Catalogue-novelty weight eta(d) in [0, 1].

    ``eta`` is zero for ``d <= r_zero`` and rises linearly to 1 over the following ``tau`` px.
    The reference family's live dose-response (44,090 -> 40,199 -> 37,654 dots at 0.2600 ->
    0.2708 -> 0.2778) is the evidence that mass inside ~2 px of the mapped catalogue earns less
    than the marginal bar; ``r_zero`` is fitted, not assumed.  Above ~3 px the profile is flat:
    the organizers state that new-fault truth may lie within 300 m of a known trace, so a hard
    exclusion as wide as the scoring kernel would discard real credit.
    """
    d = np.asarray(d_cat, dtype=np.float32)
    return np.clip((d - r_zero) / max(tau, 1e-6), 0.0, 1.0)


def skeleton_endpoints(mask: np.ndarray) -> np.ndarray:
    """Skeleton pixels with exactly one 8-connected neighbour -> candidate fault tips."""
    k = np.ones((3, 3), dtype=np.uint8)
    k[1, 1] = 0
    nb = ndimage.convolve(mask.astype(np.uint8), k, mode="constant")
    return mask & (nb == 1)


def tip_continuation_field(shape, tips, survival, max_extend: int = 24) -> np.ndarray:
    """Kaplan-Meier tip hazard: P(the fault continues a further e px past this tip) = S(e).

    This is the brief's core instruction done literally.  The product-limit estimate ``S`` of the
    catalogue's own along-strike relay-gap distribution is a survival function, so ``S(e)`` *is*
    the probability that a fault is still going at distance ``e`` beyond its mapped tip.  The
    field is the maximum of that probability over all tips, painted only at integer steps along
    each tip's outward strike.
    """
    out = np.zeros(shape, dtype=np.float32)
    h, w = shape
    for t in tips:
        for e in range(1, max_extend + 1):
            p = float(survival(np.array([float(e)]))[0])
            if p <= 1e-6:
                break
            r = int(round(t.row + e * t.urow))
            c = int(round(t.col + e * t.ucol))
            if not (0 <= r < h and 0 <= c < w):
                break
            if p > out[r, c]:
                out[r, c] = p
    return out
