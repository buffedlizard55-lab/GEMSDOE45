"""Which of the 19 supplied bands can localise a fault to the metric's 300 m kernel?

This module answers that question by measurement rather than by reading band names, because the
answer is not obvious and it is decisive for the design.  Three diagnostics are computed for every
band on the hash-pinned ``training_features.tif``:

1. ``rho_lag10`` -- Pearson correlation between a pixel and the pixel 10 px (1 km) away.  A band
   with rho close to 1 is regionally smooth and carries no information at the kernel's scale.
2. ``mean_abs_grad`` -- mean |d(band)/dy| in units of the band's own standard deviation per pixel.
3. ``auc_fault`` -- ROC AUC of the band value against the provided USGS/INGENIOUS catalogue
   (a *biased* reference: the catalogue is masked out of competition scoring, so a high AUC means
   "this band sees mapped faults", not "this band sees hidden faults").

The measured outcome (evidence/band_information.json) splits the 19 bands into
  * LOCATOR-class  (rho_lag10 < 0.90): rtp, tmi_hg, tmi_vg, tmi, det_elev_slope, iso_grav_anom_vg,
    mag_anom(0.944), iso_grav_anom_slope(0.910)
  * WEIGHT-class   (rho_lag10 ~> 0.95): geod_2ndinv, geod_shearrate, geod_dilaterate, ieq, deq,
    iso_grav_anom, depth_to_base_surf, cond_surf, det_elev
Only LOCATOR-class bands may place mass; WEIGHT-class bands may modulate it.  The three geodetic
strain-rate bands in particular are regional fields (rho(1 km) >= 0.98, mean |grad| <= 0.004
sigma/px), which is why a strain-rate "ridge detector" cannot satisfy the 300 m kernel.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage

from .grid import BAND_NAMES


def _auc(score: np.ndarray, label: np.ndarray, rng: np.random.Generator, n: int = 400_000) -> float:
    idx = rng.choice(score.size, size=min(n, score.size), replace=False)
    s, y = score[idx], label[idx]
    n1 = int(y.sum())
    n0 = int(len(y) - n1)
    if n1 == 0 or n0 == 0:
        return float("nan")
    order = np.argsort(s, kind="stable")
    ranks = np.empty(len(s), dtype=np.float64)
    ranks[order] = np.arange(len(s), dtype=np.float64)
    return float((ranks[y].sum() - n1 * (n1 - 1) / 2.0) / (n1 * n0))


def band_information(stack: np.ndarray, footprint: np.ndarray, known: np.ndarray,
                     seed: int = 20261006) -> list[dict]:
    rng = np.random.default_rng(seed)
    out = []
    for i, name in enumerate(BAND_NAMES):
        b = stack[i].astype(np.float64)
        v = b[footprint]
        mu, sd = float(v.mean()), float(v.std())
        z = (b - mu) / (sd + 1e-12)
        fy = footprint[:-10, 1:-1]
        r10 = float(np.corrcoef(z[:-10, 1:-1][fy], z[10:, 1:-1][fy])[0, 1])
        grad = np.zeros_like(z)
        grad[:-1, :] = np.abs(np.diff(z, axis=0))
        mg = float(grad[footprint].mean())
        auc_v = _auc(z[footprint], known[footprint], rng)
        gmag = np.hypot(*np.gradient(ndimage.gaussian_filter(z, 1.0)))
        auc_g = _auc(gmag[footprint], known[footprint], rng)
        lap = ndimage.laplace(ndimage.gaussian_filter(z, 1.0))
        auc_l = _auc(-np.abs(lap)[footprint], known[footprint], rng)
        out.append({
            "index": i + 1, "band": name,
            "rho_lag10px": round(r10, 4), "mean_abs_grad_sigma": round(mg, 5),
            "auc_value_vs_catalogue": round(auc_v, 4),
            "auc_gradmag_vs_catalogue": round(auc_g, 4),
            "auc_abs_laplace_vs_catalogue": round(auc_l, 4),
            "band_class": "LOCATOR" if r10 < 0.90 else ("INTERMEDIATE" if r10 < 0.95 else "WEIGHT"),
        })
    return out


LOCATOR_BANDS = ["rtp", "tmi_hg", "tmi_vg", "tmi", "det_elev_slope", "iso_grav_anom_vg",
                 "iso_grav_anom_hg", "mag_anom", "tc"]
WEIGHT_BANDS = ["geod_2ndinv", "geod_shearrate", "geod_dilaterate", "ieq_n100a15", "deq_n100a15",
                "iso_grav_anom", "depth_to_base_surf", "cond_surf", "det_elev",
                "iso_grav_anom_slope"]
