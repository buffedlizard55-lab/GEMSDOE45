"""Geophysical edge / lineament operators used by the GEMSDOE45 hypotheses.

These are clean-room, dependency-light implementations of standard geophysical
edge detectors.  The "worm" survival operator is a documented, pragmatic
approximation of the upward-continuation worming method (Archibald, Gow &
Boschetti, 1999; Hornby et al., 1999): a real fault produces an edge that
*persists* across many smoothing (upward-continuation) scales, whereas a shallow
cultural or volcanic anomaly flickers out at the first smoothing.  We approximate
continuation heights with Gaussian smoothing scales and count how many scales a
pixel is an edge — a pixel that is an edge at many scales is a "worm".
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi


def robust_normalize(arr: np.ndarray, low: float = 1.0, high: float = 99.0) -> np.ndarray:
    """Scale finite values to [0,1] using percentile clipping. NaN -> NaN."""
    a = np.asarray(arr, dtype=np.float64)
    finite = a[np.isfinite(a)]
    if finite.size == 0:
        return np.zeros_like(a)
    lo = np.percentile(finite, low)
    hi = np.percentile(finite, high)
    if hi <= lo:
        hi = lo + 1.0
    out = np.clip((a - lo) / (hi - lo), 0.0, 1.0)
    out[~np.isfinite(a)] = np.nan
    return out


def _finite(arr: np.ndarray) -> np.ndarray:
    a = np.asarray(arr, dtype=np.float64)
    a = np.nan_to_num(a, nan=0.0, posinf=0.0, neginf=0.0)
    return a


def gradient_magnitude(field: np.ndarray, sigma: float = 0.0) -> np.ndarray:
    """Total horizontal gradient magnitude (edge strength) of a field."""
    f = _finite(field)
    if sigma > 0:
        f = ndi.gaussian_filter(f, sigma=sigma)
    gy, gx = np.gradient(f)  # note: np.gradient returns (rows, cols)
    return np.hypot(gx, gy)


def worm_survival(
    field: np.ndarray,
    sigmas=(0.0, 1.0, 2.0, 3.0, 4.0),
    percentile: float = 98.0,
) -> np.ndarray:
    """Multiscale edge-survival ("worm") map in [0,1].

    For each smoothing scale the field's gradient magnitude is thresholded at a
    high percentile to mark edge pixels; the per-pixel *count* of scales at which
    it is an edge is normalised to [0,1].  Persistent edges -> ~1, flickering
    edges -> low.
    """
    f = _finite(field)
    total = np.zeros_like(f, dtype=np.float64)
    for s in sigmas:
        sm = ndi.gaussian_filter(f, sigma=s) if s > 0 else f
        mag = gradient_magnitude(sm, sigma=0.0)
        thr = np.percentile(mag, percentile)
        if not np.isfinite(thr) or thr <= 0:
            continue
        total += (mag >= thr).astype(np.float64)
    return total / max(1, len(sigmas))


def dilate(mask: np.ndarray, radius_px: int) -> np.ndarray:
    """Binary dilation of a thresholded mask by a pixel radius."""
    if radius_px <= 0:
        return mask.astype(bool)
    struct = ndi.generate_binary_structure(2, 2)
    return ndi.binary_dilation(mask > 0, structure=struct, iterations=radius_px)


def soft_intersect(a: np.ndarray, b: np.ndarray, radius_px: int = 2) -> np.ndarray:
    """a is kept only where b (a corroborating edge) is within radius_px.

    Returns a * b_dilated, i.e. the primary signal is retained only where it is
    corroborated by the secondary signal within a small neighbourhood.  This is
    the geometric expression of "a fault is where two independent geophysical
    signatures coincide", which suppresses isolated false positives.
    """
    a = np.nan_to_num(a, nan=0.0)
    b = np.nan_to_num(b, nan=0.0)
    bd = dilate(b > 0, radius_px).astype(np.float64)
    return a * bd


def thin_to_skeleton(prob: np.ndarray, threshold: float = 0.5):
    """Binarize and morphologically thin to a 1-pixel skeleton (lineament)."""
    from skimage.morphology import skeletonize, remove_small_objects
    from skimage.filters import threshold_otsu

    p = np.nan_to_num(prob, nan=0.0)
    if threshold is None:
        try:
            threshold = float(threshold_otsu(p[p > 0])) if (p > 0).any() else 0.5
        except Exception:
            threshold = 0.5
    binary = p >= threshold
    binary = remove_small_objects(binary, min_size=4)
    if not binary.any():
        return binary
    try:
        skel = skeletonize(binary)
    except Exception:
        skel = binary
    return skel
