"""Memory-bounded feature extraction for the supervised lineament detector.

The grid is 12,279,160 cells and the detector uses 37 feature planes.  Materialising all of them at
once needs about 1.8 GB, which this sandbox does not have, so the extractor streams one feature at a
time and either (a) samples it at pre-chosen pixel indices, or (b) accumulates its contribution to a
running linear predictor and discards it.  Both paths keep at most three 12.3 M-cell float64 planes
alive (about 300 MB).
"""
from __future__ import annotations

from collections.abc import Iterator

import numpy as np
from scipy import ndimage

from .grid import BAND_NAMES

GRAD_BANDS = ["det_elev_slope", "iso_grav_anom_hg", "iso_grav_anom_slope", "iso_grav_anom",
              "geod_2ndinv", "cond_surf", "tmi_hg", "rtp", "tmi_vg", "det_elev", "tc",
              "iso_grav_anom_vg", "depth_to_base_surf", "geod_dilaterate"]
VAL_BANDS = ["det_elev_slope", "cond_surf", "depth_to_base_surf", "iso_grav_anom", "det_elev",
             "geod_2ndinv", "tc"]
SIGMAS = (1.2, 2.5)


def feature_names() -> list[str]:
    names = [f"|grad{s}|_{n}" for s in SIGMAS for n in GRAD_BANDS]
    names += [f"val_{n}" for n in VAL_BANDS]
    names += ["|laplace|_det_elev_slope", "laplace_det_elev_slope"]
    return names


def _z(b: np.ndarray, footprint: np.ndarray) -> np.ndarray:
    v = b[footprint]
    return (b - v.mean()) / (v.std() + 1e-12)


def iter_features(stack: np.ndarray, footprint: np.ndarray) -> Iterator[tuple[str, np.ndarray]]:
    """Yield ``(name, full_grid_float64_plane)`` one at a time."""
    idx = {n: i for i, n in enumerate(BAND_NAMES)}
    for s in SIGMAS:
        for n in GRAD_BANDS:
            z = ndimage.gaussian_filter(_z(stack[idx[n]].astype(np.float64), footprint), s)
            gy, gx = np.gradient(z)
            yield f"|grad{s}|_{n}", np.hypot(gx, gy)
            del z, gy, gx
    for n in VAL_BANDS:
        yield f"val_{n}", _z(stack[idx[n]].astype(np.float64), footprint)
    ds = ndimage.gaussian_filter(_z(stack[idx["det_elev_slope"]].astype(np.float64), footprint), 1.2)
    lap = ndimage.laplace(ds)
    yield "|laplace|_det_elev_slope", np.abs(lap)
    yield "laplace_det_elev_slope", lap


def sample_matrix(stack: np.ndarray, footprint: np.ndarray, flat_idx: np.ndarray) -> np.ndarray:
    """Feature matrix of shape ``(len(flat_idx), n_features)`` without ever holding a full stack."""
    names = feature_names()
    X = np.empty((len(flat_idx), len(names)), dtype=np.float32)
    for j, (nm, plane) in enumerate(iter_features(stack, footprint)):
        X[:, j] = plane.ravel()[flat_idx]
        if j % 8 == 0:
            del plane
    return X


def full_field(stack: np.ndarray, footprint: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Sigmoid linear predictor over the whole grid, accumulated feature by feature."""
    z = np.full(stack.shape[1:], float(w[-1]), dtype=np.float64)
    for j, (nm, plane) in enumerate(iter_features(stack, footprint)):
        z += w[j] * plane
        del plane
    out = 1.0 / (1.0 + np.exp(-z))
    return out.astype(np.float32)
