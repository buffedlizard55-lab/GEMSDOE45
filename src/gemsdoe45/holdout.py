"""Spatially-blocked proxy holdout for candidate ranking.

The competition's *true* test set is hidden (new faults not in the public
catalogue).  As an internal, fully-reproducible proxy we use the supplied
``existing_faults.tif`` (the public catalogue) as surrogate ground truth and
evaluate candidates on a spatially-blocked 4-fold holdout.  Four contiguous
quadrants of the survey are held out in turn; a 300 m (3 px) guard is eroded from
each block boundary so that predictions trained/derived in one block cannot leak
across the split line into the held-out block.

IMPORTANT (verified limitation, see registry note IRR-01 in the prior family): a
sibling project measured the correlation between this kind of catalogue-holdout
DTI and the live DrivenData leaderboard score at Spearman rho ~ 0.09 (p ~ 0.80),
i.e. essentially uncorrelated.  The holdout is therefore used as a *relative*
ranking and sanity instrument inside this sandbox, NOT as a predictor of the
official score.  We state this plainly rather than over-claim.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi

from . import metric


def four_spatial_fold_masks(footprint: np.ndarray, guard_pixels: int = 3):
    """Return 4 evaluation masks, one per contiguous quadrant (guard-eroded)."""
    fp = footprint.astype(bool)
    rows = np.any(fp, axis=1)
    cols = np.any(fp, axis=0)
    r0, r1 = int(np.argmax(rows)), int(fp.shape[0] - 1 - np.argmax(rows[::-1]))
    c0, c1 = int(np.argmax(cols)), int(fp.shape[1] - 1 - np.argmax(cols[::-1]))
    rm, cm = (r0 + r1) // 2, (c0 + c1) // 2
    quads = [
        (slice(r0, rm), slice(c0, cm)),
        (slice(r0, rm), slice(cm, c1 + 1)),
        (slice(rm, r1 + 1), slice(c0, cm)),
        (slice(rm, r1 + 1), slice(cm, c1 + 1)),
    ]
    masks = []
    struct = ndi.generate_binary_structure(2, 2)
    for sl in quads:
        m = np.zeros_like(fp)
        m[sl] = fp[sl]
        if guard_pixels > 0:
            m = ndi.binary_erosion(m, structure=struct, iterations=guard_pixels)
        masks.append(m)
    return masks


def _restrict(pred: np.ndarray, truth: np.ndarray, region: np.ndarray):
    """Restrict truth to region and zero pred outside region (guard handling)."""
    truth_r = (truth > 0) & region
    pred_r = np.where(region, np.nan_to_num(pred, nan=0.0), 0.0)
    return pred_r, truth_r


def evaluate_holdout(
    pred: np.ndarray,
    faults: np.ndarray,
    footprint: np.ndarray,
    guard_pixels: int = 3,
    return_folds: bool = False,
):
    """4-fold spatially-blocked proxy-holdout DTI.

    Returns dict with per-fold DTI, mean, std, and the all-data DTI reference.
    """
    masks = four_spatial_fold_masks(footprint, guard_pixels=guard_pixels)
    folds = []
    for m in masks:
        pred_r, truth_r = _restrict(pred, faults, m)
        dti, comp = metric.distance_weighted_tversky(pred_r, truth_r, return_components=True)
        folds.append({"dti": dti, **comp})
    all_dti, all_comp = metric.distance_weighted_tversky(pred, faults, return_components=True)
    out = {
        "n_folds": len(folds),
        "guard_pixels": guard_pixels,
        "fold_dti": [f["dti"] for f in folds],
        "mean_dti": float(np.mean([f["dti"] for f in folds])),
        "std_dti": float(np.std([f["dti"] for f in folds])),
        "alldata_dti": all_dti,
        "alldata_components": all_comp,
    }
    if return_folds:
        out["folds"] = folds
    return out
