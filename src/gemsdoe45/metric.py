"""Distance-weighted Tversky index (DTI) — exact competition metric.

Reference (verbatim equations from the official problem page
https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/#performance-metric):

Let p(x) in [0,1] be the predicted fault probability and g(x) the ground-truth
label at pixel x.  The triangular kernel is k(d) = max(1 - d/R, 0) with support
R = 300 m (== 3 pixels at 100 m resolution).

    TP_w = sum_{g in G}  max_{x : d(x,g) <= R}  p(x) * k(d(x,g))
    FP_w = sum_{x : p(x) > 0}  p(x) * [ 1 - max_{g in G} k(d(x,g)) ]
    FN_w = sum_{g in G}  [ 1 - max_{x : d(x,g) <= R} p(x) * k(d(x,g)) ]

    DTI(a,b) = TP_w / ( TP_w + a*FP_w + b*FN_w + eps )

with a = 0.2, b = 0.8.

Key algebraic simplifications used here (exact, not approximations):
  * Because k(d) is monotonically decreasing in d, max_{g} k(d(x,g)) is achieved
    at the nearest ground-truth pixel.  Hence the FP_w bracket equals
    min(1, d_nearest(x)/R_px) where d_nearest is the Euclidean pixel distance to
    the nearest truth pixel.  This lets FP_w be computed with a single distance
    transform instead of an O(predicted x truth) loop.
  * For each truth pixel g, let best_g = max_{x} p(x) k(d(x,g)).  Then
    TP_w = sum_g best_g  and  FN_w = N_truth - TP_w  (the kernel support is only
    3 px, so best_g is computed exactly over the 7x7 disk around each g).
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi

# Competition-fixed constants (official problem page).
ALPHA = 0.2  # false-positive penalty
BETA = 0.8   # false-negative penalty
RADIUS_M = 300.0  # triangular kernel support
PIXEL_SIZE_M = 100.0  # competition grid resolution
EPS = 1e-9


def _disk_offsets(radius_px: float):
    """Return ((di, dj), k) pairs for the Euclidean disk of given pixel radius."""
    r = int(np.ceil(radius_px))
    out = []
    for di in range(-r, r + 1):
        for dj in range(-r, r + 1):
            d = np.hypot(di, dj)
            if d <= radius_px + 1e-9:
                k = max(1.0 - d / radius_px, 0.0)
                out.append((di, dj, k))
    return out


def distance_weighted_tversky(
    pred: np.ndarray,
    truth: np.ndarray,
    pixel_size_m: float = PIXEL_SIZE_M,
    radius_m: float = RADIUS_M,
    alpha: float = ALPHA,
    beta: float = BETA,
    eps: float = EPS,
    return_components: bool = False,
):
    """Compute DTI between a predicted probability map and binary truth.

    Parameters
    ----------
    pred : 2D float array
        Predicted probabilities in [0,1].  NaN / non-finite values are treated as
        0 (outside the competition footprint there is no prediction).
    truth : 2D array
        Binary ground-truth labels (truth > 0 == fault pixel).
    pixel_size_m, radius_m : floats
        Resolution and kernel support.  radius_m / pixel_size_m defines the disk
        radius in pixels (3.0 for the competition).
    alpha, beta, eps : floats
        Metric coefficients (0.2 / 0.8 / 1e-9 per the competition).

    Returns
    -------
    dti : float
        The distance-weighted Tversky index, or 0.0 if undefined.
    components : dict (only if return_components=True)
        {'tp_w', 'fp_w', 'fn_w', 'dti', 'n_truth'}
    """
    pred = np.asarray(pred, dtype=np.float64)
    truth = np.asarray(truth, dtype=np.float64)
    if pred.shape != truth.shape:
        raise ValueError(f"pred shape {pred.shape} != truth shape {truth.shape}")

    pred = np.nan_to_num(pred, nan=0.0, posinf=0.0, neginf=0.0)
    truth_mask = truth > 0
    n_truth = int(truth_mask.sum())
    if n_truth == 0:
        # No ground truth: DTI is undefined -> return 0 (nothing to discover).
        comp = dict(tp_w=0.0, fp_w=float(pred.sum()), fn_w=0.0, dti=0.0, n_truth=0)
        return (0.0, comp) if return_components else 0.0

    radius_px = radius_m / pixel_size_m
    offsets = _disk_offsets(radius_px)

    # --- TP_w and FN_w: per-truth-pixel max over the 3-px disk ---------------
    gi, gj = np.nonzero(truth_mask)
    best = np.zeros(n_truth, dtype=np.float64)
    h, w = pred.shape
    for di, dj, k in offsets:
        # predicted value sampled at the position that would contribute to truth
        # pixel (gi, gj) from offset (di, dj); i.e. predicted pixel (gi+di, gj+dj)
        pi = gi + di
        pj = gj + dj
        inside = (pi >= 0) & (pi < h) & (pj >= 0) & (pj < w)
        vals = np.zeros(n_truth, dtype=np.float64)
        vals[inside] = pred[pi[inside], pj[inside]] * k
        best = np.maximum(best, vals)
    tp_w = float(best.sum())
    fn_w = float(n_truth - tp_w)

    # --- FP_w: distance transform to nearest truth pixel ----------------------
    # distance_transform_edt on (truth==0) yields Euclidean px distance to the
    # nearest truth pixel for every pixel.
    dist_to_truth = ndi.distance_transform_edt(~truth_mask)
    fp_bracket = np.minimum(1.0, dist_to_truth / radius_px)
    fp_w = float(np.sum(pred * fp_bracket))  # pred>0 already implied (pred>=0)

    denom = tp_w + alpha * fp_w + beta * fn_w + eps
    dti = tp_w / denom if denom > 0 else 0.0

    comp = dict(tp_w=tp_w, fp_w=fp_w, fn_w=fn_w, dti=dti, n_truth=n_truth)
    return (dti, comp) if return_components else dti


def self_test() -> bool:
    """Verify the metric against analytic ground-truth cases. No external data."""
    # Case 1: perfect prediction (pred == truth line) -> DTI == 1.0
    g = np.zeros((21, 21))
    g[:, 10] = 1.0
    p = g.copy()
    dti1, c1 = distance_weighted_tversky(p, g, return_components=True)
    assert abs(dti1 - 1.0) < 1e-9, f"perfect case failed: {dti1}"

    # Case 2: empty prediction -> DTI == 0.0 (TP=0, FP=0)
    dti2, c2 = distance_weighted_tversky(np.zeros_like(g), g, return_components=True)
    assert abs(dti2 - 0.0) < 1e-9, f"empty case failed: {dti2}"

    # Case 3: prediction shifted 1 px away from truth line.
    p3 = np.zeros_like(g)
    p3[:, 11] = 1.0
    dti3, c3 = distance_weighted_tversky(p3, g, return_components=True)
    # truth pixel at distance 1px -> k=1-1/3=0.667 -> TP_w ~ N*0.667, FN ~ N*0.333,
    # FP: predicted line at dist 1 -> bracket 1/3 -> FP_w ~ N*1*(1/3)
    # DTI = 0.667N / (0.667N + 0.2*(0.333N) + 0.8*(0.333N))
    #     = 0.667 / (0.667 + 0.0667 + 0.2667) = 0.667/1.0 = 0.667
    assert abs(dti3 - 0.6667) < 1e-2, f"shift-1 case failed: {dti3} comp={c3}"

    # Case 4: all-ones prediction (broad field) -> DTI much lower than focused.
    p4 = np.ones_like(g)
    dti4, c4 = distance_weighted_tversky(p4, g, return_components=True)
    assert dti4 < dti3, f"broad field should score lower than focused line: {dti4} vs {dti3}"

    # Case 5: a probability of 0.5 on the line (perfect coverage, half confidence).
    # TP_w = 0.5N, FN_w = N - 0.5N = 0.5N, FP_w = 0 (predicted line sits on truth).
    # DTI = 0.5N / (0.5N + 0.2*0 + 0.8*0.5N) = 0.5 / 0.9 = 0.5556.
    p5 = g.copy() * 0.5
    dti5, c5 = distance_weighted_tversky(p5, g, return_components=True)
    assert abs(dti5 - 0.5556) < 1e-2, f"half-conf case failed: {dti5}"

    return True


if __name__ == "__main__":
    ok = self_test()
    print("metric self_test:", "PASS" if ok else "FAIL")
