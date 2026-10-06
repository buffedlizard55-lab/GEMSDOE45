"""Emission of a fault-prediction raster from tip extensions.

The emitter is deliberately simple and auditable: it paints integer-valued dots at full
probability (1.0) and nothing anywhere else.  Two reasons.

1.  The published metric charges ``alpha = 0.2`` per unit of predicted probability mass that is not
    covered by the kernel, independently of the value's magnitude: ``FP_w = sum_x p(x)(1 - w(x))``.
    Emitting 1.0 rather than, say, 0.4 on the same pixel therefore does not change the false-positive
    penalty at all when the pixel is uncovered -- it is charged once either way -- while it strictly
    increases the credit that pixel can collect, because ``TP_w`` uses ``max_x p(x) k(d)``.  For a
    fixed *support*, the all-ones prediction dominates every other assignment in [0, 1].
2.  A graded field is not readable by eye against the catalogue, which matters because the whole
    artifact is meant to be reviewed by a geologist before it is trusted.
"""
from __future__ import annotations

import numpy as np


def emit_tip_extensions(shape: tuple[int, int], tips, lengths: np.ndarray,
                        strike_gate: np.ndarray | None = None,
                        gate_min: float = 0.0) -> np.ndarray:
    """Paint dots along each tip's outward strike for ``min(1, e) / 1`` up to its own extension.

    ``strike_gate`` (optional) is a per-pixel [0, 1] line-strength field; a dot is only painted when
    ``strike_gate[r, c] >= gate_min``.
    """
    out = np.zeros(shape, dtype=np.float32)
    h, w = shape
    for t, L in zip(tips, lengths):
        for e in range(1, int(L) + 1):
            r = int(round(t.row + e * t.urow))
            c = int(round(t.col + e * t.ucol))
            if not (0 <= r < h and 0 <= c < w):
                break
            if strike_gate is not None and strike_gate[r, c] < gate_min:
                continue
            out[r, c] = 1.0
    return out


def emit_step_over_bridges(shape: tuple[int, int], tips, gaps: np.ndarray, events: np.ndarray,
                           lengths: np.ndarray, match: np.ndarray | None = None) -> np.ndarray:
    """Paint the relay-ramp segment between each tip and its along-strike counterpart.

    For a tip whose next along-trend structure is at distance G, the two faults are predicted to
    link at the midpoint, so the bridge is drawn over the central half of the gap --
    ``[0.25 G, 0.75 G]`` -- rather than the whole gap, which would double the mass for no extra
    credit.  Bridge pixels are painted at the same value as tip pixels.
    """
    out = np.zeros(shape, dtype=np.float32)
    h, w = shape
    for t, L, G, ev in zip(tips, lengths, gaps, events):
        if ev != 1:
            continue
        a, b = 0.25 * G, 0.75 * G
        e = int(np.floor(a))
        while e <= int(np.ceil(b)):
            r = int(round(t.row + e * t.urow))
            c = int(round(t.col + e * t.ucol))
            if 0 <= r < h and 0 <= c < w:
                out[r, c] = 1.0
            e += 1
    return out


def thin_by_distance(mask: np.ndarray, min_sep: float = 3.0, order: np.ndarray | None = None) -> np.ndarray:
    """Keep pixels so that no two survivors are closer than ``min_sep`` px (Poisson-disk).

    ``order`` gives the priority (higher first); ties are broken by raster order.  This is the
    metric's own exclusion radius: two predictions closer than R = 3 px can only cover the same
    300 m disc, so the second one adds false-positive mass for no new credit.
    """
    from scipy.spatial import cKDTree
    rr, cc = np.nonzero(mask)
    if len(rr) == 0:
        return mask
    if order is None:
        order = mask[rr, cc]
    pri = order[rr, cc] if order.shape == mask.shape else np.asarray(order)
    idx = np.lexsort((-pri, cc, rr))
    pts = np.stack([rr[idx], cc[idx]], axis=1).astype(np.float64)
    tree = cKDTree(pts)
    keep = np.zeros(len(pts), dtype=bool)
    taken = []
    for i in range(len(pts)):
        if not taken:
            keep[i] = True
            taken.append(pts[i])
            continue
        d, _ = tree.query(pts[i], k=1, distance_upper_bound=min_sep - 1e-9)
        # explicit re-check against accepted points only
        acc = np.asarray(taken)
        if np.min(np.hypot(acc[:, 0] - pts[i, 0], acc[:, 1] - pts[i, 1])) >= min_sep:
            keep[i] = True
            taken.append(pts[i])
    out = np.zeros_like(mask, dtype=np.float32)
    out[rr[idx[keep]], cc[idx[keep]]] = 1.0
    return out
