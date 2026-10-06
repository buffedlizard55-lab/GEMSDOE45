"""H51 emission: exact-decision-theory dot selection for the distance-weighted Tversky index.

WHY THIS MODULE EXISTS
----------------------
The official metric is

    DTI = TP_w / (TP_w + 0.2*FP_w + 0.8*FN_w)

with (page 967 of the competition, quote-verified) the 300 m triangular kernel

    k(d) = max(1 - d/R, 0),  R = 300 m = 3 px
    TP_w = sum_g max_{x: d(x,g)<=R} p(x) k(d(x,g))
    FP_w = sum_{x: p(x)>0} p(x) [1 - max_{g} k(d(x,g))]
    FN_w = sum_g [1 - max_{x: d(x,g)<=R} p(x) k(d(x,g))]

Three consequences are algebra, not modelling, and they drive the whole selection:

1. **Binary is optimal.**  Scaling a fixed support by lambda scales TP_w and FP_w by lambda and
   leaves the FN_w term's truth set unchanged, so ``DTI(lambda*p) <= DTI(1{p>0})``.  All 46
   live-scored owner rasters profiled in ``research/live_score_features.csv`` are exactly binary.

2. **The marginal bar is ``w(x) > 0.2*DTI``**, where ``w(x)`` is the expected realised kernel weight
   of a dot at x.  A dot whose expected weight is below that bar is *always* harmful; above it,
   *always* helpful.  At DTI = 0.2778 the bar is 0.0556, i.e. 283 m of hidden truth.

3. **Selection is a budgeted maximum-coverage problem with a submodular objective**, so greedy is
   its standard (1 - 1/e) approximation.  Greediness itself is not new: a sibling project
   (GEMSDOE28, arm H37-1) ran lazy-greedy coverage of its own detector field and reported that the
   gain did not survive a leave-one-fault-system-out far-field test (mean paired dDTI -0.000037).
   That negative result is recorded here rather than hidden, and it is why this artifact is
   published as an UNSCORED CANDIDATE.

MODEL
-----
``psi`` is the expected hidden-truth density (its integral is the modelled truth size ``K``).
The emulator keeps, exactly as the metric does,

    W(g) = max over placed dots x of k(|g - x|)        ->  T = sum_cells psi * W
    own(x) = sum_delta psi(x + delta) k(delta)         ->  F = #dots - sum over dots of own(x)

so its T and F are the metric's own quantities under a modelled truth, and the greedy rule below is
the metric's exact break-even for that model.  ``evidence/h51_live_fit.json`` measures how far the
model is from the live leaderboard; it does **not** predict it (Pearson +0.18 over 46 live-scored
rasters), so the emulator's absolute DTI is never quoted as a score and the stopping point is taken
from the live data instead (see ``dti_operating``).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .metric import R_PX


def kernel_offsets(r_px: float = R_PX) -> np.ndarray:
    """(N, 3) array of (dy, dx, k(d)) inside the triangular kernel support, k(d) = 1 - d/R."""
    rad = int(np.ceil(r_px))
    rows = []
    for dy in range(-rad, rad + 1):
        for dx in range(-rad, rad + 1):
            d = float(np.hypot(dy, dx))
            if d <= r_px:
                rows.append((dy, dx, 1.0 - d / r_px))
    return np.asarray(rows, dtype=np.float64)


OFFS = kernel_offsets()


def own_weight(psi: np.ndarray, rr: np.ndarray, cc: np.ndarray) -> np.ndarray:
    """own(x) = E[TP_w contributed by one dot at x] = sum_delta psi(x+delta) k(delta)."""
    h, w = psi.shape
    out = np.zeros(len(rr), dtype=np.float64)
    for dy, dx, kk in OFFS:
        r = rr + int(dy)
        c = cc + int(dx)
        ok = (r >= 0) & (r < h) & (c >= 0) & (c < w)
        out[ok] += kk * psi[np.clip(r[ok], 0, h - 1), np.clip(c[ok], 0, w - 1)]
    return out


def sparse_marginal(psi_c: np.ndarray, cr: np.ndarray, cc: np.ndarray,
                    W: np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    """Exact marginal coverage gain of a dot at every pixel, against coverage field W.

    gain(x) = sum_delta psi(x+delta) * max(k_delta - W(x+delta), 0).

    Only cells with ``psi > 0`` can contribute, so the sum runs over those cells and scatters back
    to the dot position that would cover them: O(29 * n_cells) instead of O(29 * H * W).
    """
    h, w = shape
    gain = np.zeros(h * w, dtype=np.float64)
    for dy, dx, kk in OFFS:
        deficit = kk - W[cr, cc]
        np.maximum(deficit, 0.0, out=deficit)
        tr = cr - int(dy)
        tc = cc - int(dx)
        ok = (tr >= 0) & (tr < h) & (tc >= 0) & (tc < w)
        gain += np.bincount(tr[ok] * w + tc[ok], weights=(psi_c * deficit)[ok], minlength=h * w)
    return gain


def marginal_bar(dti: float) -> float:
    """The metric's own marginal bar (src/gems45/metric.py::marginal_bar), re-exported."""
    return 0.2 * dti


@dataclass
class CoverageResult:
    dti_operating: float
    rows: np.ndarray
    cols: np.ndarray
    marginal: np.ndarray
    own: np.ndarray
    n_candidates: int
    n_selected: int
    dti_pred: float
    t_pred: float
    f_pred: float
    k_pred: float
    round_size: int
    stopped_by: str
    history: list = field(default_factory=list)

    def as_dict(self) -> dict:
        return dict(n_candidates=self.n_candidates, n_selected=self.n_selected,
                    dti_operating=self.dti_operating,
                    dti_pred=round(self.dti_pred, 6), t_pred=round(self.t_pred, 3),
                    f_pred=round(self.f_pred, 3), k_pred=round(self.k_pred, 3),
                    round_size=self.round_size, stopped_by=self.stopped_by,
                    history=self.history)


def greedy_coverage(psi: np.ndarray, candidates: np.ndarray, K: float,
                    max_dots: int, round_size: int = 1000,
                    dti_operating: float | None = None) -> CoverageResult:
    """Round-batched greedy maximum coverage of the modelled truth.

    Stopping has two modes.

    * ``dti_operating is None`` -- the emulator's self-consistent break-even
      ``gain * (0.2*F + 0.8*K) > 0.2 * T * (1 - own)`` evaluated on its own evolving T and F.
    * ``dti_operating = d`` -- the same algebra with the operating point pinned to a *measured*
      value: ``dT*(1 - 0.2*d) > 0.2*d*(1 - own)``.  This mode exists because the emulator's
      absolute T is not trustworthy (Pearson +0.18 against 46 live scores, and a predicted DTI of
      0.83 where the live top is 0.3345).  The ranking of gains is useful; the scale is not.
    """
    psi = np.ascontiguousarray(psi, dtype=np.float64)
    h, w = psi.shape
    if candidates.shape != (h, w):
        raise ValueError("psi and candidates must share shape")
    if abs(float(psi.sum()) - K) > 1e-6 * max(K, 1.0):
        raise ValueError("psi must integrate to K")
    allowed = np.asarray(candidates, dtype=bool)
    allowed_flat = allowed.ravel()
    cr, cc = np.nonzero(allowed & (psi > 0))
    psi_c = psi[cr, cc]
    W = np.zeros((h, w), dtype=np.float64)
    T = F = 0.0
    rows: list[int] = []
    cols: list[int] = []
    mar: list[float] = []
    own_l: list[float] = []
    hist: list[dict] = []
    stopped = "cap"

    while len(rows) < max_dots:
        gain = sparse_marginal(psi_c, cr, cc, W, (h, w))
        take = min(round_size, max_dots - len(rows))
        # dots may only be placed on candidate cells: the candidate set *is* the design's
        # hypothesis space (tip continuations, splays, ridge-corroborated, off the mapped flank),
        # and a dot anywhere else would silently break those guarantees.
        nz = np.flatnonzero((gain > 0.0) & allowed_flat)
        if nz.size == 0:
            stopped = "fully-covered"
            break
        sel = (nz[np.argpartition(-gain[nz], take - 1)[:take]] if nz.size > take
               else nz[np.argsort(-gain[nz])])
        accepted = 0
        for idx in sel:
            r, c = divmod(int(idx), w)
            g = float(gain[idx])
            if g <= 0.0:
                continue
            own = float(sum(kk * psi[r + int(dy), c + int(dx)]
                            for dy, dx, kk in OFFS
                            if 0 <= r + int(dy) < h and 0 <= c + int(dx) < w))
            own_capped = min(1.0, own)          # E[max match] is a weight; it cannot exceed 1
            if dti_operating is not None:
                if g * (1.0 - 0.2 * float(dti_operating)) <= \
                        0.2 * float(dti_operating) * (1.0 - own_capped):
                    continue
            elif len(rows) > 0 and (0.2 * F + 0.8 * K) > 0 and \
                    g * (0.2 * F + 0.8 * K) <= 0.2 * T * (1.0 - own_capped):
                continue
            for dy, dx, kk in OFFS:
                y, x = r + int(dy), c + int(dx)
                if 0 <= y < h and 0 <= x < w and kk > W[y, x]:
                    W[y, x] = kk
            F += 1.0 - own_capped
            rows.append(r)
            cols.append(c)
            mar.append(g)
            own_l.append(own)
            accepted += 1
        # T is recomputed exactly at the end of every round from the coverage field W, so it can
        # never exceed the modelled truth size K.  Per-dot `gain` values are exact against the W of
        # the previous round but stale inside a round; T is not.
        T = float((psi_c * W[cr, cc]).sum())
        if accepted and len(rows) // 5000 != (len(rows) - accepted) // 5000:
            hist.append(dict(n=len(rows), T=round(T, 2), F=round(F, 2),
                             dti=round(T / (0.2 * T + 0.2 * F + 0.8 * K), 6)))
        if accepted == 0:
            stopped = "marginal-bar" if dti_operating is not None else "break-even"
            break
    D = 0.2 * T + 0.2 * F + 0.8 * K
    return CoverageResult(
        dti_operating=float(dti_operating) if dti_operating is not None else float("nan"),
        rows=np.asarray(rows, dtype=np.int32), cols=np.asarray(cols, dtype=np.int32),
        marginal=np.asarray(mar, dtype=np.float64), own=np.asarray(own_l, dtype=np.float64),
        n_candidates=int(allowed_flat.sum()), n_selected=len(rows),
        dti_pred=float(T / D) if D > 0 else 0.0, t_pred=float(T), f_pred=float(F), k_pred=float(K),
        round_size=int(round_size), stopped_by=stopped, history=hist)


def score_emission(dots_mask: np.ndarray, psi: np.ndarray, K: float) -> dict:
    """Evaluate an ARBITRARY existing dot set under the same emulator (used for live calibration)."""
    psi = np.asarray(psi, dtype=np.float64)
    rr, cc = np.nonzero(dots_mask)
    own = np.minimum(own_weight(psi, rr, cc), 1.0)
    W = np.zeros_like(psi)
    for r, c in zip(rr, cc):
        for dy, dx, kk in OFFS:
            y, x = r + int(dy), c + int(dx)
            if 0 <= y < psi.shape[0] and 0 <= x < psi.shape[1] and kk > W[y, x]:
                W[y, x] = kk
    T = float((psi * W).sum())
    F = float(len(rr) - own.sum())
    D = 0.2 * T + 0.2 * F + 0.8 * K
    return dict(n=int(len(rr)), T=T, F=F, K=float(K), t_over_k=T / K,
                dti_pred=float(T / D) if D > 0 else 0.0)
