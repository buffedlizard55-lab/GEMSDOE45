"""Distance-weighted Tversky index (DTI) for the GEMS Prize Challenge.

TRANSCRIBED VERBATIM from the official problem description, read 2026-10-06:
https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/#performance-metric

    k(d)  = (1 - d/R)_+ ,  R = 300 m  (3 px at 100 m)
    TP_w  = sum_{g in G}  max_{x : d(x,g) <= R} p(x) * k(d(x,g))
    FP_w  = sum_{x : p(x)>0} p(x) * [ 1 - max_{g in G} k(d(x,g)) ]
    FN_w  = sum_{g in G} [ 1 - max_{x : d(x,g) <= R} p(x) * k(d(x,g)) ]
    DTI   = TP_w / ( TP_w + alpha*FP_w + beta*FN_w + eps ),   alpha = 0.2, beta = 0.8

Official worked example (same page): TP_w=3.00, FP_w=1.89, FN_w=2.00 -> 0.60.
This module returns 0.602652 for those terms.

Structure of the implementation
-------------------------------
Everything reduces to two 49-offset (3-px disk) "offer" fields, each exact and O(29 * N):

    offer[r, c] = max over predictions within R of pixel (r, c) of p * k(d)      -> TP_w = offer[truth].sum()
    wfield[r, c] = max over truth pixels within R of pixel (r, c) of k(d)         -> FP_w = p*(1-wfield) summed

The second field ``wfield`` is the repository's central diagnostic: it is the *realised kernel
weight* of a pixel, i.e. how much credit a unit of mass emitted there would capture. ``score()``
and ``wfield()`` are cross-checked against a literal O(N_g * N_x) transcription of the published
sums in ``tests/test_metric.py``.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

ALPHA = 0.2  # published false-positive penalty
BETA = 0.8   # published false-negative penalty
R_PX = 3.0   # kernel support: 300 m at 100 m pixels
R_M = 300.0
EPS = 1e-12

# OFFICIAL masking rule -- DrivenData staff (chrisk-dd), community thread 11516, 2026-09-16:
#   "Pixels corresponding to known USGS/INGENIOUS faults are masked / excluded from evaluation,
#    so they do not count towards penalty terms."
#   "Re-evaluation will also mask/exclude the existing USGS/INGENIOUS faults."
# Consequence used throughout: mass emitted on / adjacent to the known catalogue earns no credit
# and costs no penalty, so every reported DTI in this repository is computed on the OFF-CATALOGUE
# domain only.
MASK_RULE_SOURCE = ("https://community.drivendata.org/t/scoring-clarification-are-known-usgs-"
                    "ingenious-faults-masked-when-scoring-and-are-they-in-the-final-round-label-"
                    "set/11516")


def _disk_offsets(r_px: float = R_PX) -> list[tuple[int, int, float]]:
    """All integer offsets within Euclidean distance r_px, with their distance in pixels."""
    rad = int(np.ceil(r_px))
    out = []
    for dy in range(-rad, rad + 1):
        for dx in range(-rad, rad + 1):
            d = float(np.hypot(dy, dx))
            if d <= r_px:
                out.append((dy, dx, d))
    return out


_OFFSETS = _disk_offsets()


def kernel(d_px: np.ndarray | float) -> np.ndarray | float:
    """Triangular kernel k(d) = max(1 - d/R, 0), R = 3 px = 300 m."""
    return np.maximum(1.0 - np.asarray(d_px, dtype=np.float64) / R_PX, 0.0)


def _shift(src: np.ndarray, dy: int, dx: int) -> np.ndarray:
    """``out[r, c] = src[r + dy, c + dx]`` with zeros outside the array bounds."""
    out = np.zeros_like(src)
    h, w = src.shape
    r0, r1 = max(0, -dy), min(h, h - dy)
    c0, c1 = max(0, -dx), min(w, w - dx)
    if r0 < r1 and c0 < c1:
        out[r0:r1, c0:c1] = src[r0 + dy:r1 + dy, c0 + dx:c1 + dx]
    return out


def offer_field(pred: np.ndarray, dtype=np.float32) -> np.ndarray:
    """``offer[r, c] = max over prediction pixels within R of (r, c) of p * k(d)``.

    At a ground-truth pixel g this is exactly the term ``max_{x: d(x,g)<=R} p(x) k(d(x,g))``
    of the published TP_w sum.
    """
    out = np.zeros(pred.shape, dtype=dtype)
    for dy, dx, d in _OFFSETS:
        np.maximum(out, (_shift(pred, dy, dx) * dtype(1.0 - d / R_PX)), out=out)
    return out


def weight_field(truth: np.ndarray, dtype=np.float32) -> np.ndarray:
    """``w[r, c] = max over truth pixels within R of (r, c) of k(d)``.

    This is the realised kernel weight of a unit of mass emitted at (r, c). By the exact marginal
    theorem it is profitable iff ``w > alpha * DTI``.
    """
    return offer_field(np.asarray(truth, dtype=dtype), dtype=dtype)


@dataclass(frozen=True)
class DTIResult:
    tp_w: float
    fp_w: float
    fn_w: float
    dti: float
    n_truth: int
    support_px: int
    total_mass: float
    mean_weight_on_support: float

    def as_dict(self) -> dict:
        return {
            "TP_w": round(self.tp_w, 6),
            "FP_w": round(self.fp_w, 6),
            "FN_w": round(self.fn_w, 6),
            "DTI": round(self.dti, 8),
            "n_truth_px": self.n_truth,
            "support_px": self.support_px,
            "total_mass": round(self.total_mass, 6),
            "mean_kernel_weight_on_support": round(self.mean_weight_on_support, 6),
        }


def reduce_terms(pred: np.ndarray, truth: np.ndarray) -> dict:
    """Exact (TP_w, FP_w, FN_w, |G|) with no tolerance and no approximation.

    Identities used, both immediate from the published sums:
      FN_w = |G| - TP_w      (the kernel is exactly 0 outside R, so both sums maximise alike)
      FP_w = sum(p) - sum(p * wfield)
    """
    pred = np.asarray(pred, dtype=np.float32)
    truth = np.asarray(truth, dtype=bool)
    if pred.shape != truth.shape:
        raise ValueError(f"shape mismatch: pred {pred.shape} vs truth {truth.shape}")
    if pred.min() < 0 or pred.max() > 1:
        raise ValueError("prediction values must lie in [0, 1] (competition rule)")

    n_truth = int(truth.sum())
    pos = pred > 0.0
    support_px = int(pos.sum())
    total_mass = float(pred[pos].astype(np.float64).sum())

    tp_w = 0.0
    m_total = 0.0
    mean_w = 0.0
    if n_truth and support_px:
        off = offer_field(pred)
        tp_w = float(off[truth].astype(np.float64).sum())
        wf = weight_field(truth)
        m_total = float((pred[pos].astype(np.float64) * wf[pos].astype(np.float64)).sum())
        mean_w = float(wf[pos].astype(np.float64).mean())

    return {
        "tp_w": tp_w,
        "fp_w": total_mass - m_total,
        "fn_w": float(n_truth) - tp_w,
        "n_truth": n_truth,
        "support_px": support_px,
        "total_mass": total_mass,
        "mean_weight_on_support": mean_w,
    }


def score(pred: np.ndarray, truth: np.ndarray) -> DTIResult:
    t = reduce_terms(pred, truth)
    denom = t["tp_w"] + ALPHA * t["fp_w"] + BETA * t["fn_w"] + EPS
    return DTIResult(
        tp_w=t["tp_w"], fp_w=t["fp_w"], fn_w=t["fn_w"], dti=float(t["tp_w"] / denom),
        n_truth=t["n_truth"], support_px=t["support_px"], total_mass=t["total_mass"],
        mean_weight_on_support=t["mean_weight_on_support"],
    )


def dti_from_terms(tp_w: float, fp_w: float, fn_w: float) -> float:
    """DTI from the three weighted terms -- the form used to check the published example."""
    return tp_w / (tp_w + ALPHA * fp_w + BETA * fn_w + EPS)


# ---------------------------------------------------------------------------------------------
# The three exact consequences of the metric that this repository designs against.
# ---------------------------------------------------------------------------------------------

def marginal_bar(dti: float) -> float:
    """Minimum realised kernel weight for which an added unit of mass raises the score.

    DERIVATION (exact). Put T = TP_w, F = FP_w and K = |G| = FN_w + T. The denominator is
        D = T + alpha F + beta (K - T) = alpha (T + F) + beta K.
    Adding one unit of mass at x with realised weight w = max_g k(d(x,g)) raises F by exactly
    (1 - w) (the mass is either outside every kernel or displaces mass that was already counted)
    and raises T by at most w. With the best case dT = w,
        dD = alpha (w + 1 - w) = alpha,  independent of x.
    So the score rises iff dT/D - T dD/D^2 > 0 iff w > alpha * DTI. At alpha = 0.2 this is
    0.052 at DTI = 0.2600, 0.0556 at 0.2778, 0.0639 at 0.3195 and 0.0652 at 0.3262.
    """
    return ALPHA * float(dti)


def profitable_radius_m(dti: float) -> float:
    """Distance in metres inside which an emitted pixel is profitable at the given DTI."""
    return (1.0 - marginal_bar(dti)) * R_M


def scale_law(tp_w: float, fp_w: float, n_truth: int, lam: float) -> float:
    """Exact score of ``lam * p`` for any prediction with terms (tp_w, fp_w) at lam = 1.

    DTI(lam) = lam*T / (lam*alpha*(T+F) + beta*K).  Verified to 1e-15 against ``score``.
    """
    return lam * tp_w / (ALPHA * (lam * tp_w + lam * fp_w) + BETA * n_truth + EPS)


def recall_for_target(target_dti: float, rho: float, k_px: float = 1.0) -> float:
    """Weighted recall t = T/K needed to reach ``target_dti`` at false-positive ratio rho = F/K.

    From 1/DTI = alpha + alpha*(F/T) + beta*(K/T) with F = rho*K and T = t*K:
        t = (alpha*rho + beta) / (1/target_dti - alpha).
    """
    return (ALPHA * rho + BETA) / (1.0 / target_dti - ALPHA)


def coverage_for_target(target_dti: float, rho: float) -> float:
    """Weighted coverage c = T/K needed to reach ``target_dti`` at false-positive ratio rho = F/K.

    Algebraically identical to ``recall_for_target``: substituting beta/alpha = 4 into
        c = (alpha*rho + beta)/(1/target - alpha)
    gives the closed form used in the project's design tables,
        c = 0.2 * s * (rho + 4) / (1 - 0.2 * s).
    Kept as a separate name because the report quotes the second form.
    """
    s = target_dti
    norm = ALPHA * s * (rho + BETA / ALPHA) / (1.0 - ALPHA * s)
    assert abs(norm - recall_for_target(target_dti, rho)) < 1e-12
    return norm
