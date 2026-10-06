"""Verification of src/gems45/metric.py against the published metric.

The brute-force functions below are a literal, uninspired transcription of the four formulas on
https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/#performance-metric
(sum over ground-truth pixels of a max over predicted pixels, sum over predicted pixels of a max
over ground-truth pixels, triangular kernel evaluated per pixel pair with a real distance). They
are deliberately O(N_truth * N_pred) so that no vectorisation trick can hide an error.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from gems45 import metric as M  # noqa: E402


def k_bruteforce(d_px: float, R: float = 3.0) -> float:
    return max(1.0 - d_px / R, 0.0)


def dti_terms_bruteforce(pred: np.ndarray, truth: np.ndarray) -> tuple[float, float, float]:
    """Literal transcription of TP_w, FP_w, FN_w. No shortcuts."""
    h, w = pred.shape
    gt = [(r, c) for r in range(h) for c in range(w) if truth[r, c]]
    px = [(r, c, float(pred[r, c])) for r in range(h) for c in range(w) if pred[r, c] > 0]

    # TP_w = sum_{g in G} max_{x: d(x,g) <= R} p(x) * k(d(x,g))
    tp_w = 0.0
    for (gr, gc) in gt:
        best = 0.0
        for (xr, xc, pv) in px:
            d = float(np.hypot(xr - gr, xc - gc))
            if d <= 3.0:
                best = max(best, pv * k_bruteforce(d))
        tp_w += best

    # FP_w = sum_{x: p(x) > 0} p(x) * [1 - max_{g in G} k(d(x,g))]
    fp_w = 0.0
    for (xr, xc, pv) in px:
        mk = 0.0
        for (gr, gc) in gt:
            d = float(np.hypot(xr - gr, xc - gc))
            mk = max(mk, k_bruteforce(d))
        fp_w += pv * (1.0 - mk)

    # FN_w = sum_{g in G} [1 - max_{x: d(x,g) <= R} p(x) * k(d(x,g))]
    fn_w = 0.0
    for (gr, gc) in gt:
        best = 0.0
        for (xr, xc, pv) in px:
            d = float(np.hypot(xr - gr, xc - gc))
            if d <= 3.0:
                best = max(best, pv * k_bruteforce(d))
        fn_w += 1.0 - best

    return tp_w, fp_w, fn_w


def test_official_published_example():
    """The organizers' own worked example: 3.00 / (3.00 + 0.2*1.89 + 0.8*2.00) = 0.60."""
    v = M.dti_from_terms(3.00, 1.89, 2.00)
    assert abs(v - 0.60) < 0.005, v
    assert abs(v - 0.602652) < 1e-6, v


def test_bruteforce_equivalence_random():
    rng = np.random.default_rng(20261006)
    for trial in range(12):
        shape = (int(rng.integers(8, 16)), int(rng.integers(8, 16)))
        truth = rng.random(shape) < 0.12
        if truth.sum() == 0:
            truth[0, 0] = True
        pred = np.where(rng.random(shape) < 0.18, rng.random(shape), 0.0).astype(np.float32)
        bt, bf, bn = dti_terms_bruteforce(pred, truth)
        got = M.reduce_terms(pred, truth)
        assert abs(got["tp_w"] - bt) < 1e-5, (trial, got["tp_w"], bt)
        assert abs(got["fp_w"] - bf) < 1e-5, (trial, got["fp_w"], bf)
        assert abs(got["fn_w"] - bn) < 1e-5, (trial, got["fn_w"], bn)


def test_bruteforce_equivalence_graded():
    """Graded (non-binary) predictions must also agree -- the metric is defined for p in [0,1]."""
    rng = np.random.default_rng(7)
    truth = np.zeros((11, 13), dtype=bool)
    truth[3, 2:9] = True
    truth[7, 5:12] = True
    pred = np.zeros((11, 13), dtype=np.float32)
    pred[3, 3:8] = rng.random(5).astype(np.float32)
    pred[7, 4:11] = rng.random(7).astype(np.float32)
    pred[5, 4:9] = 0.3
    bt, bf, bn = dti_terms_bruteforce(pred, truth)
    got = M.reduce_terms(pred, truth)
    assert abs(got["tp_w"] - bt) < 1e-5
    assert abs(got["fp_w"] - bf) < 1e-5
    assert abs(got["fn_w"] - bn) < 1e-5


def test_identity_fn_equals_g_minus_tp():
    rng = np.random.default_rng(11)
    truth = rng.random((20, 20)) < 0.1
    pred = (rng.random((20, 20)) < 0.15).astype(np.float32)
    got = M.reduce_terms(pred, truth)
    assert abs(got["fn_w"] - (got["n_truth"] - got["tp_w"])) < 1e-9


def test_perfect_prediction_scores_one():
    truth = np.zeros((20, 20), dtype=bool)
    truth[5, 3:15] = True
    r = M.score(truth.astype(np.float32), truth)
    assert abs(r.dti - 1.0) < 1e-6, r.dti


def test_empty_prediction_scores_zero():
    truth = np.zeros((20, 20), dtype=bool)
    truth[5, 3:15] = True
    r = M.score(np.zeros((20, 20), dtype=np.float32), truth)
    assert r.dti == 0.0
    assert r.tp_w == 0.0


def test_scale_law_matches_direct_scoring():
    """DTI(lam*p) = lam*T / (lam*alpha*(T+F) + beta*K), exactly."""
    rng = np.random.default_rng(3)
    truth = rng.random((24, 24)) < 0.08
    pred = (rng.random((24, 24)) < 0.12).astype(np.float32)
    base = M.reduce_terms(pred, truth)
    for lam in (0.05, 0.5, 1.0):
        direct = M.score((pred * lam).astype(np.float32), truth).dti
        analytic = M.scale_law(base["tp_w"], base["fp_w"], base["n_truth"], lam)
        assert abs(direct - analytic) < 1e-6, (lam, direct, analytic)


def test_marginal_theorem_holds():
    """The marginal rule w > alpha*DTI predicts, exactly, whether adding mass helps."""
    rng = np.random.default_rng(5)
    truth = np.zeros((40, 40), dtype=bool)
    truth[20, 5:35] = True
    base = np.zeros((40, 40), dtype=np.float32)
    base[20, 6:20] = 1.0
    s0 = M.score(base, truth).dti
    bar = M.marginal_bar(s0)
    wf = M.weight_field(truth)
    for (r, c) in [(20, 22), (20, 30), (21, 22), (22, 25), (24, 25), (5, 5), (0, 0)]:
        cand = base.copy()
        cand[r, c] = 1.0
        s1 = M.score(cand, truth).dti
        predicted_helpful = wf[r, c] > bar
        actually_helpful = s1 > s0
        if abs(wf[r, c] - bar) > 1e-6:
            assert predicted_helpful == actually_helpful, (r, c, wf[r, c], bar, s0, s1)


def test_profitable_radius_values_quoted_in_docs():
    assert abs(M.marginal_bar(0.2600) - 0.052) < 1e-12
    assert abs(M.marginal_bar(0.2778) - 0.05556) < 1e-9
    assert abs(M.marginal_bar(0.3195) - 0.0639) < 1e-12
    assert abs(M.marginal_bar(0.3262) - 0.06524) < 1e-9
    assert abs(M.profitable_radius_m(0.2600) - 284.4) < 1e-9
    assert abs(M.profitable_radius_m(0.3262) - 280.428) < 1e-3


def test_coverage_inverts_recall():
    for s in (0.2, 0.2600, 0.3195, 0.3262):
        for rho in (0.0, 0.5, 1.0, 2.0):
            c = M.coverage_for_target(s, rho)
            # c must reproduce s when fed back through the design equation
            back = c / (M.ALPHA * (c + rho) + M.BETA)
            assert abs(back - s) < 1e-12, (s, rho, c, back)


def test_recall_for_target_matches_leader_arithmetic():
    """At DTI=0.3262: 27.9% weighted recall at zero false-positive mass."""
    assert abs(M.recall_for_target(0.3262, 0.0) - 0.2793) < 5e-4
