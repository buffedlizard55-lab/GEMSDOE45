"""Kaplan-Meier estimator: closed-form checks against hand-computable cases."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from gems45 import survival  # noqa: E402


def test_km_no_censoring_matches_ecdf():
    t = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    e = np.ones(5, dtype=int)
    km = survival.kaplan_meier(t, e)
    # S(t) after each event: 4/5, 3/5, 2/5, 1/5, 0
    assert km.survival == pytest.approx([0.8, 0.6, 0.4, 0.2, 0.0])
    # median = smallest t with S(t) <= 0.5; S drops to 0.4 at t = 3
    assert km.quantile(0.5) == pytest.approx(3.0)
    assert km.n_events == 5 and km.n_censored == 0


def test_km_censoring_raises_survival_above_naive():
    """The whole reason for using KM: censored (long) observations must not be dropped."""
    observed = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    censored = np.array([100.0, 100.0, 100.0, 100.0, 100.0])
    t = np.concatenate([observed, censored])
    e = np.concatenate([np.ones(5, dtype=int), np.zeros(5, dtype=int)])
    km = survival.kaplan_meier(t, e)
    assert km.n_censored == 5
    assert km.mean() > float(observed.mean())


def test_km_survival_is_non_increasing_and_starts_at_one():
    rng = np.random.default_rng(1)
    t = rng.integers(1, 40, 300).astype(float)
    e = (rng.random(300) < 0.8).astype(int)
    km = survival.kaplan_meier(t, e)
    assert km(np.array([0.0]))[0] == pytest.approx(1.0)
    assert np.all(np.diff(km.survival) <= 1e-12)
    assert (km.survival >= 0).all() and (km.survival <= 1).all()


def test_km_sample_respects_support():
    t = np.arange(1, 21, dtype=float)
    km = survival.kaplan_meier(t, np.ones(20, dtype=int))
    draws = km.sample(500, np.random.default_rng(2))
    assert draws.min() >= 1.0 and draws.max() <= t.max() + 1e-9


def test_per_tip_extension_varies():
    """The brief's explicit requirement: lengths must not collapse to a constant."""
    from gems45 import detector
    t = np.arange(1, 41, dtype=float)
    km = survival.kaplan_meier(t, np.ones(40, dtype=int))
    gaps = np.array([1, 2, 4, 6, 10, 20, 40, 61.0])
    events = np.array([1, 1, 1, 1, 1, 1, 1, 0])
    L = detector.per_tip_extension(gaps, events, km, e_cap=11)
    assert len(np.unique(L)) > 1
    assert L.min() >= 1 and L.max() <= 11


def test_expected_weight_profile_shapes_and_is_nonnegative():
    """For a uniform gap law the hazard is constant, so the profile is flat except at the very
    start where the kernel window is clipped at zero.  Sanity: right shape, no negatives."""
    from gems45 import detector
    t = np.arange(1, 121, dtype=float)
    km = survival.kaplan_meier(t, np.ones(120, dtype=int))
    prof = detector.expected_weight_profile(km, e_max=20)
    assert prof.shape == (20,)
    assert (prof >= 0).all()
    assert prof[0] < prof[5]
    assert prof[5] == pytest.approx(prof[-1])


def test_expected_weight_profile_decreases_for_a_decreasing_gap_density():
    """With an exponentially distributed gap the hazard rises, so the profile must fall."""
    from gems45 import detector
    rng = np.random.default_rng(9)
    t = np.minimum(rng.exponential(8.0, 5000), 120.0)
    e = (t < 120.0).astype(int)
    km = survival.kaplan_meier(t[t > 2], e[t > 2])
    prof = detector.expected_weight_profile(km, e_max=20)
    assert prof[2] > prof[-1]
    assert (prof >= 0).all()
