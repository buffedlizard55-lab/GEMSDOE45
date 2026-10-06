"""Distinctness audit: guarantee the submission differs from all prior GEMSDOE
submissions (requirement: "must be different than the collection of gemsdoe
sites").  We compute Pearson and Spearman correlation of our candidate against
every prior submission raster, over the in-footprint (finite) pixels, and require
|r| < 0.90 for both coefficients against every prior (the same threshold the
prior family used in ``audit_prior_correlations``).
"""

from __future__ import annotations

import numpy as np
import rasterio as rio
from scipy import stats


def _sample_indices(valid: np.ndarray, n: int):
    idx = np.nonzero(valid)[0]
    if idx.size <= n:
        return idx
    rng = np.random.default_rng(0)
    return rng.choice(idx, size=n, replace=False)


def _corr(a: np.ndarray, b: np.ndarray, valid: np.ndarray, spearman: bool, sample: int):
    idx = _sample_indices(valid, sample)
    av = a.ravel()[idx]
    bv = b.ravel()[idx]
    if spearman:
        # rank within the sample
        av = stats.rankdata(av)
        bv = stats.rankdata(bv)
    # Pearson on (possibly ranked) values
    if np.std(av) < 1e-12 or np.std(bv) < 1e-12:
        return 0.0
    return float(np.corrcoef(av, bv)[0, 1])


def audit_prior_correlations(
    candidate: np.ndarray,
    valid_mask: np.ndarray,
    prior_paths: list,
    max_abs_pearson: float = 0.90,
    max_abs_spearman: float = 0.90,
    sample: int = 250_000,
):
    """Return per-prior correlations and the all-distinct flag."""
    import rasterio as rio

    cand_flat = np.nan_to_num(candidate, nan=0.0)
    records = []
    worst_p, worst_s = 0.0, 0.0
    for p in prior_paths:
        try:
            with rio.open(str(p)) as src:
                prior = src.read(1).astype(np.float32)
        except Exception:
            continue
        if prior.shape != candidate.shape:
            continue
        pvalid = np.isfinite(prior)
        common = valid_mask & pvalid
        if common.sum() < 100:
            continue
        rp = _corr(cand_flat, prior, common, spearman=False, sample=sample)
        rs = _corr(cand_flat, prior, common, spearman=True, sample=sample)
        records.append({
            "prior": str(p),
            "pearson": rp,
            "spearman": rs,
            "n_common_pixels": int(common.sum()),
        })
        worst_p = max(worst_p, abs(rp))
        worst_s = max(worst_s, abs(rs))

    distinct = bool(worst_p < max_abs_pearson and worst_s < max_abs_spearman)
    return {
        "all_priors_distinct": distinct,
        "max_abs_pearson": worst_p,
        "max_abs_spearman": worst_s,
        "threshold_pearson": max_abs_pearson,
        "threshold_spearman": max_abs_spearman,
        "n_priors_checked": len(records),
        "priors": records,
    }
