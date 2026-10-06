#!/usr/bin/env python
"""Spatially-blocked supervised lineament detector.

Protocol
--------
The grid is cut into B x B blocks.  For each fold the model is trained on every block except the
held-out one and scored on the held-out one, so no pixel's neighbours in the same block can leak
into its own prediction.  Reported per fold: ROC AUC, and the *enrichment* that matters for the
competition -- the fraction of the top-N ranked pixels inside blocks that lie within the metric's
300 m kernel of a real fault.

Features deliberately exclude any distance-to-catalogue term.  Including one would let the model
predict "where the catalogue runs", which the organizers mask out of scoring, and would produce a
field whose apparent skill is an artefact of the proxy.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gems45 import grid, detfeatures  # noqa: E402
from gems45.grid import BAND_NAMES  # noqa: E402

GRAD_BANDS = ["det_elev_slope", "iso_grav_anom_hg", "iso_grav_anom_slope", "iso_grav_anom",
              "geod_2ndinv", "cond_surf", "tmi_hg", "rtp", "tmi_vg", "det_elev", "tc",
              "iso_grav_anom_vg", "depth_to_base_surf", "geod_dilaterate"]
VAL_BANDS = ["det_elev_slope", "cond_surf", "depth_to_base_surf", "iso_grav_anom", "det_elev",
             "geod_2ndinv", "tc"]
SIGMAS = (1.2, 2.5)


def build_features(stack: np.ndarray, footprint: np.ndarray) -> tuple[np.ndarray, list[str]]:
    idx = {n: i for i, n in enumerate(BAND_NAMES)}
    feats, names = [], []
    for s in SIGMAS:
        for n in GRAD_BANDS:
            b = stack[idx[n]].astype(np.float64)
            v = b[footprint]
            z = ndimage.gaussian_filter((b - v.mean()) / (v.std() + 1e-12), s)
            g = np.hypot(*np.gradient(z))
            feats.append(g.astype(np.float32))
            names.append(f"|grad{s}|_{n}")
    for n in VAL_BANDS:
        b = stack[idx[n]].astype(np.float64)
        v = b[footprint]
        feats.append(((b - v.mean()) / (v.std() + 1e-12)).astype(np.float32))
        names.append(f"val_{n}")
    # convergence/divergence of the detrended-elevation slope field: textural context
    ds = stack[idx["det_elev_slope"]].astype(np.float64)
    ds = ndimage.gaussian_filter((ds - ds[footprint].mean()) / (ds[footprint].std() + 1e-12), 1.2)
    lap = ndimage.laplace(ds)
    feats.append(np.abs(lap).astype(np.float32)); names.append("|laplace|_det_elev_slope")
    feats.append(lap.astype(np.float32)); names.append("laplace_det_elev_slope")
    return np.stack(feats), names


def fit_logreg(X: np.ndarray, y: np.ndarray, l2: float = 1.0, iters: int = 220) -> np.ndarray:
    from scipy.optimize import minimize
    Xb = np.hstack([X, np.ones((len(X), 1), dtype=X.dtype)])
    d = Xb.shape[1]

    def nll(w):
        z = Xb @ w
        # stable log-loss
        ll = np.sum(np.where(y > 0, -np.logaddexp(0, -z), -np.logaddexp(0, z)))
        ll += 0.5 * l2 * np.sum(w[:-1] ** 2)
        p = 1.0 / (1.0 + np.exp(-z))
        g = Xb.T @ (p - y)
        g[:-1] += l2 * w[:-1]
        return ll / len(Xb), g / len(Xb)

    w0 = np.zeros(d)
    res = minimize(nll, w0, jac=True, method="L-BFGS-B",
                   options={"maxiter": iters, "ftol": 1e-10})
    return res.x


def auc(y: np.ndarray, s: np.ndarray) -> float:
    n1 = int(y.sum()); n0 = int(len(y) - n1)
    if n1 == 0 or n0 == 0:
        return float("nan")
    order = np.argsort(s, kind="stable")
    ranks = np.empty(len(s)); ranks[order] = np.arange(len(s))
    return float((ranks[y > 0].sum() - n1 * (n1 - 1) / 2) / (n1 * n0))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--blocks", type=int, default=4)
    ap.add_argument("--sample", type=int, default=400_000)
    ap.add_argument("--topn", type=int, default=40000)
    ap.add_argument("--out", default="evidence/detector.json")
    a = ap.parse_args()
    t0 = time.time()

    stack, foot = grid.read_bands(ROOT / "data" / "training_features.tif")
    known = grid.read_labels(ROOT / "data" / "labels.tif")
    valid = foot.copy()
    flat_ok = np.nonzero(valid.ravel())[0]
    rng = np.random.default_rng(45)
    idx = np.sort(rng.choice(flat_ok, size=min(a.sample, len(flat_ok)), replace=False))
    y = known.ravel()[idx].astype(np.float64)
    X = detfeatures.sample_matrix(stack, foot, idx).astype(np.float64)
    names = detfeatures.feature_names()
    print("features", X.shape, "pos", int(y.sum()), "in", round(time.time() - t0, 1), "s")

    h, w = known.shape
    br, bc = h // a.blocks, w // a.blocks
    r = idx // w
    c = idx % w
    blk = (r // br) * a.blocks + (c // bc)
    folds, oof = [], np.full(len(idx), np.nan)
    for k in range(a.blocks * a.blocks):
        te = blk == k
        tr = ~te
        if te.sum() < 500 or y[tr].sum() < 50 or y[te].sum() < 20:
            continue
        wgt = fit_logreg(X[tr], y[tr])
        s = X[te] @ wgt[:-1] + wgt[-1]
        oof[te] = 1.0 / (1.0 + np.exp(-s))
        folds.append({"block": int(k), "n_test": int(te.sum()), "pos_test": int(y[te].sum()),
                      "auc": round(auc(y[te], s), 5)})
        print(json.dumps(folds[-1]))
    pooled = auc(y, oof)
    print("pooled OOF AUC", round(pooled, 5), "| mean fold AUC",
          round(float(np.mean([f["auc"] for f in folds])), 5))

    # enrichment of the top-N within held-out blocks only (honest budget test)
    enr = {}
    for k in range(a.blocks * a.blocks):
        te = blk == k
        if te.sum() < 500 or y[te].sum() < 20:
            continue
        budget = int(a.topn * te.mean())
        s = oof[te]
        order = np.argsort(-s)[:budget]
        ii = idx[te][order]
        rr, cc = ii // w, ii % w
        d = ndimage.distance_transform_edt(~known)
        dd = d[rr, cc]
        enr[int(k)] = {"budget": budget, "frac_within_300m": round(float((dd <= 3).mean()), 4),
                       "median_px": round(float(np.median(dd)), 2)}
    base = float((ndimage.distance_transform_edt(~known).ravel()[idx] <= 3).mean())
    out = {"features": names, "folds": folds, "pooled_oof_auc": round(pooled, 5),
           "enrichment_top_n": enr,
           "mean_frac_within_300m": round(float(np.mean([v["frac_within_300m"] for v in enr.values()])), 4),
           "footprint_base_rate_within_300m": round(base, 4),
           "n_folds": len(folds), "seconds": round(time.time() - t0, 1)}

    # final model on everything (for emission, reported with its own OOF AUC as the honest skill)
    w_full = fit_logreg(X, y)
    np.savez_compressed(ROOT / "evidence" / "detector_model.npz", w=w_full, names=np.array(names))
    (ROOT / a.out).write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1)[:2500])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
