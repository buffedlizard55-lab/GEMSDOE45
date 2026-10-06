#!/usr/bin/env python
"""Build the submission GeoTIFF from cached tips + a re-fitted Kaplan-Meier curve.

Fast to re-run (no thinning, no re-skeletonisation), so the emission rule can be iterated and
re-audited independently of the slow catalogue extraction.

Emission rule (exact form of the metric's first-order condition, metric.marginal_bar):
  * every tip contributes dots at along-strike distances e = 1..L_i, where L_i is that tip's own
    Kaplan-Meier extension length (see detector.per_tip_extension) capped at the distance where the
    KM-implied expected weight falls to alpha*DTI;
  * relay bridges are drawn between facing tips only when the measured gap is short enough that a
    bridge dot's expected weight still clears the same bar (G < 50 px);
  * dots landing on known-catalogue pixels are dropped: those pixels are masked out of scoring by
    the organizers, so they can neither earn credit nor incur a penalty.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gems45 import detector, emission, grid, metric, survival  # noqa: E402


def build(out_name: str, target_dti: float, bridge_max_gap: int, mask_catalogue: bool = True,
          zeros_twin: bool = True) -> dict:
    d = np.load(ROOT / "evidence" / "catalogue_cache.npz")
    known = grid.read_labels(ROOT / "data" / "labels.tif")
    tmpl = grid.read_template_mask(ROOT / "data" / "sample_submission.tif")

    rows, cols = d["rows"], d["cols"]
    ur, uc = d["ur"], d["uc"]
    gaps, events = d["gaps"], d["events"]

    class T:  # minimal duck-typed tip
        __slots__ = ("row", "col", "urow", "ucol")

        def __init__(self, r, c, a, b):
            self.row, self.col, self.urow, self.ucol = int(r), int(c), float(a), float(b)

    tips = [T(*t) for t in zip(rows, cols, ur, uc)]
    km_gen = survival.kaplan_meier(gaps[(events == 1) & (gaps > detector.G_FRAGMENT_MAX)],
                                   events[(events == 1) & (gaps > detector.G_FRAGMENT_MAX)])
    budget = detector.extension_budget(km_gen, target_dti)
    L = detector.per_tip_extension(gaps, events, km_gen, e_cap=budget["e_star_px"])

    tip_mask = emission.emit_tip_extensions(known.shape, tips, L)
    # Relay bridges, gated: a bridge dot's expected kernel weight under a uniform-within-gap
    # position model is ~0.5 * min(1, 2R/G), which clears alpha*DTI only while G < 50 px.
    bridge = np.zeros(known.shape, dtype=np.float32)
    for t, L_i, G, ev in zip(tips, L, gaps, events):
        if ev != 1 or G <= detector.G_FRAGMENT_MAX or G >= bridge_max_gap:
            continue
        e = int(np.floor(0.25 * G))
        while e <= int(np.ceil(0.75 * G)):
            r = int(round(t.row + e * t.urow))
            c = int(round(t.col + e * t.ucol))
            if 0 <= r < known.shape[0] and 0 <= c < known.shape[1]:
                bridge[r, c] = 1.0
            e += 1

    pred = np.maximum(tip_mask, bridge)
    n_before = int((pred > 0).sum())
    if mask_catalogue:
        pred[known] = 0.0
    pred = np.clip(pred, 0.0, 1.0).astype(np.float32)

    out = ROOT / "docs" / "downloads" / out_name
    grid.write_submission(pred, tmpl, out, outside="nan")
    a = grid.audit(out, tmpl)
    res = {
        "target_dti": target_dti, "bridge_max_gap_px": bridge_max_gap,
        "e_star_px": budget["e_star_px"], "marginal_bar": budget["marginal_bar"],
        "extension_lengths": {
            "min": int(L.min()), "max": int(L.max()), "mean": round(float(L.mean()), 3),
            "median": float(np.median(L)), "n_distinct": int(len(np.unique(L))),
            "histogram": {str(k): int((L == k).sum()) for k in range(int(L.max()) + 1)},
        },
        "tip_px": int((tip_mask > 0).sum()), "bridge_px": int((bridge > 0).sum()),
        "union_before_catalogue_mask_px": n_before,
        "n_positive": a["n_positive"],
        "mass": a["n_positive"],
        "audit": a,
        "zeros_twin": None,
    }
    res["km_genuine"] = km_gen.as_dict()
    if zeros_twin:
        z = np.where(tmpl, pred, np.float32(0.0))
        out_z = out.with_name(out.name.replace(".tif", "-zeros.tif"))
        grid.write_submission(z, np.ones_like(tmpl), out_z, outside="zero")
        res["zeros_twin"] = grid.audit(out_z, tmpl)
    return res


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", default="gems45-h46-kmtip-survival-20261006.tif")
    ap.add_argument("--target-dti", type=float, default=0.30)
    ap.add_argument("--bridge-max-gap", type=int, default=50)
    ap.add_argument("--out-json", default="evidence/submission_build.json")
    a = ap.parse_args()
    res = build(a.name, a.target_dti, a.bridge_max_gap)
    (ROOT / a.out_json).write_text(json.dumps(res, indent=1))
    print(json.dumps({k: v for k, v in res.items() if k != "audit"}, indent=1)[:3000])
    print("\nAUDIT:", json.dumps({k: res["audit"][k] for k in
          ("filename", "bytes", "sha256", "shape", "dtype", "crs", "inside_all_in_0_1",
           "outside_all_nan", "sentinel_anywhere", "n_positive")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
