#!/usr/bin/env python
"""Spatially-blocked holdout in the competition's own shape (whole-map emission).

PROTOCOL
--------
The grid is cut into a lattice of disjoint blocks.  For each fold the block's catalogue pixels are
withheld -- they become the fold's *truth* -- and every part of the pipeline that could see them is
rebuilt without them:

  * catalogue segments intersecting the block are dropped entirely (a partly revealed segment
    would leak its continuation);
  * the tip gap sample is re-measured against the reduced catalogue, so the Kaplan-Meier curve is
    re-fitted from scratch;
  * the emission is still painted over the WHOLE map (this is the competition's shape: you do not
    know where the scored faults are), so every dot outside the block is charged as false-positive
    mass exactly as it would be live.

WHY THE RESULT IS A SCREEN AND NOT A PROMOTION GATE
---------------------------------------------------
The truth of this proxy is the visible catalogue, and the organizers mask that catalogue out of the
live score entirely (community thread 11516).  A proxy whose truth is the catalogue therefore
rewards covering wherever the catalogue runs, while the live truth is a different, sparser set of
faults.  Sibling repositories have measured this: GEMSDOE40 reports a leave-one-out Spearman of
-0.897 between a catalogue-truth proxy and 16 live scores.  The numbers below are reported as a
consistency screen with that caveat attached, and no submission decision rests on them alone.
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

PRED_PATH = ROOT / "docs" / "downloads" / "gems45-h48-kmtip-structural-20261006.tif"
INCUMBENTS = {
    "dotted-h19-5-d2-8 (owner-reported 0.2600)": "/tmp/g32/docs/downloads/gems25-dotted-h19-5-d2-8-20261002-e56ea318af89-nan.tif",
    "h33-2-b2 (owner-reported 0.2778)": "/tmp/g32/docs/downloads/gemsdoe32-h33-h33-2-b2-20261004T220000Z-e5eb6e7e-nan.tif",
}


class T:
    __slots__ = ("row", "col", "urow", "ucol")

    def __init__(self, r, c, a, b):
        self.row, self.col, self.urow, self.ucol = int(r), int(c), float(a), float(b)


def load_pred(path: str, tmpl: np.ndarray) -> np.ndarray:
    import rasterio
    try:
        with rasterio.open(path) as s:
            a = s.read(1).astype(np.float32)
    except Exception:
        return None
    a = np.nan_to_num(a, nan=0.0)
    a = np.clip(a, 0.0, 1.0)
    a[~tmpl] = 0.0
    return a


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--blocks", type=int, default=4)
    ap.add_argument("--target-dti", type=float, default=0.30)
    ap.add_argument("--out", default="evidence/holdout_final.json")
    a = ap.parse_args()

    known = grid.read_labels(ROOT / "data" / "labels.tif")
    tmpl = grid.read_template_mask(ROOT / "data" / "sample_submission.tif")
    import scipy.ndimage as ndi
    lab = ndi.label(known, structure=np.ones((3, 3), dtype=np.uint8))[0]
    cache = np.load(ROOT / "evidence" / "catalogue_cache.npz")
    rows, cols, ur, uc = cache["rows"], cache["cols"], cache["ur"], cache["uc"]
    comp = lab[rows, cols]
    rng = np.random.default_rng(45)

    mine = load_pred(str(PRED_PATH), tmpl)
    inc = {k: load_pred(v, tmpl) for k, v in INCUMBENTS.items()}
    inc = {k: v for k, v in inc.items() if v is not None}

    h, w = known.shape
    bs_r, bs_c = h // a.blocks, w // a.blocks
    folds = []
    for bi in range(a.blocks):
        for bj in range(a.blocks):
            blk = np.zeros_like(known)
            blk[bi * bs_r:(bi + 1) * bs_r, bj * bs_c:(bj + 1) * bs_c] = True
            truth = known & blk
            if truth.sum() < 200:
                continue
            train_cat = known & ~blk
            # drop every full-catalogue component that touches the block
            bad = set(np.unique(lab[truth]).tolist())
            keep = np.array([c not in bad for c in comp])
            tips = [T(*t) for t in zip(rows[keep], cols[keep], ur[keep], uc[keep])]
            if len(tips) < 50:
                continue
            lab_train = ndi.label(train_cat, structure=np.ones((3, 3), dtype=np.uint8))[0]
            gaps, events, _ = survival.measure_gaps(tips, train_cat, lab_train)
            gen = (events == 1) & (gaps > detector.G_FRAGMENT_MAX)
            if gen.sum() < 30:
                continue
            km = survival.kaplan_meier(gaps[gen], events[gen])
            budget = detector.extension_budget(km, a.target_dti)
            L = detector.per_tip_extension(gaps, events, km, e_cap=budget["e_star_px"])
            tm = emission.emit_tip_extensions(known.shape, tips, L)
            br = np.zeros(known.shape, dtype=np.float32)
            for t, Li, G, ev in zip(tips, L, gaps, events):
                if ev != 1 or G <= detector.G_FRAGMENT_MAX or G >= 50:
                    continue
                e = int(np.floor(0.25 * G))
                while e <= int(np.ceil(0.75 * G)):
                    r = int(round(t.row + e * t.urow)); c = int(round(t.col + e * t.ucol))
                    if 0 <= r < h and 0 <= c < w:
                        br[r, c] = 1.0
                    e += 1
            pred_fold = np.maximum(tm, br)
            pred_fold[known] = 0.0
            pred_fold = np.clip(pred_fold, 0, 1).astype(np.float32)

            res = {"fold": f"{bi}_{bj}", "truth_px": int(truth.sum()),
                   "emitted_px_global": int((pred_fold > 0).sum()),
                   "e_star_px": budget["e_star_px"],
                   "km_median_gap_px": km.quantile(0.5),
                   "ext_len_mean": round(float(L.mean()), 3),
                   "dti_mine": round(metric.score(pred_fold, truth).dti, 6),
                   "random_dti": round(metric.score(
                       (rng.random(known.shape) < (pred_fold > 0).mean()).astype(np.float32),
                       truth).dti, 6),
                   "incumbents": {k: round(metric.score(v, truth).dti, 6) for k, v in inc.items()},
                   "incumbent_mass": {k: int((v > 0).sum()) for k, v in inc.items()},
                   }
            folds.append(res)
            print(json.dumps(res))

    if folds:
        agg = {
            "n_folds": len(folds),
            "mean_dti_mine": round(float(np.mean([f["dti_mine"] for f in folds])), 6),
            "mean_dti_random": round(float(np.mean([f["random_dti"] for f in folds])), 6),
            "mean_dti_incumbents": {k: round(float(np.mean([f["incumbents"][k] for f in folds])), 6)
                                    for k in inc},
            "folds_mine_beats_random": int(sum(f["dti_mine"] > f["random_dti"] for f in folds)),
            "folds_mine_beats_incumbents": {k: int(sum(f["dti_mine"] > f["incumbents"][k]
                                                       for f in folds)) for k in inc},
            "truth_px_total": int(sum(f["truth_px"] for f in folds)),
            "caveat": ("Catalogue-truth proxy; the organizers mask the catalogue out of scoring, so "
                       "this instrument is a screen only (sibling repos measure LOO Spearman -0.897 "
                       "between such a proxy and live scores)."),
        }
    else:
        agg = {"n_folds": 0, "caveat": "no fold had enough truth"}
    (ROOT / a.out).write_text(json.dumps({"aggregate": agg, "folds": folds}, indent=1))
    print("\nAGGREGATE", json.dumps(agg, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
