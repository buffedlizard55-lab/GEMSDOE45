"""Build the H51 submission: Kaplan-Meier fault-zone emission at the metric's marginal bar.

Pipeline
--------
1. Thin the provided catalogue to a skeleton; estimate the local strike everywhere on it.
2. Fit the Kaplan-Meier product-limit curve to the along-strike distance from each of the
   catalogue's own tips to the nearest *other* mapped segment (the catalogue's relay-gap sample).
   That curve is, verbatim, the probability that a fault is still going e px past a mapped tip.
3. Draw a per-tip extension length by inverse transform from that curve, and enumerate the
   candidates the organizers' own definition of "new fault" names: along-strike continuations past
   tips, and splays / parallel strands at perpendicular offsets inside the same fault zone.
4. Gate every candidate: corroboration by a topographic-scarp or potential-field ridge response,
   at least `R_ZERO` px from the mapped catalogue (the live-verified flank exclusion), inside the
   template footprint, and never on a mapped pixel.
5. Select with greedy maximum expected coverage under the official 300 m triangular kernel -- the
   metric's own submodular objective -- stopping at its exact marginal bar evaluated at the
   portfolio's best live-observed DTI.
6. Write, re-read and audit both encodings; report the measured spread of extension lengths.

Run:  .venv/bin/python scripts/build_h51_submission.py
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
import zipfile
from pathlib import Path

import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.gems45.grid import (read_bands, read_labels, read_template_mask,  # noqa: E402
                             write_submission, audit, BAND_NAMES)
from src.gems45.catalogue import extract_tips  # noqa: E402
from src.gems45.survival import kaplan_meier, measure_gaps, S_MAX_PX  # noqa: E402
from src.gems45.faultzone import (skeleton, strike_field, tip_continuation_candidates,  # noqa: E402
                                  parallel_strand_candidates, evidence_field)
from src.gems45.h51 import novelty_profile  # noqa: E402
from src.gems45.coverage import greedy_coverage, own_weight  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
EV = ROOT / "evidence"
OUT = ROOT / "docs" / "downloads"
NAME = "gems45-h51-km-faultzone-20261006"
DTI_TARGET = 0.2778      # best live-observed DTI anywhere in this family (GEMSDOE32 H33-2-B2)
R_ZERO = 2.0             # catalogue-exclusion radius, calibrated on the family's live dose-response
K_MODEL = 12000.0        # modelled hidden-truth pixel count (independent estimates 8,936 / 12,632)
MAX_DOTS = 37654         # the largest dot count ever observed to be optimal in this family
SEED = 4506
MAX_EXTEND = 30
LENGTH_GAIN = 0.35
STRAND_OFFSETS = (-3, -2, 2, 3)


def _tips(dataset: np.ndarray) -> list:
    """Tip extraction is the single slowest step (~10 min); cache it next to the other evidence."""
    cache = EV / "h51_tips.npy"
    from src.gems45.catalogue import Tip
    if cache.exists():
        try:
            arr = np.load(cache, allow_pickle=True)
            if arr.dtype == object and hasattr(arr.flat[0], "urow"):
                return [t for t in arr]                      # already Tip objects
            arr = np.atleast_2d(np.asarray(arr, dtype=float))
            return [Tip(int(r[0]), int(r[1]), float(r[2]), float(r[3]), int(r[4]),
                        int(r[5]) if arr.shape[1] > 5 else -1) for r in arr]
        except Exception as exc:                              # noqa: BLE001 - rebuild below
            print(f"tip cache unreadable ({exc}); re-extracting (~10 min)")
    tips = extract_tips(dataset)
    np.save(cache, np.array([[t.row, t.col, t.urow, t.ucol, t.seg_len, t.nbr] for t in tips],
                            dtype=float))
    return tips


def draw_lengths(tips, survival, rng, length_gain=LENGTH_GAIN, cap=MAX_EXTEND):
    """Inverse-transform draws of the along-strike extension length from the fitted KM curve."""
    u = rng.random(len(tips))
    draws = np.array([float(survival.quantile(1.0 - ui)) for ui in u])
    lat = 1.0 + length_gain * np.log1p(np.array([t.seg_len for t in tips], dtype=float) / 12.0)
    return draws, lat, np.minimum(np.round(draws * lat), cap).astype(int)


def main() -> int:
    t0 = time.time()
    stack, _ = read_bands(DATA / "training_features.tif")
    catalog = read_labels(DATA / "labels.tif")
    tmask = read_template_mask(DATA / "sample_submission.tif")
    footprint = tmask & ~catalog
    shape = catalog.shape

    lab, _ = ndimage.label(catalog, structure=np.ones((3, 3), dtype=np.uint8))
    tips = _tips(catalog)
    gaps, events, _ = measure_gaps(tips, catalog, lab)
    km = kaplan_meier(gaps, events, s_max=float(S_MAX_PX))
    print(f"tips {len(tips)}  KM events {km.n_events} censored {km.n_censored} "
          f"median gap {km.quantile(0.5):.0f} px  restricted mean {km.mean():.1f} px")

    rng = np.random.default_rng(SEED)
    draws, lat, L = draw_lengths(tips, km, rng)
    print("per-tip drawn extension: mean %.2f px  median %.1f  distinct %d  nonzero %d"
          % (draws.mean(), np.median(draws), len(np.unique(draws)), int((draws > 0).sum())))

    # ---- candidates ---------------------------------------------------------------------
    tr_r, tr_c, tr_w, tr_id = tip_continuation_candidates(shape, tips, km, MAX_EXTEND)

    skel = skeleton(catalog)
    s_tr, s_tc = strike_field(catalog, sigma=3.0)
    st_r, st_c, st_w = parallel_strand_candidates(skel, s_tr, s_tc, offsets=STRAND_OFFSETS)
    print(f"skeleton px {int(skel.sum())}  tip candidates {len(tr_r)}  strand candidates {len(st_r)}")

    rr = np.concatenate([tr_r, st_r])
    cc = np.concatenate([tr_c, st_c])
    ww = np.concatenate([tr_w, st_w])
    src = np.concatenate([np.zeros(len(tr_r), np.int8), np.ones(len(st_r), np.int8)])

    # soft multi-physics gate: the candidate must sit on *something* the physics sees
    gate = evidence_field(stack, tmask, BAND_NAMES.index("det_elev_slope"), (1.2, 2.5), 99.0)
    grav = evidence_field(stack, tmask, BAND_NAMES.index("iso_grav_anom_hg"), (2.5,), 99.0)
    ww = ww * (0.25 + 0.75 * np.maximum(gate[rr, cc], grav[rr, cc]))

    # hard catalogue-novelty gate at the live-calibrated radius
    d_cat = ndimage.distance_transform_edt(~catalog).astype(np.float32)
    eta = novelty_profile(d_cat, R_ZERO, 2.0)
    ww = ww * eta[rr, cc].astype(np.float64)
    ok = footprint[rr, cc] & (ww > 0)
    rr, cc, ww, src = rr[ok], cc[ok], ww[ok], src[ok]
    print(f"candidates after gates {len(rr)}")

    weights = np.zeros(shape, dtype=np.float64)
    np.add.at(weights, (rr, cc), ww)
    cand_mask = weights > 0
    psi = (weights / weights.sum() * K_MODEL).astype(np.float64)
    print(f"candidate cells {int(cand_mask.sum())}")

    # the metric's own marginal bar, applied to the expected weight of each dot
    own = own_weight(psi, rr, cc)
    bar = 0.2 * DTI_TARGET
    print(f"own(x): median {np.median(own):.4f}  max {own.max():.4f}  bar {bar:.4f}  "
          f"above bar {int((own >= bar).sum())} of {len(own)}")
    sel_ok = own >= bar
    c2 = np.zeros(shape, dtype=bool)
    c2[rr[sel_ok], cc[sel_ok]] = True

    res = greedy_coverage(psi, c2, K=float(psi.sum()), max_dots=MAX_DOTS,
                          round_size=1000, dti_operating=DTI_TARGET)
    print("greedy:", json.dumps(res.as_dict())[:300])
    assert res.f_pred >= 0.0, "emulator F must be a non-negative mass"
    assert res.t_pred <= res.k_pred + 1e-6, "emulator T cannot exceed the modelled truth size"
    dots = np.zeros(shape, dtype=bool)
    dots[res.rows, res.cols] = True
    sel_src = np.zeros(shape, dtype=np.int8)
    sel_src[rr[sel_ok], cc[sel_ok]] = (src + 1)[sel_ok]   # 1 = tip continuation, 2 = strand
    emitted = sel_src[res.rows, res.cols]
    comp = dict(tip_continuation=int((emitted == 1).sum()),
                strand_offset=int((emitted == 2).sum()))
    print(f"emitted {int(dots.sum())} dots in {time.time() - t0:.0f}s  {comp}")

    # ---- realised extension per fault ---------------------------------------------------
    realised = np.zeros(len(tips), dtype=int)
    first_hit = np.zeros(len(tips), dtype=int)
    for i, t in enumerate(tips):
        first = -1
        for e in range(1, MAX_EXTEND + 1):
            r = int(round(t.row + e * t.urow))
            c = int(round(t.col + e * t.ucol))
            if not (0 <= r < shape[0] and 0 <= c < shape[1]):
                break
            if dots[r, c]:
                realised[i] = e
                if first < 0:
                    first = e
        first_hit[i] = first

    def hist(a):
        v, n = np.unique(a, return_counts=True)
        return {str(int(k)): int(x) for k, x in zip(v, n)}

    ext = dict(n_tips=len(tips),
               drawn_mean=round(float(draws.mean()), 3),
               drawn_median=float(np.median(draws)),
               drawn_distinct=int(len(np.unique(np.round(draws, 3)))),
               drawn_histogram=hist(L),
               realised_mean=round(float(realised.mean()), 3),
               realised_median=float(np.median(realised)),
               realised_distinct=int(len(np.unique(realised))),
               realised_histogram=hist(realised),
               n_tips_with_realised_extension=int((realised > 0).sum()),
               first_hit_median=float(np.median(first_hit[first_hit > 0])) if (first_hit > 0).any() else 0.0)
    print("drawn   :", json.dumps({k: ext[k] for k in ("drawn_mean", "drawn_median", "drawn_distinct")}))
    print("realised:", json.dumps({k: ext[k] for k in
          ("realised_mean", "realised_median", "realised_distinct", "n_tips_with_realised_extension")}))
    assert ext["realised_distinct"] >= 5, "extension lengths collapsed to a constant"
    assert ext["drawn_distinct"] >= 50, "per-tip draws collapsed"

    # ---- corridor coverage diagnostic ---------------------------------------------------
    from src.gems45.coverage import OFFS
    vals = np.zeros(shape, dtype=bool)
    for dy, dx, _ in OFFS:
        dy, dx = int(dy), int(dx)
        y0, y1 = max(0, -dy), min(shape[0], shape[0] - dy)
        x0, x1 = max(0, -dx), min(shape[1], shape[1] - dx)
        if y0 < y1 and x0 < x1:
            vals[y0:y1, x0:x1] |= dots[y0 + dy:y1 + dy, x0 + dx:x1 + dx]
    cover = float((vals & cand_mask).sum()) / max(1, int(cand_mask.sum()))

    # ---- write ---------------------------------------------------------------------------
    pred = dots.astype(np.float32)
    f_nan = OUT / f"{NAME}.tif"
    f_z = OUT / f"{NAME}-zeros.tif"
    write_submission(pred, tmask, f_nan, outside="nan")
    write_submission(pred, tmask, f_z, outside="zero")
    a_nan = audit(f_nan, tmask)
    a_z = audit(f_z, tmask)
    zpath = OUT / f"{NAME}-zeros.zip"
    with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
        z.write(f_z, f_z.name)
    for a, label in ((a_nan, "NaN-outside"), (a_z, "all-finite")):
        print(label, json.dumps({k: a[k] for k in
              ("bytes", "sha256", "inside_min", "inside_max", "nan_inside", "sentinel_anywhere",
               "n_positive")}))

    report = dict(
        name=NAME, artifact="H51 Kaplan-Meier fault-zone emission",
        created="2026-10-06", seed=SEED, status="UNSCORED CANDIDATE", slot_approved=False,
        emulator_is_not_a_score=True, score_claim=None,
        km=dict(n_events=int(km.n_events), n_censored=int(km.n_censored),
                median_gap_px=km.quantile(0.5), restricted_mean_gap_px=round(km.mean(), 3),
                survival_at={str(d): round(float(km(np.array([d]))[0]), 5)
                             for d in (1, 2, 3, 5, 8, 12, 20, 30)}),
        candidates=dict(tips=int(len(tips)), skeleton_px=int(skel.sum()),
                        tip_candidates=int(len(tr_r)), strand_candidates=int(len(st_r)),
                        after_gates=int(len(rr)), candidate_cells=int(cand_mask.sum())),
        rules=dict(marginal_bar=bar, dti_target=DTI_TARGET, dti_operating=DTI_TARGET, r_zero_px=R_ZERO,
                   k_model=K_MODEL, max_dots=MAX_DOTS, length_gain=LENGTH_GAIN,
                   strand_offsets_px=list(STRAND_OFFSETS),
                   stopping="greedy max-coverage; keep while expected marginal gain clears "
                            "dT*(1-0.2*DTI_op) > 0.2*DTI_op*(1-own)"),
        greedy=res.as_dict(), n_dots=int(dots.sum()), corridor_covered_fraction=round(cover, 4),
        emission_composition=comp,
        live_anchored_ceiling_dots=MAX_DOTS,
        live_evidence=("Within this family the public leaderboard shows a monotone dose-response as "
                       "mass is removed: 44,090 px -> 0.2600, 40,199 px -> 0.2708, 37,654 px -> "
                       "0.2778 (research/live_score_features.csv). Those are owner-reported scores "
                       "with no organizer receipt."),
        extension_lengths=ext,
        model_caveat=("T/F/DTI in `greedy` are outputs of the expected-coverage emulator under an "
                      "UNVALIDATED candidate density field. They are NOT score predictions: the "
                      "field was measured to have no usable rank correlation with live scores "
                      "(Pearson +0.176 over 46 live-scored owner rasters, evidence/h51_live_fit.json)."),
        files=dict(nan_outside=a_nan, all_finite=a_z,
                   zip=dict(name=zpath.name, bytes=zpath.stat().st_size,
                            sha256=hashlib.sha256(zpath.read_bytes()).hexdigest())),
        provenance=dict(
            catalogue="data/labels.tif sha256 7ba308ccdc4418b31a178f4f1ef21aaa6e152e4028f2f6f64b01f7eb25ae4093",
            features="data/training_features.tif sha256 4371c82e3b8339b807bdffcf4ef59a225520fe2988d521be208ae33743123bc5",
            template="data/sample_submission.tif sha256 2176d08e485aa2cd2860ce8df539db4faf4d76163b38a4dd8c30a40454d35cbc"),
        runtime_seconds=round(time.time() - t0, 1),
    )
    (EV / "h51_submission.json").write_text(json.dumps(report, indent=1, default=str) + "\n")
    print("wrote evidence/h51_submission.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
