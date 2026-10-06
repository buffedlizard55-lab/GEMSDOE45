"""Build the H51 expected-truth field psi and the Kaplan-Meier tip-continuation hazard.

psi(x) = normalised [ (multi-physics line corroboration) * (catalogue novelty)
                      + lambda_tip * (Kaplan-Meier tip continuation hazard) ]

Everything is built from the supplied 19-band stack and the provided catalogue only.  Expensive
stages are cached under evidence/ (git-ignored) so that re-runs after a parameter change cost
seconds rather than a quarter of an hour.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.gems45.grid import read_bands, read_labels, read_template_mask  # noqa: E402
from src.gems45.h51 import corroboration, novelty_profile, tip_continuation_field  # noqa: E402
from src.gems45.survival import kaplan_meier, measure_gaps, S_MAX_PX  # noqa: E402
from src.gems45.catalogue import extract_tips, catalogue_stats  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
EV = ROOT / "evidence"
EV.mkdir(exist_ok=True)


def cached(name, fn):
    p = EV / name
    if p.exists():
        return np.load(p, allow_pickle=True)
    t0 = time.time()
    v = fn()
    np.save(p, np.asarray(v, dtype=object), allow_pickle=True)
    print(f"  computed {name} in {time.time() - t0:.1f}s")
    return np.asarray(v, dtype=object)


def main() -> int:
    catalog = read_labels(DATA / "labels.tif")
    mask = read_template_mask(DATA / "sample_submission.tif")
    footprint = mask

    # ---- Kaplan-Meier tip survival from the catalogue's own relay gaps -------------------
    lab, ncomp = ndimage.label(catalog, structure=np.ones((3, 3), dtype=np.uint8))
    tips = cached("h51_tips.npy", lambda: np.asarray(
        [(t.row, t.col, t.urow, t.ucol, t.seg_len) for t in extract_tips(catalog)], dtype=np.float64))
    from src.gems45.catalogue import Tip
    tiplist = [Tip(int(t[0]), int(t[1]), float(t[2]), float(t[3]), int(t[4]), -1) for t in tips]
    gaps_ev = cached("h51_gaps.npy", lambda: np.asarray(
        measure_gaps(tiplist, catalog, lab)[:2], dtype=np.float64))
    gaps = gaps_ev[0].astype(float)
    events = gaps_ev[1].astype(int)
    km = kaplan_meier(gaps, events, s_max=float(S_MAX_PX))
    print("KM events %d censored %d median %.1f px" % (km.n_events, km.n_censored, km.quantile(0.5)))

    tip_field = cached("h51_tiphazard.npy",
                       lambda: tip_continuation_field(catalog.shape, tiplist, km, max_extend=30))
    tip_field = np.asarray(tip_field, dtype=np.float32)
    print("tip hazard nonzero %d max %.3f" % (int((tip_field > 0).sum()), float(tip_field.max())))

    # ---- corroboration ---------------------------------------------------------------------
    def _corr():
        stack, _ = read_bands(DATA / "training_features.tif")
        c, t, q = corroboration(stack, footprint)
        for nm, arr in (("h51_corr.npy", c), ("h51_topo.npy", t), ("h51_pot.npy", q)):
            np.save(EV / nm, np.asarray(arr, dtype=np.float32))
        return np.zeros(1)
    cached("h51_corr_done.npy", _corr)
    corr = np.load(EV / "h51_corr.npy").astype(np.float32)
    topo = np.load(EV / "h51_topo.npy").astype(np.float32)
    pot = np.load(EV / "h51_pot.npy").astype(np.float32)
    print("corroboration max %.4f mean %.5f" % (float(corr.max()), float(corr[footprint].mean())))

    d_cat = cached("h51_dcat.npy", lambda: ndimage.distance_transform_edt(~catalog).astype(np.float32))
    d_cat = np.asarray(d_cat, dtype=np.float32)

    (EV / "h51_km.json").write_text(json.dumps(
        dict(tips=len(tiplist), events=int(km.n_events), censored=int(km.n_censored),
             curve=km.as_dict(), catalogue=catalogue_stats(catalog), components=int(ncomp)), indent=1) + "\n")
    print("wrote evidence/h51_km.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
