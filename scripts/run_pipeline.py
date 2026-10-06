#!/usr/bin/env python
"""End-to-end: data -> band screen -> catalogue tips -> Kaplan-Meier -> emission -> audit.

Writes every intermediate number into ``evidence/`` so that each claim in the site can be traced
back to a file produced by this script on hash-pinned inputs.
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gems45 import bands, catalogue, detector, emission, grid, metric, survival  # noqa: E402

DATA = ROOT / "data"
EV = ROOT / "evidence"
EV.mkdir(exist_ok=True)


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    t_start = time.time()
    log: dict = {"started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}

    # ---------------------------------------------------------------- inputs
    log["inputs"] = {p.name: {"bytes": p.stat().st_size, "sha256": sha(p)}
                     for p in sorted(DATA.glob("*.tif"))}
    stack, footprint = grid.read_bands(DATA / "training_features.tif")
    known = grid.read_labels(DATA / "labels.tif")
    tmpl = grid.read_template_mask(DATA / "sample_submission.tif")
    log["geometry"] = {
        "shape": list(stack.shape[1:]),
        "n_bands": int(stack.shape[0]),
        "footprint_band1_px": int(footprint.sum()),
        "footprint_template_px": int(tmpl.sum()),
        "footprint_in_all_bands_px": int(footprint.sum()),
        "known_fault_px": int(known.sum()),
        "irregularity": ("band-1 finite mask (5,165,852), all-band finite mask (5,165,840) and the "
                         "sample submission's finite mask (5,167,373) differ; the writer uses the "
                         "sample submission mask because the problem description designates it as "
                         "the format template.  Band 6 `tc` carries 12 sentinel pixels inside the "
                         "band-1 mask.  See IR-45-001."),
    }

    # ---------------------------------------------------------------- band screen
    t0 = time.time()
    info = bands.band_information(stack, footprint, known)
    log["band_information"] = info
    log["band_screen_seconds"] = round(time.time() - t0, 1)

    # ---------------------------------------------------------------- catalogue + tips
    t0 = time.time()
    stats = catalogue.catalogue_stats(known)
    tips = catalogue.extract_tips(known)
    log["catalogue"] = stats
    log["n_tips"] = len(tips)
    log["tip_extraction_seconds"] = round(time.time() - t0, 1)

    # ---------------------------------------------------------------- Kaplan-Meier
    t0 = time.time()
    lab = __import__("scipy.ndimage", fromlist=["label"]).label(
        known, structure=np.ones((3, 3), dtype=np.uint8))[0]
    gaps, events, _ = survival.measure_gaps(tips, known, lab)
    km = survival.kaplan_meier(gaps, events)
    log["kaplan_meier"] = km.as_dict()
    gen = (events == 1) & (gaps > detector.G_FRAGMENT_MAX)
    km_genuine = survival.kaplan_meier(gaps[gen], events[gen])
    log["kaplan_meier_genuine"] = km_genuine.as_dict()
    log["gap_event_histogram_px"] = {
        str(int(b)): int(((gaps >= b) & (gaps < b + 5) & (events == 1)).sum()) for b in range(0, 65, 5)
    }
    log["km_seconds"] = round(time.time() - t0, 1)

    # ---------------------------------------------------------------- budget from the metric
    budget = detector.extension_budget(km_genuine, detector.TARGET_DTI)
    log["extension_budget"] = budget

    # ---------------------------------------------------------------- per-tip extension
    L = detector.per_tip_extension(gaps, events, km_genuine, e_cap=budget["e_star_px"])
    log["extension_lengths_px"] = {
        "n_tips": int(len(L)),
        "min": int(L.min()), "max": int(L.max()), "mean": round(float(L.mean()), 3),
        "median": float(np.median(L)),
        "histogram": {str(k): int((L == k).sum()) for k in range(0, int(L.max()) + 1)},
        "n_distinct_values": int(len(np.unique(L))),
        "is_constant": bool(len(np.unique(L)) == 1),
        "total_m": float(L.sum() * 100.0),
    }
    assert len(np.unique(L)) > 1, "extension lengths collapsed to a constant"

    # ---------------------------------------------------------------- emission
    t0 = time.time()
    strike = detector.strike_agreement_field(stack, footprint)
    tip_mask = emission.emit_tip_extensions(known.shape, tips, L)
    bridge = emission.emit_step_over_bridges(known.shape, tips, gaps, events, L)
    raw = np.maximum(tip_mask, bridge)
    log["emission"] = {
        "tip_px": int((tip_mask > 0).sum()),
        "bridge_px": int((bridge > 0).sum()),
        "raw_union_px": int((raw > 0).sum()),
        "mass_allowed_outside_catalogue": int((raw > 0).sum()),
    }
    # Poisson-disk thinning at the metric's own exclusion radius, priority = strike evidence
    pri = strike + 1e-6 * raw
    thinned = emission.thin_by_distance(raw > 0, min_sep=metric.R_PX, order=pri)
    log["emission"]["after_thinning_px"] = int((thinned > 0).sum())
    log["emission_seconds"] = round(time.time() - t0, 1)

    # ---------------------------------------------------------------- write + audit
    pred = np.clip(thinned, 0.0, 1.0).astype(np.float32)
    out = ROOT / "docs" / "downloads" / "gems45-h46-kmtip-eN-20261006.tif"
    grid.write_submission(pred, tmpl, out, outside="nan")
    a = grid.audit(out, tmpl)
    log["submission"] = a
    out_z = out.with_name(out.name.replace(".tif", "-zeros.tif"))
    z = np.where(tmpl, pred, np.float32(0.0))
    grid.write_submission(z, np.ones_like(tmpl), out_z, outside="zero")
    log["submission_zeros_twin"] = grid.audit(out_z, tmpl)

    log["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    log["total_seconds"] = round(time.time() - t_start, 1)
    (EV / "pipeline.json").write_text(json.dumps(log, indent=1))
    print(json.dumps({k: log[k] for k in
                      ("geometry", "catalogue", "n_tips", "kaplan_meier_genuine",
                       "extension_budget", "extension_lengths_px", "emission", "submission")},
                     indent=1)[:6000])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
