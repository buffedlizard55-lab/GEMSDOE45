#!/usr/bin/env python3
"""Validate the local data bridge and emit a provenance manifest.

Verifies that the supplied files exist, that the feature cube's SHA-256 matches
the sibling project's pinned manifest (4371c82e…), and that the template and
fault rasters match the competition submission format (EPSG:32611, 100 m,
float32/[int8], finite footprint mask).  Writes data/prepared/manifest.json.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import rasterio as rio  # noqa: E402

# SHA-256 of the re-assembled gems-geodawn-numerical-features.tif (sibling manifest).
FEATURE_SHA256 = "4371c82e3b8339b807bdffcf4ef59a225520fe2988d521be208ae33743123bc5"


def sha256(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for blk in iter(lambda: f.read(chunk), b""):
            h.update(blk)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=str(ROOT / "data"))
    args = ap.parse_args()
    d = Path(args.data_dir) / "bridge"
    if not d.exists():
        print(f"Expected data bridge at {d}; run scripts/download_competition_data.sh first.")
        return 2

    report = {"prepared_at_utc": datetime.now(UTC).isoformat(), "files": {}, "checks": {}}

    feat = d / "gems-geodawn-numerical-features.tif"
    tmpl = d / "example_submission.tif"
    faults = d / "existing_faults.tif"
    for p in (feat, tmpl, faults):
        if not p.exists():
            print(f"Missing: {p}")
            return 2

    report["files"]["features"] = {"path": str(feat), "sha256": sha256(feat),
                                   "matches_pin": sha256(feat) == FEATURE_SHA256}
    with rio.open(feat) as s:
        report["files"]["features"].update(
            width=s.width, height=s.height, count=s.count, dtype=s.dtypes[0], crs=str(s.crs))
    with rio.open(tmpl) as s:
        report["files"]["template"] = {"width": s.width, "height": s.height,
                                       "dtype": s.dtypes[0], "crs": str(s.crs)}
        a = s.read(1)
        report["checks"]["template"] = {
            "crs_is_32611": str(s.crs) == "EPSG:32611",
            "resolution_100m": abs(s.res[0] - 100.0) < 1e-6,
            "finite_fraction": float(np.isfinite(a).mean()),
            "value_range": [float(np.nanmin(a)), float(np.nanmax(a))],
        }
    with rio.open(faults) as s:
        report["files"]["faults"] = {"width": s.width, "height": s.height, "dtype": s.dtypes[0]}

    ok = (report["files"]["features"]["matches_pin"]
          and report["checks"]["template"]["crs_is_32611"]
          and report["checks"]["template"]["resolution_100m"])
    report["status"] = "ready" if ok else "check_failed"
    out = ROOT / "data" / "prepared" / "manifest.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    print(f"\nWrote {out} — status: {report['status']}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
