"""Score every stored candidate *field* with the exact instrument used for the bands.

The site quotes enrichment numbers for the band table and for the two ensembles (the H47 consensus
field and the shipped H48 rank-ensemble field).  If those come from two different instruments they
are not comparable, and a reader cannot check either.  This script runs one instrument over all of
them and writes ``evidence/field_screen.json``.

The instrument is inherited from ``scripts/screen_bands.py``: the strongest k pixels by field value,
the fraction of them within 3 px (300 m) of the catalogue, against the measured random baseline --
plus the fragmentation test that separates a real localiser from a regional blob.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from gems45 import grid  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
KERNEL_PX = 3
MASSES = (20_000, 40_000, 80_000)


def screen(field: np.ndarray, foot: np.ndarray, dist: np.ndarray, baseline: float) -> dict:
    out = {}
    v = np.where(foot & np.isfinite(field), field.astype(np.float64), -np.inf).ravel()
    for k in MASSES:
        cut = np.partition(v, -k)[-k]
        sel = np.flatnonzero(v >= cut)
        if sel.size > k:
            sel = sel[:k]
        m = np.zeros(v.size, dtype=bool)
        m[sel] = True
        m = m.reshape(field.shape)
        d = dist[m]
        frac = float((d <= KERNEL_PX).mean())
        lab, n = ndimage.label(m, structure=np.ones((3, 3), dtype=int))
        sizes = np.bincount(lab.ravel())[1:]
        out[f"top_{k//1000}k"] = {
            "mass_px": int(m.sum()),
            "within_300m_frac": round(frac, 4),
            "within_300m_pct": round(100.0 * frac, 1),
            "enrichment_vs_random": round(frac / baseline, 2),
            "n_components": int(n),
            "median_component_px": int(np.median(sizes)),
            "perimeter_over_area": round(float((m & ~ndimage.binary_erosion(m)).sum()) / int(m.sum()), 3),
        }
    return out


def main() -> int:
    z, foot = grid.read_bands(ROOT / "data/training_features.tif")
    known = grid.read_labels(ROOT / "data/labels.tif") & foot
    dist = ndimage.distance_transform_edt(~known)
    baseline = float(((dist <= KERNEL_PX) & foot).sum()) / float(foot.sum())

    fields = {}
    for name, path in (("H48_rank_ensemble", "evidence/field.npy"),
                       ("H47_consensus", "evidence/consensus.npy")):
        p = ROOT / path
        if not p.exists():
            print(f"  skip {name}: {path} absent")
            continue
        fields[name] = {"source": path, **screen(np.load(p), foot, dist, baseline)}
        print(f"  {name}: top-40k {fields[name]['top_40k']['within_300m_pct']} % "
              f"({fields[name]['top_40k']['enrichment_vs_random']}x)")

    # The shipped file, read back from disk: the ground truth for every claim about it.
    shipped = ROOT / "docs/downloads/gems45-h48-kmtip-structural-20261006.tif"
    if shipped.exists():
        pred = grid.read_template_mask(ROOT / "data/sample_submission.tif")
        with __import__("rasterio").open(shipped) as s:
            a = np.nan_to_num(s.read(1), nan=0.0)
        m = a > 0
        d = dist[m]
        frac = float((d <= KERNEL_PX).mean())
        lab, n = ndimage.label(m, structure=np.ones((3, 3), dtype=int))
        sizes = np.bincount(lab.ravel())[1:]
        fields["H48_shipped_file"] = {
            "source": str(shipped.relative_to(ROOT)), "mass_px": int(m.sum()),
            "within_300m_frac": round(frac, 4), "within_300m_pct": round(100.0 * frac, 1),
            "enrichment_vs_random": round(frac / baseline, 2),
            "n_components": int(n), "median_component_px": int(np.median(sizes)),
            "perimeter_over_area": round(float((m & ~ndimage.binary_erosion(m)).sum()) / int(m.sum()), 3),
        }
        print(f"  H48_shipped_file: {fields['H48_shipped_file']['within_300m_pct']} % "
              f"({fields['H48_shipped_file']['enrichment_vs_random']}x)")

    out = {"kernel_px": KERNEL_PX, "random_within_300m_pct": round(100.0 * baseline, 2),
           "catalogue_px": int(known.sum()), "footprint_px": int(foot.sum()), "fields": fields}
    p = ROOT / "evidence/field_screen.json"
    p.write_text(json.dumps(out, indent=1))
    print(f"wrote {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
