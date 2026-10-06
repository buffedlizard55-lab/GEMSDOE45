"""Measure, for every band *and* every derived plane, whether its strongest 40,000 pixels land on
the catalogue.

This is the screen that decides LOCATOR vs WEIGHT in ``bands.py``, and it is the screen a reader of
the site will check.  It is a script rather than an inline snippet because the site publishes its
numbers, and an unbacked number in a table is a liability.

Two design points matter.

1. **Like-for-like on mass.**  Every plane is reduced to the same 40,000 positive pixels as the
   shipped submission, so the column is comparable across planes and comparable to the submission.
   A plane that only looks good at 200,000 pixels has not been measured here.

2. **The baseline is computed, not asserted.**  The control is the fraction of *footprint* pixels
   that lie within 3 px (300 m) of a catalogue pixel.  Quoting a plane's fraction without that
   baseline is meaningless.

The derivative planes are taken from ``gems45.detfeatures.iter_features`` -- the same generator the
shipped field is built from -- so this document cannot drift from the artifact.

Output: ``evidence/band_screen.json``.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from gems45 import bands, detfeatures, grid  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TARGET_MASS = 40_000
KERNEL_PX = 3  # 300 m at 100 m/px


def topk_mask(plane: np.ndarray, foot: np.ndarray, k: int) -> np.ndarray:
    """Boolean mask of the k largest finite-in-footprint values, deterministic under ties."""
    v = np.where(foot & np.isfinite(plane), np.abs(plane).astype(np.float64), -np.inf).ravel()
    cut = np.partition(v, -k)[-k]
    sel = np.flatnonzero(v >= cut)
    if sel.size > k:  # ties at the cut: keep the lowest indices
        sel = sel[:k]
    m = np.zeros(v.size, dtype=bool)
    m[sel] = True
    return m.reshape(plane.shape)


def score_plane(plane, foot, dist, k, baseline):
    """Skill of a plane's strongest k pixels, plus a test of whether that skill is *localised*.

    A regionally smooth field can post a good ``within_300m`` fraction for an uninteresting reason:
    its top k pixels are then one or two large connected blobs, and a blob that happens to sit over
    the fault-rich half of the map inherits the catalogue's own density.  That is regional
    information, not fault localisation, and a 300 m kernel cannot use it.  ``n_components`` and the
    perimeter-to-area ratio separate the two cases: a kernel-scale localiser paints many small,
    elongated pieces (high ratio), a regional field paints one fat blob (ratio < 0.05).
    """
    m = topk_mask(plane, foot, k)
    d = dist[m]
    frac = float((d <= KERNEL_PX).mean())
    lab, n = ndimage.label(m, structure=np.ones((3, 3), dtype=int))
    sizes = np.bincount(lab.ravel())[1:]
    per_area = float((m & ~ndimage.binary_erosion(m)).sum()) / int(m.sum())
    return {
        "mass_px": int(m.sum()),
        "within_300m_frac": round(frac, 4),
        "within_300m_pct": round(100.0 * frac, 1),
        "enrichment_vs_random": round(frac / baseline, 2),
        "median_dist_px": round(float(np.median(d)), 1),
        "n_components": int(n),
        "median_component_px": int(np.median(sizes)),
        "largest_component_px": int(sizes.max()),
        "perimeter_over_area": round(per_area, 3),
        "localised": bool(n >= 100 and per_area >= 0.30),
    }


def main() -> int:
    t0 = time.time()
    z, foot = grid.read_bands(ROOT / "data/training_features.tif")
    known = grid.read_labels(ROOT / "data/labels.tif") & foot
    names = list(grid.BAND_NAMES)

    dist = ndimage.distance_transform_edt(~known)
    baseline = float(((dist <= KERNEL_PX) & foot).sum()) / float(foot.sum())

    # ---- raw band values -------------------------------------------------------------------
    raw = []
    for i, name in enumerate(names):
        r = {"band": name, "index": i + 1, **score_plane(z[i], foot, dist, TARGET_MASS, baseline)}
        raw.append(r)
    print(f"raw band values screened ({len(raw)} rows)")

    # ---- derived planes, canonical definitions from detfeatures -----------------------------
    planes = []
    for name, plane in detfeatures.iter_features(z, foot):
        r = {"plane": name, **score_plane(plane, foot, dist, TARGET_MASS, baseline)}
        planes.append(r)
        if len(planes) % 10 == 0:
            print(f"  ... {len(planes)} planes")
        del plane
    print(f"derived planes screened ({len(planes)} rows)")

    # ---- reconcile: the measured class rule vs the explicit mass-placement lists --------------
    info = {d["band"]: d for d in bands.band_information(z, foot, known)}
    elig = set(bands.LOCATOR_BANDS)
    for r in raw:
        r.update(rho_lag10px=info[r["band"]]["rho_lag10px"],
                 mean_abs_grad_sigma=info[r["band"]]["mean_abs_grad_sigma"],
                 auc_value_vs_catalogue=info[r["band"]]["auc_value_vs_catalogue"],
                 band_class=info[r["band"]]["band_class"],
                 mass_eligible=r["band"] in elig)
    disagree = [{"band": r["band"], "measured_class": r["band_class"], "mass_eligible": r["mass_eligible"],
                 "within_300m_pct": r["within_300m_pct"]}
                for r in raw if r["mass_eligible"] != (r["band_class"] == "LOCATOR")]

    out = {
        "target_mass_px": TARGET_MASS,
        "kernel_px": KERNEL_PX,
        "catalogue_px": int(known.sum()),
        "footprint_px": int(foot.sum()),
        "random_within_300m_pct": round(100.0 * baseline, 2),
        "raw_bands": sorted(raw, key=lambda r: -r["within_300m_frac"]),
        "derived_planes": sorted(planes, key=lambda r: -r["within_300m_frac"]),
        "class_disagreements": disagree,
        "seconds": round(time.time() - t0, 1),
    }
    p = ROOT / "evidence/band_screen.json"
    p.write_text(json.dumps(out, indent=1))

    print(f"\nrandom baseline {100.0*baseline:.2f} %   catalogue {int(known.sum())} px   "
          f"footprint {int(foot.sum())} px")
    print("\ntop 8 raw bands by top-40k within 300 m:")
    for r in out["raw_bands"][:8]:
        print(f"  {r['band']:<22}{r['within_300m_pct']:>6} %  {r['enrichment_vs_random']:>5}x")
    print("\ntop 8 derived planes by top-40k within 300 m:")
    for r in out["derived_planes"][:8]:
        print(f"  {r['plane']:<30}{r['within_300m_pct']:>6} %  {r['enrichment_vs_random']:>5}x")
    print(f"\nwrote {p} in {out['seconds']} s")
    for d in disagree:
        print(f"  IRREGULARITY: {d['band']} measured {d['measured_class']} but "
              f"mass_eligible={d['mass_eligible']} ({d['within_300m_pct']} % within 300 m)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
