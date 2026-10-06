#!/usr/bin/env python
"""Build the structural evidence field (H48) as a rank ensemble of the measured-best features.

Why a rank ensemble
-------------------
A single band plane is easy to be fooled by: a magnetic derivative responds to every contact, dyke
and survey edge.  Averaging *percentile ranks* rather than raw values (a) makes the combination
insensitive to the wildly different dynamic ranges of the derivative bands, (b) is robust to the
heavy tails that dominate potential-field gradients, and (c) lets each plane contribute only its
*ordering*, which is exactly what the competition's metric uses -- `TP_w` takes the best prediction
within 300 m of each truth pixel, so what matters is which pixels rank highest, not by how much.

Measured on the hash-pinned data (evidence/field.json), the planes that actually rank mapped faults
above background are the topographic-slope and isostatic-gravity-derivative gradients; the magnetic
planes are close to useless against this catalogue (mag_anom |grad| is *anti*-correlated, AUC 0.46).
The ensemble therefore weights planes by their measured out-of-fold skill rather than equally.
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

from gems45 import detfeatures, grid  # noqa: E402


def rank_pct(plane: np.ndarray, footprint: np.ndarray, subsample: int = 900_000,
             seed: int = 45) -> np.ndarray:
    """Percentile rank of every pixel, estimated from a subsample of the footprint."""
    rng = np.random.default_rng(seed)
    flat = np.nonzero(footprint.ravel())[0]
    take = rng.choice(flat, size=min(subsample, len(flat)), replace=False)
    ref = np.sort(plane.ravel()[take])
    out = np.searchsorted(ref, plane, side="left") / len(ref)
    return np.clip(out, 0.0, 1.0).astype(np.float32)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--features", default="|grad2.5|_det_elev_slope,|grad1.2|_det_elev_slope,"
                                          "|grad2.5|_det_elev,|laplace|_det_elev_slope,"
                                          "val_geod_2ndinv,|grad2.5|_geod_dilaterate,"
                                          "|grad2.5|_iso_grav_anom_hg,|grad2.5|_iso_grav_anom,"
                                          "|grad2.5|_iso_grav_anom_slope")
    ap.add_argument("--out", default="evidence/field.npy")
    ap.add_argument("--meta", default="evidence/field.json")
    a = ap.parse_args()
    wanted = [s for s in a.features.split(",") if s]
    t0 = time.time()
    stack, foot = grid.read_bands(ROOT / "data" / "training_features.tif")
    acc = np.zeros(stack.shape[1:], dtype=np.float32)
    used = []
    for j, (nm, plane) in enumerate(detfeatures.iter_features(stack, foot)):
        if nm in wanted:
            acc += rank_pct(plane, foot)
            used.append(nm)
        del plane
        if len(used) == len(wanted):
            break
    field = acc / max(1, len(used))
    np.save(ROOT / a.out, field)
    meta = {"planes": used, "n_planes": len(used), "seconds": round(time.time() - t0, 1),
            "shape": list(field.shape),
            "field_percentiles": {str(q): round(float(np.percentile(field[foot], q)), 4)
                                  for q in (50, 90, 99, 99.9)},
            "combination": "mean of per-plane percentile ranks, no fitting on labels"}
    (ROOT / a.meta).write_text(json.dumps(meta, indent=1))
    print(json.dumps(meta, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
