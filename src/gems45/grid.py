"""Competition grid geometry and raster IO.

Geometry measured first-hand from the hash-pinned competition files (see registry/data_manifest.json):

    EPSG:32611 (UTM 11N), 100 m, 3292 cols x 3730 rows, origin (243350, 4508550)
    training_features.tif : 19 bands float32, nodata sentinel -3.4028235e+38
    labels.tif            : int8, values {-1 outside, 0 no-fault, 1 fault}, 60,988 fault px
    sample_submission.tif : float32, nodata NaN  <-- the organizers' own format template

Footprint, measured (IR-45-001): the finite mask of band 1 is 5,165,852 px, the finite mask of the
sample submission is 5,167,373 px, and only 5,165,840 px are finite in all 19 bands (band 6, `tc`,
carries 12 sentinel pixels inside the band-1 mask). The three masks are NOT identical. Because the
problem description designates the sample submission as the format template ("You can use this as a
template to ensure that your submission is correctly formatted"), the writer below uses the sample
submission's mask as the authority and additionally guarantees that every finite cell lies in [0, 1]
and that no sentinel value survives anywhere.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import rasterio

NODATA_F32 = np.float32(-3.4028234663852886e38)
HEIGHT, WIDTH = 3730, 3292
CRS = "EPSG:32611"
TRANSFORM = (100.0, 0.0, 243350.0, 0.0, -100.0, 4508550.0)

BAND_NAMES = [
    "mag_anom", "rtp", "tmi_hg", "geod_2ndinv", "iso_grav_anom_slope", "tc",
    "geod_shearrate", "geod_dilaterate", "tmi_vg", "deq_n100a15", "iso_grav_anom_vg",
    "det_elev", "iso_grav_anom", "tmi", "depth_to_base_surf", "ieq_n100a15",
    "cond_surf", "iso_grav_anom_hg", "det_elev_slope",
]


def read_bands(path: str | Path) -> np.ndarray:
    """Return the 19-band stack as float32 with sentinels replaced by the per-band median."""
    with rasterio.open(path) as s:
        a = s.read().astype(np.float32)
        nod = float(s.nodata)
    finite = np.ones(a.shape[1:], dtype=bool)
    for i in range(a.shape[0]):
        finite &= np.isfinite(a[i]) & (a[i] > nod * 0.5)
    for i in range(a.shape[0]):
        a[i][~finite] = np.float32(np.median(a[i][finite]))
    return a, finite


def read_labels(path: str | Path) -> np.ndarray:
    """True where the provided USGS/INGENIOUS catalogue maps a fault."""
    with rasterio.open(path) as s:
        return s.read(1) >= 1


def read_template_mask(path: str | Path) -> np.ndarray:
    """The sample submission's finite mask -- the authority for 'inside the footprint'."""
    with rasterio.open(path) as s:
        return np.isfinite(s.read(1))


def write_submission(pred: np.ndarray, template_mask: np.ndarray, out_path: str | Path,
                     outside: str = "nan") -> Path:
    """Write a submission GeoTIFF on the exact template grid.

    Guarantees, asserted before the file is written:
      * shape (3730, 3292), single band, float32, EPSG:32611, template transform
      * every finite cell lies in [0, 1] -- the portal error "Predicted values must be in
        range [0, 1]" cannot occur
      * no -3.4028235e38 sentinel survives anywhere in the file
      * outside is 'nan' (matches the organizers' template) or 'zero' (every cell in [0, 1])
    """
    pred = np.asarray(pred, dtype=np.float32)
    if pred.shape != (HEIGHT, WIDTH):
        raise ValueError(f"shape {pred.shape} != {(HEIGHT, WIDTH)}")
    mask = np.asarray(template_mask, dtype=bool)
    if mask.shape != (HEIGHT, WIDTH):
        raise ValueError(f"template mask shape {mask.shape} != {(HEIGHT, WIDTH)}")
    vals = pred[mask]
    if not np.isfinite(vals).all():
        raise ValueError("non-finite value inside the template footprint")
    if vals.min() < 0.0 or vals.max() > 1.0:
        raise ValueError(f"value outside [0, 1]: [{vals.min()}, {vals.max()}]")
    out = np.where(mask, pred, np.nan if outside == "nan" else np.float32(0.0))
    assert not np.any(out == NODATA_F32), "sentinel leaked into the submission"
    if outside == "nan":
        assert np.isnan(out[~mask]).all()

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    meta = dict(driver="GTiff", dtype="float32", count=1, crs=CRS, transform=TRANSFORM,
                width=WIDTH, height=HEIGHT, nodata=np.nan, compress="deflate", tiled=True)
    with rasterio.open(out_path, "w", **meta) as dst:
        dst.write(out, 1)
    return out_path


def audit(path: str | Path, template_mask: np.ndarray) -> dict:
    """Independent re-read audit of a written submission."""
    import hashlib
    p = Path(path)
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    with rasterio.open(p) as s:
        a = s.read(1)
        meta = dict(shape=s.shape, count=s.count, dtype=s.dtypes[0], crs=str(s.crs),
                    transform=tuple(float(x) for x in s.transform)[:6], nodata=str(s.nodata))
    m = np.asarray(template_mask, dtype=bool)
    inside = a[m]
    return {
        "filename": p.name,
        "bytes": p.stat().st_size,
        "sha256": h.hexdigest(),
        **meta,
        "bytes_match": p.stat().st_size,
        "finite_inside": int(np.isfinite(inside).sum()),
        "nan_inside": int(np.isnan(inside).sum()),
        "sentinel_anywhere": int((a == NODATA_F32).sum()),
        "inside_min": float(np.nanmin(inside)),
        "inside_max": float(np.nanmax(inside)),
        "inside_all_in_0_1": bool(np.all((inside >= 0.0) & (inside <= 1.0))),
        "outside_all_nan": bool(np.all(np.isnan(a[~m]))),
        "outside_all_in_0_1_or_nan": bool(np.all(np.isnan(a[~m]) | ((a[~m] >= 0) & (a[~m] <= 1)))),
        "n_positive": int((inside > 0).sum()),
        "n_at_one": int((inside >= 1.0).sum()),
    }
