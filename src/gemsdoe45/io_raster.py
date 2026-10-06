"""Raster I/O and submission-format validation for the GEMS Prize.

The competition submission must (verified from the official submission-format
section of the problem page, https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/#submission-format):

  * Same CRS as training data: UTM zone 11N, EPSG:32611
  * Same resolution: 100 m
  * Same bounds; data outside the bounds is null or NaN
  * Single layer, 32-bit float (float32), values in [0,1]
  * A sample "total fault absence" submission is provided as the template

The feature cube (``training_features.tif`` / ``gems-geodawn-numerical-features.tif``)
stores its no-data sentinel as -3.4028e+38 (float32 minimum).  The reference
solution replaces values below -1e38 with NaN; we follow that exactly so that
edge/curvature operators are not corrupted by the sentinel.
"""

from __future__ import annotations

import json
import zipfile
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import rasterio as rio
from rasterio.enums import Compression

# float32 no-data sentinel used by the competition feature/label rasters.
SENTINEL = -3.4028234663852886e38
SENTINEL_FLOOR = -1.0e38

# Competition-fixed format (official submission-format section).
COMPETITION_CRS = "EPSG:32611"
COMPETITION_PIXEL_M = 100.0


def _to_nan(arr: np.ndarray) -> np.ndarray:
    """Replace the float32 sentinel with NaN (no copy if not needed)."""
    out = arr.astype(np.float32, copy=False)
    out[out < SENTINEL_FLOOR] = np.nan
    return out


def load_features(path: str | Path):
    """Load the multiband feature cube.

    Returns
    -------
    data : (H, W, C) float32 ndarray with the sentinel replaced by NaN
    descriptions : list[str]  band descriptions (may be empty)
    meta : dict  {'height','width','crs','transform'} copied from the file
    """
    path = str(path)
    with rio.open(path) as src:
        bands = src.read()  # (C, H, W)
        descriptions = []
        for i in range(1, src.count + 1):
            try:
                descriptions.append(src.tags(i).get("description", f"band{i}"))
            except Exception:
                descriptions.append(f"band{i}")
        meta = {
            "height": src.height,
            "width": src.width,
            "crs": str(src.crs),
            "transform": src.transform,
        }
    data = np.moveaxis(bands, 0, -1).astype(np.float32)
    data = _to_nan(data)
    return data, descriptions, meta


def load_template(path: str | Path):
    """Load the sample submission template.

    Returns
    -------
    template : (H, W) float32 array (NaNs outside the competition footprint)
    valid_mask : bool (H, W)  finite (in-footprint) pixels
    meta : dict
    """
    path = str(path)
    with rio.open(path) as src:
        arr = src.read(1).astype(np.float32)
        meta = {
            "height": src.height,
            "width": src.width,
            "crs": str(src.crs),
            "transform": src.transform,
        }
    valid = np.isfinite(arr)
    return arr, valid, meta


def load_faults(path: str | Path):
    """Load the known-fault raster as a binary mask (value == 1 -> fault).

    The file is int8 with nodata = -1; -1 and 0 are treated as non-fault.
    """
    path = str(path)
    with rio.open(path) as src:
        arr = src.read(1)
    return arr == 1


def validate_prediction_tiff(
    pred: np.ndarray,
    template_meta: dict,
    template_valid: np.ndarray,
    label: str = "candidate",
) -> dict:
    """Validate a predicted array against the competition submission format.

    This is the gate that prevents the ``Predicted values must be in range [0,1]``
    error: it checks the value range, the grid geometry, and the footprint mask.

    Returns a receipt dict with ``ok`` (bool) and per-check booleans.
    """
    pred = np.asarray(pred)
    receipt = {"label": label, "ok": True, "checks": {}}

    # 1. finite (in-footprint) value range must be within [0,1]
    inside = np.isfinite(pred)
    if inside.any():
        vmin = float(np.nanmin(pred[inside]))
        vmax = float(np.nanmax(pred[inside]))
    else:
        vmin, vmax = 0.0, 0.0
    range_ok = (vmin >= -1e-6) and (vmax <= 1.0 + 1e-6)
    receipt["checks"]["value_range_in_0_1"] = bool(range_ok)
    receipt["value_min"] = vmin
    receipt["value_max"] = vmax
    if not range_ok:
        receipt["ok"] = False

    # 2. shape matches the template
    shape_ok = (pred.shape == (template_meta["height"], template_meta["width"]))
    receipt["checks"]["shape_matches_template"] = bool(shape_ok)
    receipt["pred_shape"] = list(pred.shape)
    if not shape_ok:
        receipt["ok"] = False

    # 3. footprint mask matches the template (finite region == in-bounds region)
    if inside.shape == template_valid.shape:
        mask_ok = bool(np.array_equal(inside, template_valid))
    else:
        mask_ok = False
    receipt["checks"]["finite_mask_matches_template"] = mask_ok
    # Do NOT fail on mask mismatch alone if outside-values are 0 (competition
    # accepts null-or-zero outside bounds), but record it.
    if not mask_ok:
        # If every outside-template pixel is 0 (not NaN), it is still acceptable.
        outside = ~template_valid
        zero_outside_ok = bool(np.all(pred[outside] == 0) if outside.any() else True)
        receipt["checks"]["outside_is_zero_or_nan"] = zero_outside_ok
        if not zero_outside_ok:
            receipt["ok"] = False

    receipt["n_predicted_pixels"] = int((pred > 0).sum())
    receipt["coverage_fraction"] = float((pred > 0).mean())
    return receipt


def write_prediction_tiff(
    out_path: str | Path,
    pred: np.ndarray,
    template_path: str | Path,
    outside: str = "nan",
    compress: bool = True,
) -> dict:
    """Write a competition-compliant submission GeoTIFF.

    The output copies the template's CRS, transform, width/height and profile,
    writes ``pred`` (clipped to [0,1]) inside the template's finite footprint,
    and sets outside-footprint pixels to NaN (or 0 when ``outside='zero'``).

    Returns a validation receipt (see ``validate_prediction_tiff``).
    """
    out_path = Path(out_path)
    template_path = str(template_path)
    pred = np.asarray(pred, dtype=np.float32)

    with rio.open(template_path) as src:
        profile = src.profile.copy()
        template_arr = src.read(1)
        valid = np.isfinite(template_arr)

    # Clip to the required [0,1] range — this is what prevents the upload error.
    pred = np.clip(pred, 0.0, 1.0)
    out = np.full(pred.shape, np.nan, dtype=np.float32)
    out[valid] = pred[valid]
    if outside == "zero":
        out[~valid] = 0.0

    profile.update(
        dtype="float32",
        count=1,
        nodata=np.nan if outside == "nan" else 0.0,
        compress="deflate" if compress else "none",
        driver="GTiff",
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with rio.open(out_path, "w", **profile) as dst:
        dst.write(out, 1)

    receipt = validate_prediction_tiff(out, _meta_from_template(template_path), valid,
                                       label=out_path.name)
    receipt["path"] = str(out_path)
    receipt["outside_mode"] = outside
    receipt["written_at_utc"] = datetime.now(UTC).isoformat()
    return receipt


def _meta_from_template(template_path: str) -> dict:
    with rio.open(template_path) as src:
        return {
            "height": src.height,
            "width": src.width,
            "crs": str(src.crs),
            "transform": src.transform,
        }


def zip_submission(tif_path: str | Path, zip_path: str | Path) -> Path:
    """Package a single GeoTIFF into a .zip (competition-accepted container)."""
    tif_path = Path(tif_path)
    zip_path = Path(zip_path)
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as z:
        z.write(tif_path, tif_path.name)
    return zip_path
