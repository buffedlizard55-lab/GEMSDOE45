"""The writer's guarantees: a file that cannot trigger the portal's range error."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from gems45 import grid, metric  # noqa: E402


def test_write_submission_roundtrip(tmp_path):
    tmpl = np.zeros((grid.HEIGHT, grid.WIDTH), dtype=bool)
    tmpl[100:140, 100:160] = True
    pred = np.zeros(tmpl.shape, dtype=np.float32)
    pred[110:114, 110:150] = 1.0
    out = grid.write_submission(pred, tmpl, tmp_path / "s.tif", outside="nan")
    a = grid.audit(out, tmpl)
    assert tuple(a["shape"]) == (grid.HEIGHT, grid.WIDTH)
    assert a["count"] == 1 and a["dtype"] == "float32"
    assert a["crs"] == grid.CRS
    assert a["inside_all_in_0_1"] is True
    assert a["outside_all_nan"] is True
    assert a["sentinel_anywhere"] == 0
    assert a["n_positive"] == 4 * 40


def test_writer_rejects_out_of_range(tmp_path):
    tmpl = np.zeros((grid.HEIGHT, grid.WIDTH), dtype=bool)
    tmpl[10, 10] = True
    pred = np.zeros(tmpl.shape, dtype=np.float32)
    pred[10, 10] = 1.5
    with pytest.raises(ValueError, match="outside"):
        grid.write_submission(pred, tmpl, tmp_path / "bad.tif")


def test_writer_rejects_sentinel(tmp_path):
    tmpl = np.zeros((grid.HEIGHT, grid.WIDTH), dtype=bool)
    tmpl[10, 10] = True
    pred = np.zeros(tmpl.shape, dtype=np.float32)
    pred[10, 10] = grid.NODATA_F32
    with pytest.raises(ValueError):
        grid.write_submission(pred, tmpl, tmp_path / "bad2.tif")


def test_zero_copy_has_no_nan_anywhere(tmp_path):
    tmpl = np.zeros((grid.HEIGHT, grid.WIDTH), dtype=bool)
    tmpl[50:60, 50:60] = True
    pred = np.zeros(tmpl.shape, dtype=np.float32)
    pred[55, 55] = 1.0
    out = grid.write_submission(pred, np.ones_like(tmpl), tmp_path / "z.tif", outside="zero")
    a = grid.audit(out, tmpl)
    assert a["outside_all_nan"] is False
    assert a["outside_all_in_0_1_or_nan"] is True
    assert a["n_positive"] == 1


def test_marginal_bar_and_radius():
    assert metric.marginal_bar(0.30) == pytest.approx(0.06)
    assert metric.profitable_radius_m(0.30) == pytest.approx(282.0)
