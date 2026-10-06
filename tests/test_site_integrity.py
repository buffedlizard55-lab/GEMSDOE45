"""Guard the claims the site makes about files on disk.

The point of these tests is that the site is *generated* data, not prose.  Every published hex digest
and every published band number is checked against the bytes or the evidence file it claims to come
from, so a silent drift fails CI instead of shipping.
"""
from __future__ import annotations

import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from conftest import require_competition_data

ROOT = Path(__file__).resolve().parents[1]
SITE = ROOT / "docs/e45/index.html"

SHIPPED = {
    "gems45-h48-kmtip-structural-20261006.tif":
        "40eb32bf995433b15052745c7ea647ed344dc4d56ca9356f8b843a4c1466a402",
    "gems45-h48b-structural-noflank-20261006.tif":
        "be7915f2f68bb077b486a3c210c6fe6ee6c831f5bc21a56761613906f532ea8c",
}


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.mark.parametrize("name,expected", list(SHIPPED.items()))
def test_published_digest_matches_the_bytes(name, expected):
    p = ROOT / "docs/downloads" / name
    if not p.exists():
        pytest.skip(f"{name} not present")
    assert _sha(p) == expected, f"{name} has changed; update the site and this test together"


@pytest.mark.parametrize("name", list(SHIPPED))
def test_every_published_digest_is_the_real_one(name):
    """The digest in the HTML must be the digest of the file, whatever the file is."""
    p = ROOT / "docs/downloads" / name
    if not p.exists():
        pytest.skip(f"{name} not present")
    assert _sha(p) in SITE.read_text(), f"{name}'s true sha256 is not published on the site"


@pytest.mark.parametrize("name", list(SHIPPED))
def test_submission_is_metric_legal(name):
    rasterio = pytest.importorskip("rasterio")
    require_competition_data()
    p = ROOT / "docs/downloads" / name
    if not p.exists():
        pytest.skip(f"{name} not present")
    sys.path.insert(0, str(ROOT / "src"))
    from gems45 import grid
    with rasterio.open(ROOT / "data/sample_submission.tif") as t:
        ref = (t.shape, str(t.crs), tuple(t.transform)[:6])
    with rasterio.open(p) as s:
        assert (s.shape, str(s.crs), tuple(s.transform)[:6]) == ref
        assert s.count == 1 and s.dtypes[0] == "float32"
        a = s.read(1)
    assert not (a == grid.NODATA_F32).any(), "sentinel present"
    fin = np.isfinite(a)
    assert fin.any()
    assert a[fin].min() >= 0.0 and a[fin].max() <= 1.0, "values outside [0, 1]"
    tmpl = grid.read_template_mask(ROOT / "data/sample_submission.tif")
    assert np.isnan(a[~tmpl]).all(), "NaN required outside the template footprint"


def test_site_band_tables_are_in_sync_with_evidence():
    r = subprocess.run([sys.executable, str(ROOT / "scripts/sync_site_tables.py"), "--check"],
                       capture_output=True, text=True, cwd=ROOT)
    assert r.returncode == 0, r.stdout + r.stderr


def test_published_enrichment_baseline_is_measured_not_guessed():
    ev = json.loads((ROOT / "evidence/band_screen.json").read_text())
    html = SITE.read_text()
    assert f'{ev["random_within_300m_pct"]} %' in html, "the measured random baseline is not on the page"
    # the old hand-typed baseline must be gone
    assert "8.3 % random" not in html


def test_flank_ab_is_published_with_both_files():
    html = SITE.read_text()
    for n in SHIPPED:
        assert n in html, f"{n} is not offered on the site"
    assert 'id="flank"' in html
