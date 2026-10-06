"""Guard the claims the H51 pages and README make about the H51 artifact.

The point of these tests is the same as ``test_site_integrity.py``'s: the site is *generated* data,
so every published digest and every published count is checked against the bytes or the evidence
file it claims to come from.  In particular the brief requires the extension lengths to be shown to
vary per fault *before* the file is downloaded, and that claim is asserted here rather than trusted.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence"
DL = ROOT / "docs" / "downloads"
SUB = json.loads((EV / "h51_submission.json").read_text())
SITE = ROOT / "docs/index.html"


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def test_extension_lengths_vary_per_fault():
    ext = SUB["extension_lengths"]
    assert ext["drawn_distinct"] >= 50, "per-tip draws are effectively constant"
    assert ext["realised_distinct"] >= 5, "realised extension lengths collapsed to a constant"
    assert ext["n_tips_with_realised_extension"] > 0.5 * ext["n_tips"]
    # a spread, not one value with one outlier
    realised = {int(k): v for k, v in ext["realised_histogram"].items() if int(k) > 0}
    assert len(realised) >= 5
    assert max(realised.values()) < sum(realised.values()), "one value dominates the histogram"


def test_both_published_digests_match_the_bytes():
    for key in ("nan_outside", "all_finite"):
        f = DL / SUB["files"][key]["filename"]
        assert f.exists(), f"{f} is not published"
        assert _sha(f) == SUB["files"][key]["sha256"], f"{f.name} has changed since the evidence was written"


def test_zip_contains_exactly_the_all_finite_tif_and_matches_it():
    import zipfile
    z = DL / SUB["files"]["zip"]["name"]
    assert _sha(z) == SUB["files"]["zip"]["sha256"]
    with zipfile.ZipFile(z) as zf:
        assert zf.namelist() == [SUB["files"]["all_finite"]["filename"]]
        assert hashlib.sha256(zf.read(zf.namelist()[0])).hexdigest() == SUB["files"]["all_finite"]["sha256"]


def test_published_digest_is_on_the_site():
    html = SITE.read_text()
    assert SUB["files"]["all_finite"]["sha256"] in html or SUB["files"]["nan_outside"]["sha256"] in html
    assert SUB["name"] in html


def test_dot_count_is_the_audited_positive_pixel_count_and_within_the_cap():
    assert SUB["n_dots"] <= SUB["rules"]["max_dots"]
    for key in ("nan_outside", "all_finite"):
        assert SUB["files"][key]["n_positive"] == SUB["n_dots"]
    assert SUB["emission_composition"]["tip_continuation"] + \
        SUB["emission_composition"]["strand_offset"] == SUB["n_dots"]


def test_no_score_is_claimed_and_the_slot_gate_is_closed():
    assert SUB["slot_approved"] is False
    assert SUB["score_claim"] is None
    assert SUB["emulator_is_not_a_score"] is True
    assert "NOT score predictions" in SUB["model_caveat"]


def test_emulator_sanity_the_assertions_that_fail_the_build():
    g = SUB["greedy"]
    assert g["f_pred"] >= 0.0, "an emulator false-positive mass cannot be negative"
    assert g["t_pred"] <= g["k_pred"] + 1e-6, "emulator T cannot exceed the modelled truth size"


def test_live_fit_does_not_claim_to_predict_the_leaderboard():
    fit = json.loads((EV / "h51_live_fit.json").read_text())
    assert fit["best"]["pearson"] < 0.5, "the emulator's field must not be presented as predictive"
    assert fit["n"] >= 40


@pytest.mark.parametrize("name", ["gems45-h51-km-faultzone-20261006.tif",
                                  "gems45-h51-km-faultzone-20261006-zeros.tif"])
def test_inside_the_catalogue_guarantee_is_actually_held(name):
    """The design claims: nothing on the mapped catalogue, nothing within 200 m of it."""
    rasterio = pytest.importorskip("rasterio")
    from conftest import require_competition_data
    require_competition_data()
    with rasterio.open(DL / name) as s:
        dots = s.read(1) > 0
    with rasterio.open(ROOT / "data/labels.tif") as s:
        catalog = s.read(1) >= 1
    assert not (dots & catalog).any(), "a dot sits on a mapped fault pixel"
    from scipy import ndimage
    d = ndimage.distance_transform_edt(~catalog)
    assert float(d[dots].min()) >= SUB["rules"]["r_zero_px"] - 1e-6, \
        "a dot sits inside the calibrated catalogue-exclusion radius"


def test_submission_is_format_legal_and_matches_the_template():
    rasterio = pytest.importorskip("rasterio")
    from conftest import require_competition_data
    require_competition_data()
    sys.path.insert(0, str(ROOT / "src"))
    from gems45 import grid
    with rasterio.open(ROOT / "data/sample_submission.tif") as t:
        ref = (t.shape, str(t.crs), tuple(t.transform)[:6])
    for name, outside_nan in ((SUB["files"]["nan_outside"]["filename"], True),
                              (SUB["files"]["all_finite"]["filename"], False)):
        with rasterio.open(DL / name) as s:
            assert (s.shape, str(s.crs), tuple(s.transform)[:6]) == ref
            assert s.count == 1 and s.dtypes[0] == "float32"
            a = s.read(1)
        assert not (a == grid.NODATA_F32).any(), "sentinel present anywhere"
        tmpl = grid.read_template_mask(ROOT / "data/sample_submission.tif")
        assert np.isnan(a[~tmpl]).all() if outside_nan else np.isfinite(a).all()
        fin = np.isfinite(a)
        assert a[fin].min() >= 0.0 and a[fin].max() <= 1.0


def test_no_prior_artifact_has_the_same_pixels():
    fp = json.loads((EV / "prior_fingerprints.json").read_text())["files"]
    comp = [f for f in fp if f["status"] == "fingerprinted"]
    report = json.loads((EV / "h51_validation.json").read_text())
    assert report["v3_uniqueness"]["compared"] == len(comp)
    assert report["v3_uniqueness"]["identical_float32_sha256"] == 0
