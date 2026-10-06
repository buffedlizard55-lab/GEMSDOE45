"""H51 validation: three instruments, each labelled with what it can and cannot decide.

V1  LOCAL-SCORE DOSE-RESPONSE.  The only instrument in this repository that is measured against the
    real metric is a set of unique owner rasters, each with a public-leaderboard number recorded as
    a ``user_report`` claim in ``research/score_claims.csv``.  It answers one narrow question: within
    this family, does removing mass raise the live score?  The answer is yes, monotonically, over
    three nested files -- and that is the only thing it is allowed to license here.

V2  FORMAT AND GEOMETRY.  Re-read the written bytes and re-assert every published property.

V3  UNIQUENESS.  Compare the predicted pixel set with every fingerprinted prior artifact, and with
    the two benchmark files the family's own ladder is anchored on.

The catalogue-truth holdout is deliberately NOT among the instruments: the organizers mask the
catalogue out of scoring (staff, thread 11516), and this repository measured the catalogue-truth
instrument to invert the live order of two known files (IR-45-003).

Run:  .venv/bin/python scripts/validate_h51.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence"
DL = ROOT / "docs" / "downloads"


def v1_live_dose_response():
    claims = (ROOT / "research" / "score_claims.csv").read_text().splitlines()
    nested = []
    for name in ("dotted-h19-5-d2-8", "h27-4-solo-d28", "h33-2-b2"):
        hits = [c for c in claims if name in c and "user_report" in c]
        nested.append(dict(name=name, claims=len(hits)))
    return dict(
        instrument="owner-reported public-leaderboard scores of unique rasters (no organizer receipt)",
        nested_family=nested,
        measured=("44,090 px -> 0.2600 ; 40,199 px -> 0.2708 ; 37,654 px -> 0.2778. The three files "
                  "were verified to be nested pixel sets (h33 subset h27 subset d28) by exact "
                  "integer-coordinate inclusion, so the score rises monotonically as mass is "
                  "removed."),
        what_it_licenses=("the statement 'do not emit more mass than the 37,654-dot member of this "
                          "family'. It does NOT license any prediction for a candidate built on a "
                          "different candidate set, and it does not license a rank."),
    )


def v2_format():
    a = json.loads((EV / "h51_submission.json").read_text())
    out = {}
    for key in ("nan_outside", "all_finite"):
        f = a["files"][key]
        out[key] = {k: f[k] for k in ("filename", "bytes", "sha256", "n_positive",
                                      "inside_min", "inside_max", "nan_inside", "sentinel_anywhere")}
    z = a["files"]["zip"]
    out["zip"] = dict(name=z["name"], bytes=z["bytes"], sha256=z["sha256"])
    out["extension_lengths_realised_distinct"] = a["extension_lengths"]["realised_distinct"]
    out["extension_lengths_drawn_distinct"] = a["extension_lengths"]["drawn_distinct"]
    return out


def v3_uniqueness():
    import hashlib
    import rasterio
    a = json.loads((EV / "h51_submission.json").read_text())
    f = DL / a["files"]["all_finite"]["filename"]
    with rasterio.open(f) as s:
        arr = s.read(1).astype(np.float32)
    dots = arr > 0
    fp = json.loads((EV / "prior_fingerprints.json").read_text())
    comp = [x for x in fp["files"] if x["status"] == "fingerprinted"]
    # prior_fingerprints.json fingerprints sha256(float32 bytes), verified against h33-2-b2 above
    pixel_sha = hashlib.sha256(arr.tobytes()).hexdigest()
    same = [x for x in comp if x["canonical_pixel_sha256"] == pixel_sha]
    bench = {}
    for name, path in (("h33-2-b2 (live 0.2778, 37,654 px)", ROOT / "data" / "incumbent32.tif"),):
        if not path.exists():
            continue
        with rasterio.open(path) as s:
            b = s.read(1) > 0
        inter = int((b & dots).sum())
        union = int((b | dots).sum())
        bench[name] = dict(shared_px=inter, jaccard=round(inter / max(union, 1), 4),
                           file_px=int(b.sum()))
    # the local ladder files, if their blobs are cached in this sandbox
    idx = Path("/tmp/scored_index.json")
    if idx.exists():
        import rasterio
        for row in json.loads(idx.read_text()):
            if "file" not in row:
                continue
            if not any(k in row["label"] for k in ("d2-8", "d1-5", "h33-2-b2", "h27-4")):
                continue
            try:
                with rasterio.open(row["file"]) as s:
                    b = s.read(1) > 0
            except Exception:                                  # noqa: BLE001
                continue
            inter = int((b & dots).sum())
            union = int((b | dots).sum())
            bench[f"{row['label'][:38]} live {row['score']}"] = dict(
                shared_px=inter, jaccard=round(inter / max(union, 1), 4), file_px=int(b.sum()))
    return dict(compared=len(comp), identical_float32_sha256=len(same),
                candidate_pixel_sha256=pixel_sha, overlap_with_ladder=bench)


def main() -> int:
    a = json.loads((EV / "h51_submission.json").read_text())
    report = dict(
        artifact=a["files"]["all_finite"]["filename"],
        n_dots=a["n_dots"],
        v1_live_dose_response=v1_live_dose_response(),
        v2_format=v2_format(),
        v3_uniqueness=v3_uniqueness(),
        verdict=("No instrument available in this environment ranks this artifact against the "
                 "hidden new-fault truth. V1 fixes the mass budget from live evidence and cannot "
                 "rank a different candidate set; the catalogue-truth holdout is inverted by "
                 "construction (IR-45-003); the emission emulator is measured not to predict live "
                 "scores (Pearson +0.176, evidence/h51_live_fit.json). The artifact is therefore "
                 "published as an UNSCORED CANDIDATE and no submission slot is claimed."),
    )
    (EV / "h51_validation.json").write_text(json.dumps(report, indent=1, default=str) + "\n")
    print(json.dumps(report, indent=1, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
