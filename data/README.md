# Local data

This directory is ignored except this notice. Run `bash scripts/download_competition_data.sh` from the repository after installing requirements. Mirrors are pinned in `research/prior_data_manifest.json`; evidence records integrity, not organizer authentication. Raw rasters, features.npy (~2.3 GB), truth/footprint arrays and prior-TIF scratch files must not be committed. Small final submission artifacts belong under docs/downloads/.


---

# data/

The competition inputs are **not tracked in git** (they are ~420 MB). They are mirrored publicly and
were fetched here with the GitHub API and verified against pinned SHA-256 digests. To reproduce:

```bash
# 425,830 B each -- labels.tif and existing_faults.tif are byte-identical
gh api "repos/buffedlizard55-lab/GEMSDOE24/contents/data/bridge/labels.tif?ref=07345ea0604953d7efb858d9cfbc21e20c7aca0b" \
  -H "Accept: application/vnd.github.raw" > data/labels.tif
# expect sha256 7ba308ccdc4418b31a178f4f1ef21aaa6e152e4028f2f6f64b01f7eb25ae4093

gh api "repos/buffedlizard55-lab/GEMSDOE24/contents/data/bridge/sample_submission.tif?ref=07345ea0604953d7efb858d9cfbc21e20c7aca0b" \
  -H "Accept: application/vnd.github.raw" > data/sample_submission.tif
# expect sha256 2176d08e485aa2cd2860ce8df539db4faf4d76163b38a4dd8c30a40454d35cbc

# 418,912,844 B, published in five parts
for i in 000 001 002 003 004; do
  gh api "repos/buffedlizard55-lab/GEMSDOE/contents/data/bridge/gems-geodawn-numerical-features.tif.part-$i?ref=c0c06ac82178f26b94fce3397036ef8f12a2f3a0" \
    -H "Accept: application/vnd.github.raw" > /tmp/part-$i
done
cat /tmp/part-000 /tmp/part-001 /tmp/part-002 /tmp/part-003 /tmp/part-004 > data/training_features.tif
# expect sha256 4371c82e3b8339b807bdffcf4ef59a225520fe2988d521be208ae33743123bc5
```

**Provenance warning.** These are owner-supplied mirrors of login-walled DrivenData competition
files. The hashes prove the bytes are the ones the mirrors publish; they do **not** constitute
organizer authentication. The organizers' own download page is
<https://www.drivendata.org/competitions/306/competition-doe-gems/data/> and requires an account.

Measured geometry: EPSG:32611, 100 m, 3730 rows x 3292 columns, origin (243350, 4508550).
`labels.tif`: int8, nodata -1, values {-1 outside, 0 no-fault, 1 fault}, 60,988 fault pixels.
`training_features.tif`: 19 float32 bands, nodata sentinel -3.4028235e+38.
`sample_submission.tif`: float32, nodata NaN -- the organizers' designated format template.
