#!/usr/bin/env bash
# Reproduce every number quoted on the site, in order. ~5 minutes on CPU, no GPU, no network.
set -euo pipefail
PY=${PY:-/home/user/.venv/bin/python}
cd "$(dirname "$0")"
$PY scripts/run_pipeline.py            # metric checks, band values, catalogue, Kaplan-Meier
$PY scripts/screen_bands.py            # per-band / per-plane localisation + fragmentation screen
$PY scripts/score_fields.py            # every stored field on one instrument -> evidence/field_screen.json
$PY scripts/sync_site_tables.py        # regenerate the site tables from the evidence
$PY scripts/build_field.py             # 9-plane rank ensemble              -> evidence/field.npy
$PY scripts/build_final_submission.py  # writes the shippable GeoTIFF and its -zeros twin
$PY scripts/run_holdout.py             # per-fold rebuild of the method (the honest refusal)
$PY scripts/sync_site_tables.py --check
$PY -m pytest tests/ -q
