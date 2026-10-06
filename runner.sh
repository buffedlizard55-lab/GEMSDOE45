#!/usr/bin/env bash
# Reproduce every number quoted on the site, in order. ~20 minutes on CPU, no GPU, no network.
set -euo pipefail
PY=${PY:-/home/user/.venv/bin/python}
cd "$(dirname "$0")"
$PY scripts/run_pipeline.py            # metric checks, band screen, catalogue, Kaplan-Meier
$PY scripts/build_field.py             # 9-plane rank ensemble            -> evidence/field.npy
$PY scripts/build_final_submission.py  # writes the shippable GeoTIFF
$PY scripts/run_holdout.py             # per-fold rebuild of the method
$PY -m pytest tests/ -q
