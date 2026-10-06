#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
# Pinned owner mirrors, not a bypass of organizer authentication.
exec .venv/bin/python scripts/download_data.py
