"""Shared test helpers.

The competition inputs are ~420 MB, are not redistributable inside git, and are reproduced on an
unrestricted machine by ``scripts/download_competition_data.sh`` + ``scripts/prepare_data.py``.  CI
therefore runs without them, so any test that needs the real rasters must *skip* rather than fail —
a failing red build for an absent 420 MB file teaches nobody anything.

``require_competition_data`` is that guard.  It is deliberately explicit: a test that silently skips
because a path moved is worse than one that fails.
"""
from __future__ import annotations

from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

REQUIRED = ("training_features.tif", "labels.tif", "sample_submission.tif")


def competition_data_available() -> bool:
    return all((ROOT / "data" / n).exists() for n in REQUIRED)


def require_competition_data() -> None:
    """Skip the calling test unless every required competition raster is present."""
    missing = [n for n in REQUIRED if not (ROOT / "data" / n).exists()]
    if missing:
        pytest.skip("competition data absent (reproduce it with scripts/download_competition_data.sh): "
                    + ", ".join(missing))
