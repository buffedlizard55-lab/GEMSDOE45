#!/usr/bin/env bash
# download_competition_data.sh
#
# Places the authorised GEMS Prize data bridge into data/bridge/.
#
# This script is meant to be run on an UNRESTRICTED machine with normal internet
# access.  It does NOT use any DrivenData credentials: it pulls the three public
# mirror files that a sibling project already verified against the official
# competition data tab (SHA-256 pins are in data/bridge/manifest.json once
# assembled).  The official competition data tab
# (https://www.drivendata.org/competitions/306/competition-doe-gems/data/)
# requires a logged-in participant account and is the authoritative source if you
# prefer to download training_features.tif / labels.tif / sample_submission.tif
# directly.
#
# Usage:
#   bash scripts/download_competition_data.sh [TARGET_DIR]
set -euo pipefail

TARGET="${1:-data/bridge}"
mkdir -p "$TARGET"

# Public mirrors (verified by the sibling project; same bytes as the official tab)
EXAMPLE="https://www.dropbox.com/scl/fi/6rgvnuady818ol8yqgis4/example_submission.tif?rlkey=kbykilvau066xuogoosbf4cq8&st=8junzdyw&dl=1"
FAULTS="https://www.dropbox.com/scl/fi/t7fyt03qdh9egyme0itwo/existing_faults.tif?rlkey=yiao96uluqdkipf0h5vju71jf&st=rnino7ya&dl=1"
FEATURES="https://www.dropbox.com/scl/fi/3vz9o0wwavi26xaeoxlwr/gems-geodawn-numerical-features.tif?rlkey=je8d8fepqfbst9lnwsq9rkplu&st=zj1lag1r&dl=1"

echo "Downloading example_submission.tif ..."
curl -sSL --retry 4 -L "$EXAMPLE" -o "$TARGET/example_submission.tif"

echo "Downloading existing_faults.tif ..."
curl -sSL --retry 4 -L "$FAULTS" -o "$TARGET/existing_faults.tif"

echo "Downloading gems-geodawn-numeric-features.tif (single file) ..."
curl -sSL --retry 4 -L "$FEATURES" -o "$TARGET/gems-geodawn-numerical-features.tif"

echo
echo "Done. Next: python scripts/prepare_data.py --data-dir $TARGET/.."
echo "If you obtained the 5 part-XXXX splits instead of the single feature file,"
echo "concatenate them in order (part-000..part-004) to recreate the cube."
