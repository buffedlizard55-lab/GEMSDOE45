"""Mirror every evidence JSON (and the score-claims CSV) into ``docs/data/``.

``scripts/verify_repo.py`` requires the site's published evidence to be byte-identical to the
evidence on disk, so that a reader can fetch exactly what the analysis produced.  This is the
standalone version of that one step, kept separate from ``scripts/build_site.py`` because that
script regenerates the earlier session's pages wholesale -- useful for them, destructive here.
"""
from __future__ import annotations

import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    out = ROOT / "docs/data"
    out.mkdir(parents=True, exist_ok=True)
    n = 0
    for p in sorted((ROOT / "evidence").glob("*.json")):
        shutil.copyfile(p, out / p.name)
        n += 1
    claims = ROOT / "research/score_claims.csv"
    if claims.exists():
        shutil.copyfile(claims, out / claims.name)
    print(f"mirrored {n} evidence JSON files (+ score_claims.csv) into docs/data/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
