# GEMSDOE45 · fault discovery, with evidence

**Read this README and the standing project prompt below at the start of every session.**

## Download the unique research submission

**[Download the GeoTIFF](docs/downloads/gems45-phase-parity-20261006-362e1338f293-research.tif)** · [ZIP](docs/downloads/gems45-phase-parity-20261006-362e1338f293-research.zip) · [Format receipt](docs/downloads/submission-audit.json)

**[Website / executive summary](https://buffedlizard55-lab.github.io/GEMSDOE45/)** · [Submission instructions](https://buffedlizard55-lab.github.io/GEMSDOE45/executive-summary.html)

- Submission name: `GEMS45-PHASE-362e1338f293`
- Note: `H45 differential phase parity; 46 channels; nested geographic calibration; all-fault target. Research only; off-catalogue gate not passed.`
- **Research only. Do not spend a weekly submission slot yet.** Local format passes; portal acceptance and hidden-fault performance are unverified.
- New model, not a copied or renamed submission. 82,678 emitted pixels; 116,332 pixels differ from the supplied GEMSDOE32 artifact. Pixel comparison against all 309 compatible prior raster blobs in our inventory found no identical predictions (two additional blobs had incompatible grids/bands).
- Four-fold spatial proxy DTI: baseline **0.193084**, phase model **0.208521**, mean improvement **+0.015437**, positive in **4/4** folds. These are public-catalogue transfer scores, not leaderboard scores. No comparable inherited holdout existed in this initially empty repository.
- Official leaderboard snapshot on 2026-10-06: **0.3345** at rank 1, not the brief's stale 0.3195. A public score of 0.2778 exists; its mapping to the named GEMSDOE32 artifact is user-reported.

## What is implemented

A CPU-runnable, deterministic, spatially blocked 31-versus-46-channel feature ablation, all-fault inference, exact published DW-Tversky formula, safe GeoTIFF writer, internal footprint mask, single-file ZIP, pixel-level prior comparison, tested download-first GitHub Pages site, daily source-health refresh, and claim/source registers. The pipeline uses gradient boosting rather than pretending that a pre-existing GPU pipeline was present.

Read [scientific analysis](research/analysis.md), [frozen hypotheses](research/hypotheses.md), [holdout results](evidence/holdout.json), [data provenance](evidence/data_receipt.json), [official-source register](evidence/official_sources.json), [prior-method inventory](evidence/prior_inventory.json), [all user-reported score observations](research/score_claims.csv), and [three-pass review](research/review.md).

## Reproduce without competition credentials

Python 3.11, approximately 4 GB RAM and 5 GB free disk. No GPU required for this model. `gh` must be available and connected to GitHub; do not put credentials in files or chat.

```bash
python -m venv .venv
.venv/bin/pip install -r requirements.txt
bash scripts/download_competition_data.sh
.venv/bin/python scripts/prepare_data.py
.venv/bin/python scripts/run_experiment.py
.venv/bin/python scripts/build_site.py
.venv/bin/pytest -q
.venv/bin/python scripts/verify_repo.py
```

The download script restores **checksum-pinned owner mirrors**, not authenticated official bytes. Raw data, feature arrays and scratch files are ignored. `sample_submission.tif` is used only for grid and valid footprint, never as predictions. Its positive values are an explicit provenance irregularity. All train/cal/test decisions are frozen in `research/hypotheses.md` (SHA recorded with results). Read that document before changing anything; do not retroactively edit it to fit results.

To refresh prior source and pixel inventories (network-heavy, not required to train):

```bash
.venv/bin/python scripts/review_prior_sites.py
.venv/bin/python scripts/fingerprint_priors.py
.venv/bin/python scripts/build_site.py
```

Source-health checks: `.venv/bin/python scripts/refresh_sources.py`. Daily GitHub Actions refreshes the deployed feed without falsely changing the dated scientific evidence. A failed check displays failure/staleness, never invents a current score. Ordinary sandbox HTTPS failed here while research-tool access and `gh` worked. GitHub Pages deployment status must be checked independently.

## Next session — ordered work

1. **Read this brief, the irregularities and the measured results first.** Do not reuse the four test outcomes for more tuning and call the same holdout fresh.
2. Obtain a genuinely independent, legally usable off-catalogue evaluation set and the incumbent's original reproducible holdout protocol. Check overlap with the supplied catalogue before scoring. Published external data can help, but merely being newer does not prove independence. Never request hidden competition test labels.
3. Resolve ambiguous `tc` and depth-band metadata against the authenticated organizer documentation. Authenticate mirror bytes through an authorized data connection or organizer-supplied checksums.
4. Investigate high-budget boundary selection: all inner calibration folds chose the largest preregistered budget (1.6%). This is a diagnostic, not permission to tune upward on the consumed outer holdout. Register a new calibration/confirmation split first.
5. Evaluate higher-resolution terrain evidence, geological hard negatives (river/lake scarps, lithologic contacts, roads) and fault-segment-level withholding. The official USGS GeoDAWN and GDR resources are research leads; verify download, coverage, licensing and independence before treating a new source as viable.
6. Only promote a weekly slot after comparable incumbent improvement and off-catalogue confirmation. No automated competition upload or eligibility certification is implemented.

## Core values

**Maximize P(Win).** Weigh evidence, risk and tradeoffs. Reserve limited submissions for validated experiments rather than score-chasing.

**Own the Outcome.** Report failed experiments, blockers and uncertainty. Own the complete data→validation→download→deployment chain. Do not hide a negative result behind a visually attractive site.

---

# Standing project prompt

The following preserves the user's project instructions in normalized formatting. Duplicate verification paragraphs and malformed nested links are consolidated, not presented as a verbatim transcript. The complete supplied score observations are preserved in `research/score_claims.csv`; prior URLs and pinned source revisions are in `evidence/prior_inventory.json`.

> Review the repo. Work on next steps from previous sessions first. The highest urgency is to **generate a unique TIF submission** for the competition. Do not copy a previous submission except for learning/education. The submission must be different from the listed GEMSDOE projects. Read the entire prompt. The goal is to place at the top of the leaderboard and ultimately contribute to scientifically grounded geothermal discovery.
>
> The earlier GEMSDOE sites are tested starting points, not files to reuse. Study the best supplied result: GEMSDOE32 `h33-h33-2-b2-20261004T220000Z-e5eb6e7e-zeros`, reported DTI **0.2778**. Explain why and how it performed best, whether a higher score is achievable, and use that analysis to generate a new TIF. The prompt also cites **0.3195** as the leading score; verify current leaderboard information rather than assuming it remains current.
>
> Before implementing, generate **3–5 candidate geological hypotheses** not already tried. Each must name layers, physical signature/transform, why it might catch a fault missing from the USGS/INGENIOUS catalogue, how it differs from prior implementation, expected DTI improvement and implementation cost. Rank the candidates. Validate the top candidate on a **spatially blocked holdout** before touching a weekly submission slot. Do not spend a slot on an idea that has not beaten the current holdout best. If new external data are needed, name a free official source and confirm obtainability before treating the idea as viable.
>
> Work autonomously without asking for manual data placement. Research official, trusted scientific sources; organize an auditable table of knowledge, calculations, methods, hypotheses, data links and evidence. Be contrarian but scientifically grounded. Distinguish observations from hypotheses, flag irregularities, and do not hallucinate. Verify code, calculations and claims rather than claiming unsupported line-by-line certainty. Explain limitations and needed access. Store knowledge for future projects and keep a current source feed.
>
> Put this prompt in the README and read it every session. Keep **Maximize P(Win)** and **Own the Outcome** central: weigh tradeoffs and risk, make evidence-led decisions, take accountability end to end, act on problems without waiting for permission, and learn from both failures and successes.
>
> Build a clean, simple, user-friendly **GitHub Pages website**. At the very beginning/executive summary, make the submission TIF download obvious. Include relevant methods, research, results and official source links for manual review. Create an **executive-summary subpage** explaining exactly how to submit, with a unique submission name and a short optional note. The portal accepts one single-band GeoTIFF or a ZIP containing one GeoTIFF, matching the CRS, shape and geotransform. Fix the previously encountered error: **“Predicted values must be in range [0, 1]”.**
>
> Follow the official competition overview, problem description, resources and rules PDF. Collect data, create/train a model, generate appropriately formatted predictions, and provide reproducible solution assets. Investigate free official external data when useful. The earlier-session claim was: “The single remaining blocker to training is data placement: run `bash scripts/download_competition_data.sh` into `data/`, then `python scripts/prepare_data.py`; the full train→inference→validate pipeline is ready (GPU needed for training).” Check whether that claim applies to this repository; complete the necessary work independently, not by repeating instructions to the user.
>
> Run **three passes**: (1) implement and verify; (2) review bugs, missing requirements, assumptions and edge cases, fixing problems; (3) recheck against the original request and improve accuracy, reliability, completeness and quality. Do not stop after pass one. Be explicit about unmet requirements rather than asserting full completion without evidence.
>
> Create a pull request and merge it to main. Suggest remaining work and limitations for the next sessions. Continue emphasizing usable, unique, verified competition artifacts rather than documentation alone.

### User-supplied official/reference links

- [Competition](https://www.drivendata.org/competitions/306/competition-doe-gems/)
- [Problem and format](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/)
- [About / resources](https://www.drivendata.org/competitions/306/competition-doe-gems/page/968/)
- [Leaderboard](https://www.drivendata.org/competitions/306/competition-doe-gems/leaderboard/)
- [Authenticated data page](https://www.drivendata.org/competitions/306/competition-doe-gems/data/)
- [DrivenData reference solution](https://github.com/drivendataorg/gems-prize-reference-solution)
- [Official September 2026 rules PDF](https://docs.nlr.gov/docs/fy26osti/96647.pdf)
- [GDR INGENIOUS compilation](https://gdr.openei.org/submissions/1391)
- [GEMSDOE32 best supplied artifact site](https://buffedlizard55-lab.github.io/GEMSDOE32/docs/index.html)

### User-supplied data mirrors (not independently organizer-authenticated)

- [GEMS_96647.pdf](https://www.dropbox.com/scl/fi/aemhtutjgcp6tr3tint94/GEMS_96647.pdf?rlkey=rek210cj2smnmzb8n0sla1vmd&dl=1)
- [example_submission.tif](https://www.dropbox.com/scl/fi/6rgvnuady818ol8yqgis4/example_submission.tif?rlkey=kbykilvau066xuogoosbf4cq8&dl=1)
- [existing_faults.tif](https://www.dropbox.com/scl/fi/t7fyt03qdh9egyme0itwo/existing_faults.tif?rlkey=yiao96uluqdkipf0h5vju71jf&dl=1)
- [gems-geodawn-numerical-features.tif](https://www.dropbox.com/scl/fi/3vz9o0wwavi26xaeoxlwr/gems-geodawn-numerical-features.tif?rlkey=je8d8fepqfbst9lnwsq9rkplu&dl=1)
- [Digital-elevation-model-links-JSON.pdf](https://www.dropbox.com/scl/fi/ig0mban712ns1atphgphe/Digital-elevation-model-links-JSON.pdf?rlkey=zm77f1vbtt2if8hlruymptnu3&dl=1)

### Prior project list from the prompt

GEMSDOE, GEMSDOE2, GEMSDOE3, GEMSDOE4, 5GEMSDOE, 6GEMSDOE, 7GEMSDOE, 8GEMSDOE, GEMSDOE9, GEMSDOE10, 11GEMSDOE, 12GEMSDOE, 13GEMSDOE, 14GEMSDOE, 15GEMSDOE, 16GEMSDOE, 17GEMSDOE, 18GEMSDOE, 19GEMSDOE, 20GEMSDOE, GEMSDOE21 through GEMSDOE41, and the names 42GEMSDOE, 43GEMSDOE, 44GEMSDOE (resolved to GEMSDOE42, GEMSDOE43, GEMSDOE44). Pages under `/docs/`, root, and `index.html` are represented by their checked repository landing source in the inventory. No blank score is interpreted as zero.
