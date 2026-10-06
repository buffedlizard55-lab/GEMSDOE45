# GEMSDOE45 · two workstreams, one repository

Both projects live here. Each has its own site, its own submission GeoTIFF, and its own evidence.
**No score on this page is an organizer receipt** — see the provenance notes inside each section.

| | Download the submission | Site | Status |
|---|---|---|---|
| **H51 — Kaplan–Meier fault-zone emission** (newest) | **[GeoTIFF](docs/downloads/gems45-h51-km-faultzone-20261006-zeros.tif)** · [NaN-outside](docs/downloads/gems45-h51-km-faultzone-20261006.tif) · [.zip](docs/downloads/gems45-h51-km-faultzone-20261006-zeros.zip) | [site](docs/h51/index.html) · [how to submit](docs/h51/executive-summary.html) | **UNSCORED CANDIDATE**, 37,654 dots, 30 distinct realised extension lengths, 309 priors compared / 0 identical |
| **Kaplan–Meier fault-tip survival** (H48b) | **[GeoTIFF](docs/downloads/gems45-h48b-structural-noflank-20261006.tif)** · [.zip](docs/downloads/gems45-h48b-structural-noflank-20261006-zeros.zip) | [site](docs/e45/index.html) · [how to submit](docs/e45/executive-summary.html) | candidate under HOLD, 40,000 dots, unique vs 69 prior rasters (max Jaccard 0.017) |
| **Phase-parity fault detector** | [GeoTIFF](docs/downloads/gems45-phase-parity-20261006-362e1338f293-research.tif) · [.zip](docs/downloads/gems45-phase-parity-20261006-362e1338f293-research.zip) | [site](docs/phase-parity.html) · [summary](docs/executive-summary.html) | research artifact with spatial validation and a uniqueness audit |

Reproduce every published number with **`./runner.sh`** (CPU only, no network, no credentials).
The tables on the Kaplan–Meier site are *generated* from `evidence/*.json`, and
`scripts/sync_site_tables.py --check` fails CI if they drift from it.

---

## H51 — Kaplan–Meier fault-zone emission (this session's submission)

**[Download the GeoTIFF](docs/downloads/gems45-h51-km-faultzone-20261006-zeros.tif)** ·
[all-finite twin](docs/downloads/gems45-h51-km-faultzone-20261006.tif) ·
[.zip](docs/downloads/gems45-h51-km-faultzone-20261006-zeros.zip) ·
[site](docs/h51/index.html) · [how to submit](docs/h51/executive-summary.html) ·
[sources](docs/h51/sources.html)

- Unique submission name: `gems45-h51-km-faultzone-20261006`.
- **What is genuinely new here.** The extension beyond a mapped fault tip is drawn **per fault** by
  inverse transform from a **Kaplan–Meier product-limit estimator** (Kaplan & Meier 1958) fitted to
  the catalogue's own along-strike relay gaps — 3,163 events /
  703 right-censored, median gap 9 px, restricted
  mean 37.0 px, S(3 px) = 0.667. It is
  a distribution, not a radius: **118 distinct drawn lengths** and
  **30 distinct realised extensions** (3,636
  of 3,866 tips carry at least one dot). The build and `tests/test_h51.py` both assert
  this spread, so a constant-length collapse fails instead of shipping.
- **Candidate set** = the organizers' own definition of a new fault ("any fault pixel not already
  captured by USGS/INGENIOUS", "can include newly mapped geometry of an existing fault system"):
  115,554 along-strike tip continuations + 236,798 splay /
  parallel-strand positions at ±2 and ±3 px; 153,447 survive the corroboration and
  catalogue-novelty gates; 28,853 tip + 8,801 strand dots
  are emitted.
- **Selection is the metric's own decision problem.** Because `TP_w` is a MAX over predictions near
  each truth pixel and `FP_w` is a SUM over emitted mass, the optimal prediction is binary, and a
  dot pays for itself exactly when its expected kernel weight exceeds `0.2·DTI`. The file is
  therefore selected by greedy maximum expected coverage — the submodular (1 − 1/e) approximation of
  the metric's own objective — stopped at the marginal bar evaluated at the portfolio's best
  live-observed DTI (0.2778, bar 0.0556).
- **Calibrated on live data, not on the catalogue.** The 2 px (200 m) catalogue exclusion is the
  radius the family's live leaderboard actually validated. The three nested files
  `dotted-h19-5-d2-8` (44,090 px → 0.2600) ⊃ `h27-4-solo-d28` (40,199 px → 0.2708) ⊃
  `h33-2-b2` (37,654 px → 0.2778) were verified by **exact integer-coordinate set inclusion**, so
  the live score rises monotonically as mass is removed. H51 adopts that radius and that mass
  ceiling (37,654 dots).
- **Format, re-read from disk.** All-finite twin: 142,899 B, sha256
  `81c90dedb7f775aa6e1731088be90e2f9f63a7ac0a9c979638abdf5a842e2a0c`, 37,654 positive px, min 0.0, max
  1.0, no NaN inside, no sentinel anywhere, single band float32, EPSG:32611,
  3730 × 3292. NaN-outside twin: 185,873 B, sha256 `b8d0bdc55b8ef20ba7ccf5a2d43773092663f4477967891ab1aab9327c80a69c`. ZIP: 108,324 B,
  sha256 `5110361c2b2f3a6ba2015bd42ce06c790c5c6d23e893c02de867e73add22d989`.
- **Uniqueness.** 309 fingerprinted prior rasters compared by
  sha256 of the float32 bytes: **0 identical**. Largest Jaccard against any member of the family's
  live ladder is 0.0074.
- **What is NOT validated, stated plainly.** No instrument available here ranks this file against the
  hidden new-fault truth. The catalogue-truth holdout is inverted by construction (IR-45-003); the
  emission emulator that placed the dots correlates with live scores at only
  **Pearson +0.176** (R² +0.031) over 46 live-scored rasters — it
  predicts 0.575 where the live top is 0.3345, and it cannot even order `h33-2-b2` (live 0.2778)
  above `h34-scatter` (live 0.0778). `slot_approved` is **false**. Read
  [the limitations](docs/h51/index.html) before spending a weekly slot.
- **Reproduce:** `.venv/bin/python scripts/build_h51_submission.py` (≈2 min with the tip cache;
  ~18 min cold) then `scripts/validate_h51.py`.

---

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


---

# GEMSDOE45 — Kaplan–Meier fault-tip survival, and an honest audit of the GEMS Prize metric

**Competition:** [The Geologic Enhanced Mapping System (GEMS) Prize Challenge](https://www.drivendata.org/competitions/306/competition-doe-gems/)
· **Region:** GeoDAWN, northwestern Great Basin, Nevada · **Grid:** EPSG:32611, 100 m, 3730 × 3292

---

## ⬇️ Download the submission

| File | Purpose | Format |
|---|---|---|
| **[`docs/downloads/gems45-h48-kmtip-structural-20261006.tif`](docs/downloads/gems45-h48-kmtip-structural-20261006.tif)** | **Primary submission** — matches the organizers' own template encoding (NaN outside) | single-band float32 GeoTIFF |
| `docs/downloads/gems45-h48-kmtip-structural-20261006-zeros.tif` | Belt-and-braces twin — **every** cell finite and in [0, 1], no NaN anywhere, for a validator that rejects NaN | single-band float32 GeoTIFF |
| `docs/downloads/gems45-h48-kmtip-structural-20261006-zeros.zip` | Same file zipped, for the portal's `.zip` option | zip containing one GeoTIFF |

**Suggested submission note** (paste into the *Note (optional)* box):

```
GEMSDOE45 H48 | KM tip-survival extension (3866 catalogue tips, product-limit gap curve,
per-tip lengths 2-11 px) unioned with a 9-plane rank-ensemble structural field; 40,000 dots
at the metric's own marginal bar, all >=200 m from the mapped catalogue
```

**SHA-256 of the primary file** — verify before uploading:

```
40eb32bf995433b15052745c7ea647ed344dc4d56ca9356f8b843a4c1466a402  gems45-h48-kmtip-structural-20261006.tif
```

**Why the portal previously said *"Predicted values must be in range [0, 1]"*.** Two independent
causes are reproduced and fixed here. (1) `training_features.tif` stores the sentinel
**−3.4028235e+38** for cells outside the survey; any raster assembled from it wholesale contains
values near −10³⁸. (2) `0 <= nan <= 1` is **False**, so a validator that range-tests the whole array
rejects a NaN-outside file even though the organizers' own `sample_submission.tif` is NaN-outside.
Every file published here is asserted before writing: single band, float32, template CRS/transform,
zero sentinel values anywhere, all finite in-footprint values in [0, 1]. Both encodings are shipped
so the choice is yours, not the writer's.

---

## Executive summary

1. **The metric is verified, not assumed.** `src/gems45/metric.py` is a transcription of the
   published formulas and reproduces the organizers' worked example to 0.602652 (page states 0.60).
   A literal O(N_truth × N_pred) brute-force transcription is checked against it on random rasters
   in `tests/test_metric.py`.

2. **The single most decision-relevant rule is a staff statement, not a page of documentation.**
   DrivenData staff confirmed on the forum that *"pixels corresponding to known USGS/INGENIOUS
   faults are masked / excluded from evaluation, so they do not count towards penalty terms"* and
   that the Final Round re-score masks them too. Every design choice here follows from that:
   mass on the catalogue is worth nothing, so the whole problem is placing mass **off**-catalogue.

3. **The metric reduces to one design equation and one emission rule.**
   With `T = TP_w`, `F = FP_w`, `K = |G|`:

   ```
   1/DTI = 0.2 + 0.2·(F/T) + 0.8·(K/T)          DTI = T / (0.2(T+F) + 0.8K)
   ```

   and adding one unit of mass at a pixel whose realised kernel weight is `w` raises the score
   **iff `w > 0.2·DTI`** — exactly, because the denominator rises by exactly `α` regardless of
   where the mass lands. At DTI = 0.30 the bar is 0.06: a dot must sit within ≈282 m of a *hidden*
   truth pixel. This is the competition's whole emission rule, and it is derived, not fitted.

4. **Nine of the nineteen supplied bands cannot localise a fault at the kernel's scale, and one of
   them is broken.** Measured (see the site's band table): the three geodetic strain-rate bands and
   `ieq`/`deq` are regional fields with ρ(1 km) ≥ 0.98 — a strain-rate *ridge detector* physically
   cannot satisfy a 300 m kernel, which is why that family of ideas under-performs. `tc` (band 6)
   additionally carries 12 nodata sentinels **inside** the nominal footprint. The usable
   fault-scale information lives in the topographic-slope and gravity-derivative gradients
   (best single-feature AUC 0.61; `det_elev_slope` is the strongest at 2.8× background enrichment).

5. **The novel contribution: tip extension as a survival function, not a search radius.**
   The supplied catalogue is not long traces but **3,199 8-connected segments with a median length
   of 12 pixels**. Successive sections of one structure are arranged en echelon along a trend,
   so for every one of the **3,866 segment tips** we measure the along-strike distance to the
   nearest pixel of a *different* segment — an event for 3,163 tips, right-censored for 703 — and
   fit the product-limit estimator of Kaplan & Meier (1958). Each tip is then extended by its *own*
   length taken from that curve: measured per-tip lengths span **2–11 px (10 distinct values,
   mean 7.24 px, 11.4 % of tips at the cap)**, so the extension demonstrably varies by fault rather
   than collapsing to a constant.

6. **What is shipped, and what is refused.** 40,000 dots: the Kaplan–Meier tip-extension and relay
   bridges, unioned with the top of a 9-plane rank ensemble of the measured-informative bands,
   every dot at least 200 m from the mapped catalogue. **No score is claimed.** On the only
   spatially-blocked instrument this environment can run (catalogue truth) the file scores 0.0068
   against a random control's 0.0127 — it does **not** beat the incumbent, and by this repository's
   own promotion rule the submission is therefore published as a **candidate under HOLD**, not as a
   validated improvement. Section 4 of the site explains, with numbers, why that instrument cannot
   settle the question either way.

---

## THE PROJECT BRIEF — read at the start of every session

> The text below is the working brief, reproduced verbatim so that every session starts from the
> same statement of what is being aimed for. Core Values to hold throughout: **Maximize P(Win)**
> and **Own the Outcome**.

**Highest urgency.** Generate a *unique* TIF submission. Do not copy a previous submission except
for learning and education — the submission must differ from every artifact in the collection of
GEMSDOE sites listed in the brief. There must be an easy-to-download submission TIF. Read the
entire brief.

**The scientific instruction.** Model how far a fault plausibly extends past its mapped tip as a
*survival function*, not a fixed search radius. The tip-and-step-over candidates
(prethin-tip-euler 0.2649, analog-tip-stepover 0.2632) are already the second-strongest family, but
they are built on a guessed extension distance. Kaplan & Meier's nonparametric survival estimator
(*Journal of the American Statistical Association*, 1958) — built for exactly the "how far does
something go before it stops" question — gives a disciplined alternative: treat the distance from
each mapped fault's last point to where it actually appears to terminate (geomorphically, or by
proximity to where other faults in the catalogue terminate) as a time-to-event sample, fit the
empirical survival curve from the catalogue's own tip distances, and extend every known fault's tip
by a distance drawn from that fitted curve rather than one fixed guess — shorter, well-attested
faults get a shorter plausible extension, longer or more structurally active ones get more
latitude, both calibrated from your own data instead of asserted. Normalize to [0, 1], write to the
required format, and confirm the extension lengths actually vary by fault rather than collapsing
back to a constant before download.

**The record to beat.** Owner- or snapshot-reported scores across the family record include
0.1563, 0.0286, 0.1193, 0.0830, 0.1152, 0.1560, 0.0343, 0.1461, 0.0107, 0.0202, 0.1294, 0.0782,
0.0020, 0.0187, 0.0297, 0.1894, 0.1922, 0.0461, 0.0921, 0.1280, 0.1839, 0.0904, 0.1855, 0.0976,
0.0360, 0.1890, 0.1859, 0.1002, 0.0748, 0.1352, 0.2477, 0.2600, 0.1223, 0.2449, 0.2708, 0.2649,
0.2600, 0.0041, 0.2778, 0.2632, 0.0778, 0.0418. The best single site result is
`h33-h33-2-b2-20261004T220000Z-e5eb6e7e` at **0.2778**. The brief states 0.3195 as the leaderboard
top; the GEMSDOE32 snapshot of 2026-10-04 records **0.3262** at rank 1 and 0.3195 at rank 3.
Both are unreceipted here — see the irregularities ledger.

**Method requirements.**
* Study the highest-scoring artifact and explain *why* it scored what it scored, with numbers.
* Generate 3–5 candidate geological hypotheses not yet tried, each naming the specific layer(s)
  involved, the physical signature being targeted, why it should catch a fault missing from the
  USGS/INGENIOUS catalogue rather than one already in it, and how it differs from anything already
  implemented. Rank by expected DTI improvement and implementation cost.
* Validate the top candidate on a **spatially-blocked holdout before touching a weekly submission
  slot**; do not spend a slot on an idea that has not beaten the current holdout best. If a
  candidate cannot be validated without new external data, name the specific free official source
  and check it is obtainable before proposing it as viable.
* Work line by line, verify against official and trusted sources, provide links for manual review,
  flag irregularities, and **no hallucinations**. Where a number is a claim rather than a receipt,
  label it as one.

**Deliverables.**
* A unique, format-valid submission TIF, downloadable in one click, with its download obvious at
  the very top of the site.
* An executive-summary subpage explaining exactly how to submit to the competition.
* A clean, user-friendly GitHub Pages site containing all relevant information with official
  verified source links, so nobody has to check anything by hand.
* A project that keeps itself current rather than needing manual re-checking.

---

## What was actually measured here (all first-hand, all reproducible)

| Quantity | Measured value | Where |
|---|---|---|
| Official metric worked example | 0.602652 (page: 0.60) | `tests/test_metric.py` |
| Brute-force vs reduced metric | agree < 1e-5 on random rasters | `tests/test_metric.py` |
| Grid | EPSG:32611, 100 m, 3730 × 3292, origin (243350, 4508550) | `evidence/pipeline.json` |
| `training_features.tif` | 19 bands, float32, nodata −3.4028235e38, sha256 `4371c82e…23bc5` | `evidence/pipeline.json` |
| `labels.tif` | int8 {−1,0,1}, 60,988 fault px, sha256 `7ba308cc…4093` | `evidence/pipeline.json` |
| `existing_faults.tif` | byte-identical to `labels.tif` | sha256 comparison |
| `sample_submission.tif` | float32, nodata NaN, sha256 `2176d08e…5cbc` | `evidence/pipeline.json` |
| Footprint **disagreement** | band-1 finite 5,165,852 · all-band finite 5,165,840 · template finite 5,167,373 | IR-45-001 |
| Band 6 `tc` | 12 nodata sentinels **inside** the band-1 finite mask | IR-45-002 |
| Catalogue structure | 3,199 8-connected segments, median 12 px, largest 360 px, 3,866 tips | `evidence/pipeline.json` |
| KM sample | 3,866 observations, 3,163 events, 703 right-censored at 120 px | `evidence/pipeline.json` |
| KM gap curve (genuine, gap > 2 px) | n = 2,127 · median 12 px · q75 39 px · restricted mean 26.9 px · Greenwood SE ≤ 0.0108 | `evidence/final_submission.json` |
| Extension budget from the metric | e* = 11 px (1,100 m) at the 0.06 bar for DTI 0.30 | `evidence/final_submission.json` |
| Per-tip extension lengths | min 2 · max 11 · mean 7.24 · median 6 · **10 distinct values** | `evidence/final_submission.json` |
| Shipped file | 40,000 positive px, 138,260 B, sha256 `40eb32bf…a402` | `evidence/final_submission.json` |
| Shipped file vs catalogue | min distance 2.0 px · median 11.4 px · 11.9 % inside 300 m | `evidence/final_submission.json` |
| Blocked holdout (catalogue truth) | **H48b 0.030215** · H48 0.006813 · random 0.012654 · d2.8 0.0634 · h33-2-b2 0.002791 | `evidence/holdout_h48b.json`, `evidence/holdout_h48.json` |
| Blocked holdout, blocks won vs random | **H48b 9/11** · H48 2/11 · H48b beats d2.8 0/11 | same |
| Proximity screen, measured random baseline | **8.60 %** of footprint is within 300 m of the catalogue | `evidence/band_screen.json` |
| Best raw band by top-40k within 300 m | `geod_2ndinv` 30.6 % — but its top 40,000 px are **one connected component** (perimeter/area 0.02): regional, not a localiser | `evidence/band_screen.json` |
| Best derived plane | `\|grad2.5\|_det_elev_slope` **42.1 %** = 4.89× random, 1,162 components, median 13 px — a genuine localiser | `evidence/band_screen.json` |
| Band value that looked best and is not | `det_elev_slope` scores **5.1 %**, *below* the 8.6 % baseline; its gradient is what scores | `evidence/band_screen.json` |
| H47 consensus field, same instrument | 9.9 % = 1.15× — rejected | `evidence/field_screen.json` |
| H48 rank ensemble, same instrument | top-20k 27.6 % (3.21×) · top-40k 24.4 % (2.83×) · top-80k 21.2 % | `evidence/field_screen.json` |
| Second published file | H48b 138,304 B, sha256 `be7915f2…2ea8c`, max Jaccard vs 69 priors **0.017** | `evidence/final_submission_noflank.json` |
| Reproducibility of the shipped bytes | rebuilding with the default flags reproduces `40eb32bf…a402` **exactly** | `scripts/build_final_submission.py` |

Reproduce everything — or just run **`./runner.sh`**, which runs the whole chain in order:

```bash
./runner.sh
# equivalently, step by step:
python scripts/run_pipeline.py            # metric checks, band values, catalogue, Kaplan-Meier
python scripts/screen_bands.py            # per-band / per-plane localisation + fragmentation
python scripts/score_fields.py            # every stored field on ONE instrument
python scripts/sync_site_tables.py        # regenerate the site tables from the evidence
python scripts/build_field.py             # the 9-plane rank ensemble
python scripts/build_final_submission.py  # writes the shippable GeoTIFF and its -zeros twin
python scripts/run_holdout.py --pred <file>   # leakage-controlled blocked holdout
python scripts/sync_site_tables.py --check    # fails if the site has drifted from the evidence
python -m pytest tests/ -q
```

`scripts/sync_site_tables.py --check`, `scripts/verify_repo.py` and the test suite both run in CI
terms; the test suite includes `tests/test_site_integrity.py`, which fails if any published sha256
stops matching the bytes it describes.

---

## Ranked hypothesis set (the brief's required 3–5), and their verdicts

Ranked by expected DTI gain per unit of implementation cost. "Tried before?" is checked against the
cloned sibling repositories (GEMSDOE32, 40, 41) and the family record, not asserted.

| # | Hypothesis | Layers / physical signature | Why it should find a *new* fault | Tried before? | Verdict |
|---|---|---|---|---|---|
| **H46** | **Kaplan–Meier tip-survival extension** — fit the product-limit curve of along-strike relay gaps and extend each tip by its own draw | catalogue geometry only; `det_elev_slope` (19) as an optional veto | an expert adds a fault by *extending* a known one or filling a relay gap far more often than by inventing a trace in blank ground; the extension zone is a 1-D uncertainty, not a 2-D search | **No** — every prior tip family used a fixed radius (d2.8, prethin-tip-euler, analog-tip-stepover) | **Implemented & shipped.** Lengths vary 2–11 px; unvalidated |
| **H47** | Multi-physics oriented-lineament consensus (second-order orientation order parameter across independent fields) | `tmi_hg`(3) `tmi_vg`(9) `rtp`(2) `tmi`(14) `det_elev_slope`(19) `iso_grav_anom_hg`(18) `iso_grav_anom_vg`(11) `tc`(6) | the catalogue over-represents faults mappable in *one* dataset; demanding that several independent physics agree on a strike selects what single-dataset mapping misses | Partly — GEMSDOE32's H33-A proposed the same object | **Tested and rejected.** Measured 9.9 % of top-40k within 300 m of truth vs 8.3 % random — a 1.2× enrichment, *worse* than the single best plane (2.8×). The geometric mean over 16 terms is too brittle |
| **H48** | Rank ensemble of the measured-informative planes (no fitting on labels) | 9 planes: gradients of `det_elev_slope`, `det_elev`, `iso_grav_anom{,_hg,_slope}`, `geod_dilaterate`, `val_geod_2ndinv`, `|laplace|` | ranking by percentile makes the combination robust to the heavy tails that dominate potential-field derivatives, and to the wild dynamic-range mismatch between magnetic and gravity planes | **No** — prior repos use per-band products or supervised scores | **Implemented & shipped** as the structural support. 3.3× background enrichment at top-20k |
| **H49** | Cross-catalogue disagreement field, gated by structural coherence | USGS SGMC (public domain) faults absent from the given catalogue, weighted by the H48 ensemble | represents *known-real* faults the provided catalogue lacks, which is the hidden-truth population's character | Partly — GEMSDOE29 shipped an unscored `sgmc-off-catalogue-44k`; GEMSDOE32 ranked it direction 2 | **Blocked.** `mrdata.usgs.gov` is unreachable from this sandbox (HTTP 000, curl). Source named, obtainability **not** verified, so not proposed as ready |
| **H50** | Sub-kilometre seismicity lineaments (epicentre density + focal-mechanism nodal strikes) | replaces `ieq_n100a15`(16) `deq_n100a15`(10) with USGS FDSN/ComCat events | instrumental seismicity is a *dynamic* inventory: a blind, low-slip, active fault lights up seismically while leaving no scarp | **No** — nobody in the record has fetched a raw catalogue | **Blocked, and the supplied proxies are provably inadequate:** `ieq` has ρ(1 km) = 0.9939 and `deq` 0.9447, i.e. `ieq` is constant by construction at the metric's scale. Free official source named (USGS FDSN event web service); obtainability from this sandbox **unverified** |

**Promotion rule applied.** Per the brief, no submission slot may be spent on a candidate that has
not beaten the current holdout best. On the blocked instrument, the holdout best is the incumbent
`dotted-h19-5-d2-8` at 0.0634. Nothing built here beat it. The honest consequence, stated on the
site rather than buried: **the shipped file is a candidate under HOLD, not a validated win.**

---

## Irregularities flagged for review

| ID | Finding | Evidence |
|---|---|---|
| IR-45-001 | Three different footprints coexist: band-1 finite **5,165,852**, all-band finite **5,165,840**, sample-submission finite **5,167,373**. 1,540 px are all-band finite but NaN in the template; 3,073 px are template-finite but nodata in every feature band. The writer uses the **template** mask (the problem description designates the sample submission as the format template) and this is asserted, not assumed | `evidence/pipeline.json` |
| IR-45-002 | Band 6 `tc` carries **12 nodata sentinels inside** the band-1 finite mask, so "inside the footprint" is not a per-band constant | `evidence/pipeline.json` |
| IR-45-003 | The catalogue-truth holdout **cannot rank candidates**: it puts `d2.8` (live 0.2600) 23× above `h33-2-b2` (live 0.2778), inverting the live order. Any "we beat the holdout" claim built on it is unsound | `evidence/holdout_shipped.json` |
| IR-45-004 | Score provenance. No file anywhere in this repository has an organizer receipt. 0.2778, 0.2600, 0.3262 and every other score in the brief are **owner or snapshot claims**. The live leaderboard is dynamic and was read only through third-party snapshots | site, Sources page |
| IR-45-005 | `deq_n100a15` reaches 4.96 × 10⁶ in a grid spanning ~329 × 373 km, which is not a plausible straight-line distance in metres. The band is used nowhere here, but the value looks like a unit or accumulation error worth reporting to the organizers | `evidence/pipeline.json` |
| IR-45-006 | **Poisson-disk thinning at the kernel radius is wrong for this metric.** A 1-px line at dot spacing 1/2/3/4 px scores DTI 1.000/0.862/0.813/0.712 — thinning *destroys* score, and a marginal dot in a 3-px lattice adds 0.222 credit for 0.044 denominator (ratio 5.0 > DTI), so the good dots are the ones already present. Two independent routines in `emission.py` are now quarantined behind `allow_thinning=False` | `emission.py` docstring, tests |
| IR-45-007 | **A table published on the site was wrong, and the way it was wrong is instructive.** It attributed "23.2 % top-40k within 300 m" to the band `det_elev_slope`; that number belongs to `\|grad1.2\|_det_elev_slope`, and the band's own value scores 5.1 %, *below* the 8.6 % baseline. Its ρ and mean-\|grad\| cells also came from a later pipeline run than its percentage cell, and its class column contradicted the class rule in `bands.py` for three bands. Found by **regenerating** the table from evidence. Fixed by generation + `--check` | `evidence/band_screen.json`, `scripts/sync_site_tables.py` |
| IR-45-008 | **Two bands look like excellent locators and are not.** `geod_2ndinv` (30.6 %) and `geod_shearrate` (27.7 %) post the two best raw-band fractions, yet each one's top 40,000 px form **a single connected component of 40,000 px**. They mark the fault-rich region, not faults. Their gradients change by 0.0023–0.0037 σ/px, so no ridge detector on them can place a dot inside a 300 m kernel. The `perimeter/area` test that exposes this is now published | `evidence/band_screen.json` |
| IR-45-010 | **H51-C is in tension with H51's own novelty gate, and the tension is recorded rather than smoothed over.** The organizers' staff statement says a new-fault ground-truth pixel may lie within 300 m of a known trace as a *correction*; the file therefore emits strands at ±2 and ±3 px of the mapped skeleton — but the live-calibrated 2 px (200 m) catalogue exclusion removes the inner part of that same annulus. Both rules are evidence-backed; together they are inconsistent for the 0–200 m band. The artifact keeps the live-calibrated gate, and the correction annulus is declared partially untested | `evidence/h51_submission.json`, `research/hypotheses_h51.md` |
| IR-45-011 | **The emission emulator is not a score, and saying so is not a disclaimer but a measurement.** Fitted to 46 live-scored rasters it reaches Pearson +0.176 / R² +0.031, predicts 0.575 where the live top is 0.3345, ranks the live-best file below a file scoring 0.0904, and cannot separate `h33-2-b2` (0.2778) from `h34-scatter` (0.0778). Every `dti_pred` value printed by the builder is therefore labelled `emulator_is_not_a_score` and no score claim is derived from it | `evidence/h51_live_fit.json`, `tests/test_h51.py` |
| IR-45-009 | **My own holdout script was scoring a reconstruction and calling it the file.** `scripts/run_holdout.py` loaded the published GeoTIFF and never used it: each fold's `dti_mine` was a per-fold re-derivation, so passing a different file changed nothing (two different files returned an identical 0.000223 — the tell). Fixed: predictions are stripped of the training catalogue before scoring (the metric's own masking rule, applied to the incumbent files too), and `dti_file` / `dti_rebuild` are now separate keys. The previously published 0.006813 was correct and reproduces exactly | `evidence/holdout_h48.json`, `evidence/holdout_h48b.json` |

---

## Limitations — what is blocking a genuinely competitive result

1. **No live scoring loop.** DrivenData forbids automated monitoring and this environment has no
   competition account, so every judgement about the leaderboard is a claim, not a measurement. The
   single highest-value addition would be a manual submission of this file with its note, returning
   one real number.
2. **No instrument that ranks.** The only truth available locally is the catalogue, which the
   organizers mask out. Every sibling repository that has measured this reaches the same
   conclusion. Until an independent fault inventory (SGMC, or a fresh expert set) can be fetched
   and validated, candidate selection is an argument, not a measurement.
3. **The supplied features are weak for this task.** Best single-plane AUC against the catalogue is
   0.61, and 9 of 19 planes carry no information at the kernel's scale. A competitive jump almost
   certainly requires the 1 m DEM (the links CSV is behind the competition login) or a real
   detection model on it — not a better blend of these 19 planes.
4. **No GPU and no training loop.** The reference solution's CNN cannot be trained here; only
   CPU-feasible detectors are represented.
5. **Skeletonisation is slow and single-threaded** (~10 min for the catalogue), and the tip
   extraction over-merges segments that touch. Both are known and documented, not hidden.
6. **The flank question was settled on an instrument with a known bias.** H48b beat H48 4.4× on the
   leakage-controlled blocked holdout, and the blocked protocol is the fairer of the two local
   instruments (held-out truth behaves like unmapped faults, and every file is stripped of the
   training catalogue first). But the gap points in the same direction as the residual bias, so it is
   strong evidence, not proof. A single real leaderboard number on either file would settle it.

---

## Repository layout

```
src/gems45/metric.py        verified DTI, the marginal theorem, the design equation
src/gems45/grid.py          grid geometry, writer with hard format assertions, audit
src/gems45/catalogue.py     segment/tip extraction, local strike estimation
src/gems45/survival.py      Kaplan-Meier product-limit estimator + relay-gap measurement
src/gems45/detector.py      extension budget from the metric's bar, per-tip lengths, gating
src/gems45/emission.py      dot painting, relay bridges, Poisson-disk thinning (with its caveat)
src/gems45/bands.py         band information screen (LOCATOR vs WEIGHT classes)
src/gems45/lineament.py     H47 consensus detector (tested and rejected -- kept as evidence)
src/gems45/detfeatures.py   memory-bounded feature streaming
scripts/                    run_pipeline, screen_bands, score_fields, sync_site_tables,
                            build_field, build_final_submission, run_holdout, sync_evidence_data
tests/                      metric verification against a literal brute-force transcription
evidence/                   every number quoted on the site, as JSON produced by the scripts
docs/                       the GitHub Pages site
data/                       hash-pinned competition inputs (not tracked -- see data/README.md)
```

## Sources for manual review

| Claim | Source |
|---|---|
| Metric formulas, α = 0.2, β = 0.8, R = 300 m, worked example 0.60, submission format | <https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/> |
| Competition structure, prize pools, two-round design, expert-labelled new faults | <https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/#competition-structure> |
| **Known USGS/INGENIOUS faults masked / excluded from evaluation (staff answer)** | <https://community.drivendata.org/t/scoring-clarification-are-known-usgs-ingenious-faults-masked-when-scoring-and-are-they-in-the-final-round-label-set/11516> |
| How the new test faults were identified (question; staff reply not yet retrievable) | <https://community.drivendata.org/t/how-were-the-new-test-faults-identified-data-sources-and-fault-types/11527> |
| Kaplan & Meier (1958), product-limit estimator | <https://doi.org/10.1080/01621459.1958.10501452> |
| GeoDAWN airborne magnetic and radiometric surveys (data DOI) | <https://doi.org/10.5066/P93LGLVQ> |
| INGENIOUS geothermal compilation (GDR submission 1391) | <https://gdr.openei.org/submissions/1391> |
| USGS State Geologic Map Compilation (H49 source, currently unreachable) | <https://mrdata.usgs.gov/geology/state/> |
| USGS FDSN earthquake web service (H50 source, currently unreachable) | <https://earthquake.usgs.gov/fdsnws/event/1/query> |
| Reference solution | <https://github.com/drivendataorg/gems-prize-reference-solution> |
