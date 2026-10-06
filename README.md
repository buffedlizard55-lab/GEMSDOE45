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
| Blocked holdout (catalogue truth) | shipped 0.0068 · random 0.0127 · d2.8 0.0634 · h33-2-b2 0.0028 | `evidence/holdout_shipped.json` |

Reproduce everything:

```bash
python scripts/run_pipeline.py            # metric checks, band screen, catalogue, KM  (~14 min)
python scripts/build_field.py             # the 9-plane rank ensemble                    (~35 s)
python scripts/build_final_submission.py  # writes the shippable GeoTIFF                (~10 s)
python scripts/run_holdout.py             # per-fold rebuild of the method               (~80 s)
python -m pytest tests/ -q
```

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
scripts/                    run_pipeline, build_field, build_final_submission, run_holdout, ...
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
