# GEMSDOE45 — unique DOE GEMS geothermal-fault submission

Autonomous, reproducible pipeline that generates a **brand-new, unique** GeoTIFF
submission for the **DOE GEMS Prize Challenge** (DrivenData competition 306:
*Geologic Enhanced Mapping System*). The submission predicts geothermal-fault
likelihood from the 19 supplied geophysical feature bands and is engineered to be
**distinct from every prior GEMSDOE submission** and fully format-compliant
(including the `[0,1]` value range that caused the earlier upload error).

> **Core values (Arena AI):** *Maximize P(Win)* and *Own the Outcome*. We ship a
> verifiable, unique artifact, validate it honestly, and document every limitation.

## What's in this repo

| Path | Purpose |
|------|---------|
| `src/gemsdoe45/metric.py` | Distance-weighted Tversky index (DTI) — exact re-derivation, self-tested |
| `src/gemsdoe45/io_raster.py` | Load the feature cube / template / faults; write **format-validated** GeoTIFF |
| `src/gemsdoe45/features.py` | Geophysical edge / "worm" survival operators |
| `src/gemsdoe45/hypotheses.py` | Candidate generators (H45-1a … H45-4) |
| `src/gemsdoe45/holdout.py` | 4-fold spatially-blocked proxy holdout |
| `src/gemsdoe45/distinctness.py` | Uniqueness audit vs prior submissions (Pearson/Spearman) |
| `scripts/build_submission.py` | End-to-end build → validated TIF + receipts |
| `docs/` | **GitHub Pages site** (landing, executive-summary/submit, methodology, sources) |
| `docs/downloads/` | The submission TIF/zip, preview, and JSON evidence |

## The deliverable (download & submit)

- **Submission file:** `docs/downloads/GEMSDOE45_submission_latest.tif` (mirror of the
  timestamped `gems45-H45-1a-…-nan.tif`).
- **How to submit:** open **`docs/executive_summary.html`** — it has the one-click
  download, the unique submission name, the short comment to paste, and the exact
  step-by-step upload instructions. Also reachable live at the GitHub Pages site.
- **Format:** single float32 band, EPSG:32611, 100 m, 3292×3730, values in `[0,1]`,
  `NaN` outside the GeoDAWN footprint — matches the official sample submission.
- **Uniqueness:** max |Pearson| / |Spearman| correlation with all 24 checked prior
  GEMSDOE submissions = **0.10** (threshold 0.90). It is genuinely new.

## Reproduce / regenerate

The pipeline needs the authorised **data bridge** (the three public mirror files).
Place them and run:

```bash
# 1) Put the bridge files where the loader expects them:
#    data/bridge/example_submission.tif
#    data/bridge/existing_faults.tif
#    data/bridge/gems-geodawn-numerical-features.tif   (re-assemble from the 5 parts)
#    See scripts/download_competition_data.sh for the official mirrors.
python scripts/prepare_data.py --data-dir data        # validate & manifest

python scripts/build_submission.py \
  --features data/bridge/gems-geodawn-numerical-features.tif \
  --template data/bridge/example_submission.tif \
  --faults   data/bridge/existing_faults.tif \
  --priors-dir /path/to/GEMSDOE32/docs/downloads \
  --out-dir  docs/downloads
```

No DrivenData login and no external downloads are required for the build itself;
everything runs from the supplied feature bands.

## Honesty & limitations (read before spending a submission slot)

- **Live score is unknown here.** The hidden test labels are not public, so the
  real DTI cannot be computed. We validate on a *proxy* holdout (the public fault
  catalogue) and on *uniqueness* vs priors.
- **IRR-01 (flagged):** a sibling project measured near-zero correlation
  (Spearman ρ≈0.09, p≈0.80) between catalogue-holdout DTI and the live leaderboard.
  The proxy is a *relative* instrument only — we do **not** claim it predicts the
  official score.
- **Fault-alignment gap:** our edge ensemble tops out near 0.07 proxy DTI on the
  catalogue; the strongest prior methods reach ≈0.18. Closing this (better
  lineament extraction / a learned ensemble, or adding native GeoDAWN radiometrics
  and a real seismicity catalogue — both blocked in this sandbox) is the highest-
  leverage next step to raise the live score. See `registry/irregularities.json`.

## The novel hypothesis (why this can catch unmapped faults)

**H45-1a** leads with the **magnetotelluric depth-to-basement (band 15) and
conductivity-surface (band 17) offset edges** — a deep fault-zone signature that
*no prior hypothesis used* — corroborated by the magnetic+gravity+topographic edge
ensemble and blended with the edge field, activity-gated by geodetic strain and
earthquake density. Full ranking of 3–5 candidate hypotheses, validation, and the
IRR-01 caveat are on the **`docs/methodology.html`** page.

---

## PROJECT PROMPT (verbatim — read at the start of every session)

~~~~
Review the repo.   
  
THE FOLLOWING IS THE HIGHEST URGENCY AND MUST BE FOLLOWED!  
  
MUST GENERATE A UNIQUE TIF SUBMISSION FOR THE COMPETITION.  DO NOT COPY A PREVIOUS SUBMISSION UNLESS IT'S FOR LEARNING AND EDUCATION.  BUT WE MUST GENERATE A UNIQUE TIF SUBMISSION.   The submission must be different than the collection of gemsdoe sites below.  
  
There should be an easy to download submission tif file as described by the prompt.  Read the entire prompt.  
  
The following sites should serve as a starting point for understanding how to generate TIF submissions.  These websites are researched, and tested and have generated TIF submissions.  But we need to generate high scoring submissions.  
  
Here are the results from submissions into the competition, separated by ....:  
  (see the full results list in the repository / project brief — GEMSDOE through GEMSDOE44, with leaderboard scores up to 0.2778 for h33-h33-2-b2 and 0.3195 as the current leaderboard top)
  
We need to study, analyze, and understand the highest score from the GEMDOE site where the submission TIF is downloaded from: h33-h33-2-b2-20261004T220000Z-e5eb6e7e-zeros: 0.2778. Why and how did this get the highest score and are we able to generate a submission that scores higher than 0.2778?  
  
Answer the question using PhD level experience, knowledge, and judgement. Then use the answer to generate a unique TIF submission into the competition.  Must be unique submission unlike any within the GEMSDOE sites above.  Verify working line by line no hallucinations.  
  
The leaderboard for the competition: https://www.drivendata.org/competitions/306/competition-doe-gems/leaderboard/  — top score 0.3195.  Design a new strategy to generate a submission that can score higher than 0.3195.  
  
Put this prompt into the repo readme and read it every time we work on the project as a starting point.  
  
Review the repo.   
  
Core Values: Maximize P(Win); Own the Outcome. Work line by line verifying from official verified trusted sources, provide links for manual review.  There should be no manual input, work on your own to complete tasks.  Flag any irregularities for review.  No hallucinations.  Verify no hallucinations.  
  
We need to focus on being able to generate a submission into the competition.  The site should be able to generate a TIF file that is required for submission.  It should be as easy as download to click a File to submit into the competition.  This needs to be in the executive summary or the very beginning of the site.  it should be obvious when you visit the site.  
  
The user tried to submit a downloaded file and got: "Predicted values must be in range [0,1]".  Provide a unique name and a short comment for the submission Note.  
  
Create an executive summary subpage that explains exactly how to make a submission into the contest.  
  
Work on the next steps from the previous sessions first.  The goal is to place top of the leaderboard.  
  
Understand the problem; download data from the data tab (requires login); create and train your own model; generate predictions that match the submission format.  
  
Before implementing, generate 3–5 candidate geological hypotheses we haven't tried yet, each naming: the specific layer(s) involved, the physical signature being targeted, why it should catch a fault missing from the USGS/INGENIOUS catalogue rather than one already in it, and how it differs from anything already implemented in this repo.  Rank them by expected DTI improvement and implementation cost.  Validate the top candidate on our spatially-blocked holdout set before touching a weekly submission slot — do not spend a submission slot on an idea that hasn't beaten the current holdout best.  If a candidate can't be validated without new external data, name the specific free, official source needed and check it's obtainable before proposing the idea as viable.  
  
Site creation: Create a github page for this repo that has clean ui, user friendly, simple and easy to use.  It should include all relevant information in an easy to read format with official verified links as sources for review.  
  
The single remaining blocker to training is data placement: run `bash scripts/download_competition_data.sh` on any unrestricted machine into `data/`, then `python scripts/prepare_data.py` — after that the full train→inference→validate pipeline is ready to run (GPU needed for training; metric/losses/validation all verified working here on CPU).  
  
Run this task through multiple passes (Pass 1 implement; Pass 2 review/fix; Pass 3 re-check).  
  
Go ahead and create a pull request and then merge the pull request onto the main. Make suggestions for what work still needs to be done and any limitations that is in the way of a successful project.  Work line by line verify everything no hallucinations.
~~~~

See `docs/` for the full site and `docs/methodology.html` for the hypothesis ranking and validation.
