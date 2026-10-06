# H51 candidate hypotheses — registered before implementation, with measured verdicts

Registered 2026-10-06 on branch `arena/3d712b92-gemsdoe45`, before any H51 code was written.
The purpose of this file is to fix the hypotheses, the ranking rule and the promotion gate *before*
seeing results, so that the results below cannot be retro-fitted. Verdicts were added afterwards
and are marked as such; nothing above a verdict line was edited.

## Official facts the hypotheses are built on (quote-verified, links in `docs/sources.html`)

| # | Fact | Source (verified by reading the page on 2026-10-06) |
|---|---|---|
| F1 | "Pixels corresponding to known USGS/INGENIOUS faults are masked / excluded from evaluation, so they do not count towards penalty terms." | chrisk-dd (DrivenData staff), 2026-09-16, [thread 11516](https://community.drivendata.org/t/scoring-clarification-are-known-usgs-ingenious-faults-masked-when-scoring-and-are-they-in-the-final-round-label-set/11516) |
| F2 | "The mask is indeed pixel-exact - it is identical to the provided set of training fault labels." | chrisk-dd, 2026-09-21, thread 11516 post 4 |
| F3 | "A predicted pixel that is near a known fault trace but far from a new-fault ground truth pixel will be fully penalized, i.e., the buffer does not apply to known faults." | chrisk-dd, 2026-09-21, thread 11516 post 4 |
| F4 | "A new-fault ground truth pixel can indeed lie within 300m of a known fault trace. Such pixels would constitute corrections or modifications to existing fault traces. Identifying these corrections is one outcome we are aiming for as part of this competition." | chrisk-dd, 2026-09-21, thread 11516 post 4 |
| F5 | "For the purposes of this competition, 'new fault' means 'any fault pixel not already captured by USGS/INGENIOUS' and can include newly mapped geometry of an existing fault system." | chrisk-dd, 2026-09-23, [thread 11536](https://community.drivendata.org/t/where-do-you-draw-the-line/11536) |
| F6 | Metric, kernel and submission format: `TI = TP/(TP + alpha*FP + beta*FN)`, `alpha = 0.2`, `beta = 0.8`, triangular kernel `k(d) = (1 - d/300m)+`. | [Problem description](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/) |

Derived here (algebra, not assumption): because `TP_w + FN_w = K` identically,
`DTI = T / (0.2T + 0.2F + 0.8K)`, and adding mass at a pixel of realised kernel weight `w`
improves the score iff `w > 0.2 * DTI`. At DTI = 0.2778 that bar is 0.0556, i.e. a dot is worth
emitting iff it is within 283 m of a hidden-truth pixel. `tests/test_metric.py` checks the identity
numerically. Consequence: **binary emission of a selected subset is optimal** — and all 46
live-scored owner rasters profiled in `research/live_score_features.csv` are exactly binary.

## Ranked hypotheses

Ranking rule, fixed before results: expected DTI gain per unit of implementation cost, with a hard
requirement that the top candidate be testable with the data already restored here.

### H51-A — Catalogue fault-zone strand completion (perpendicular offsets) — **RANK 1**
* **Layers.** `det_elev_slope` (band 19) ridge response at sigma = 1.2 and 2.5 px; `iso_grav_anom_hg`
  (band 18) and `tmi_hg` (band 3) edge response; plus catalogue geometry (skeleton + local strike).
* **Physical signature.** A 1-px line response sampled on lines drawn *parallel* to the mapped
  trace at perpendicular offsets of +/-2 and +/-3 px (200-300 m), i.e. the hanging-wall/footwall
  strands of the same fault zone.
* **Why it should catch a fault the catalogue lacks.** F5 says a newly mapped parallel strand
  inside an existing fault zone *is* a new fault. Reconnaissance mapping (the USGS/INGENIOUS
  compilation) draws the master trace; splays and synthetic/antithetic strands of the same zone are
  the classic omission. No pixel of such a strand is in the catalogue, so F1 does not mask it, and
  F3 says it will be penalised if no truth pixel is nearby — so precision matters and the ridge
  gate is essential.
* **How it differs from everything already implemented here or in the sibling repositories.** The
  GEMSDOE family has shipped along-strike *tip* families (H32-1 prethin-tip-euler, H33-D
  analog-tip-stepover) and dotted unions of a ridge surface (H19-4/H19-5/H25-1/H27-4/H33-2-B2).
  None enumerates *perpendicular* strand positions inside the fault zone. The nearest neighbour is
  GEMSDOE33's "stepover" arm, which still emits along strike from tips.
* **Expected gain / cost.** Highest of the five; medium CPU cost.

### H51-B — Kaplan-Meier along-strike continuation, per fault — **RANK 2**
* **Layers.** Catalogue geometry alone (tip position, terminal strike, segment length); KM survival
  curve of the catalogue's own relay-gap sample; `det_elev_slope` as a veto.
* **Physical signature.** The product-limit estimator `S(e)` of the catalogue's own along-strike
  tip-to-nearest-other-segment distances is, verbatim, the probability that a structure is still a
  fault `e` px beyond a mapped tip. It is fitted here from 3,866 tips: 3,163 events, 703
  right-censored, median 9 px, restricted mean 37.0 px, S(1) = 0.896, S(3) = 0.667, S(9) = 0.489,
  S(20) = 0.392.
* **Why it should catch a fault the catalogue lacks.** F5 says an unmapped continuation past a
  mapped endpoint is a new fault. Tips are where mapping stops, not where faults stop.
* **How it differs.** Every prior tip arm used a *fixed* radius or a fixed packing (B = 2 px,
  d = 2.8 px, r = 30 px). Here each tip draws its own length from the fitted curve.
* **Expected gain / cost.** Medium-high; low cost. Shares its emission engine with H51-A.

### H51-C — Correction annulus (misaligned-trace recovery) — **RANK 3**
* **Layers.** The same derivative fields, sampled in the 1-3 px annulus *around* the mapped trace.
* **Physical signature.** The problem description states that "portions of the existing fault data
  may be misaligned from the true location of the surface fault"; F4 says corrections count as new
  truth. A misaligned trace shows a gradient edge offset by 1-3 px from the drawn line, on the same
  strike.
* **Why it should catch a fault the catalogue lacks.** It is by construction the population F4
  describes.
* **How it differs.** H27-4/H33-2-B2 *deleted* the 1-2 px annulus on live evidence that it was
  unproductive; nobody has attempted to *re-place* mass there at a better estimate of the true
  position.
* **Expected gain / cost.** Uncertain-to-medium; medium cost.
* **Registered risk (flagged before the result).** This hypothesis is in direct tension with the
  live dose-response that motivated the 2 px novelty gate in H51-A/B. Having both in one artifact is
  inconsistent, and the inconsistency is recorded rather than smoothed over.

### H51-D — Blind-lattice truth-density calibration — **RANK 4, NOT NOVEL**
* **Idea.** A content-blind uniform lattice earns the same expected credit wherever the truth is,
  so one live lattice score inverts to the hidden truth density.
* **Why it is ranked last.** It is a *measurement* arm, not an emission arm, and it has already
  been done and published in this family: GEMSDOE24 states the spacing-5 lattice scored 0.0904 and
  inverts to "0.24 % of cells" (the r13-lattice file is in `research/scored_blob_manifest.json` at
  206,895 px and live 0.0904). Re-deriving it would be repetition, so it is used here only as an
  independent order-of-magnitude cross-check on the modelled `K`.
* **Cost.** Low.

### H51-E — Independent off-catalogue fault inventory — **RANK 5, BLOCKED**
* **Idea.** Score against, and emit from, a fault inventory that is genuinely independent of the
  supplied catalogue: USGS Quaternary Fault and Fold Database
  (<https://earthquake.usgs.gov/ws/qfault/>) and the USGS State Geologic Map Compilation
  (<https://mrdata.usgs.gov/geology/state/>).
* **Why it is not proposed as viable here.** The candidate cannot be validated without that data,
  and the named free official sources were checked for obtainability **from this environment**:
  `curl` to `mrdata.usgs.gov` and `earthquake.usgs.gov` both returned HTTP 000 (TLS handshake
  failure) on 2026-10-06, the same failure this repository already recorded for
  `www.drivendata.org`, `gdr.openei.org` and `www.dropbox.com`. Only GitHub and PyPI were
  reachable. Obtainability is therefore **not** verified and the hypothesis is not viable from
  here. It is also not novel: GEMSDOE29 already shipped an unscored `sgmc-off-catalogue-44k`.
* **Cost.** High; needs an unrestricted network.

## Promotion gate (fixed before results)

A candidate may not consume a weekly submission slot unless it beats the current holdout best on a
spatially blocked instrument that is *not* the catalogue in disguise. The gate is CLOSED for this
session, for the reason recorded in `research/analysis.md` §"Instruments": the catalogue-truth
holdout inverts the live order, and the multi-physics field fitted here was measured to have no
rank correlation with live scores. The artifact is therefore published as an **unscored candidate
under HOLD**, with its format audit, its uniqueness audit and its limitations, and no score claim.

## Verdicts (added after the runs; nothing above this line was edited)

| ID | Verdict | Measured |
|---|---|---|
| A | **Implemented and shipped** (part of the H51 artifact) | 236,798 strand candidates at offsets -3,-2,+2,+3 px from a 59,209-px skeleton before gating; see `evidence/h51_submission.json` |
| B | **Implemented and shipped** (same artifact) | per-tip draws from the fitted curve: see `extension_lengths` in `evidence/h51_submission.json` |
| C | **Partially implemented, and in tension with the novelty gate** | the +-2/+-3 px offsets are the annulus; the live-derived 2 px novelty gate removes the 1-2 px part. Recorded as IR-45-010 |
| D | **Used as a cross-check only** | modelled K = 12,000 px sits between the two independent estimates 8,936 (lattice inversion, this repo's algebra with E[k] = 0.5) and 12,632 (GEMSDOE24's published 0.24 % of cells) |
| E | **Not viable from this environment** | HTTP 000 measured for both named sources on 2026-10-06 |
