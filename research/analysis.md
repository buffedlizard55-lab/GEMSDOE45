# Scientific assessment: why 0.2778, and what could improve it?

## Evidence classes

**Official fact:** the task, scoring formula, data description and rules were read from DrivenData and the complete 19-page September 2026 rules PDF on 2026-10-06. The public leaderboard showed 0.3345 at rank 1 and 0.2778 at rank 13. The public table does not identify the submitted file responsible for 0.2778. The mapping of that score to GEMSDOE32's named artifact is **user-reported**, not independently organizer-authenticated.

**Measured here:** restored file hashes, raster geometry, finite/range/mask tests, freshly computed spatial holdout results, and candidate/prior pixel comparisons.

**Hypothesis:** geological interpretation of our differential features and causal reasons for a leaderboard improvement. Hidden truth and per-submission organizer receipts are unavailable. No honest method can derive their exact performance from a list of aggregate scores.

## The useful lesson from GEMSDOE32

Its landing source and downloaded raster describe H33-2-B2: the existing 0.2708 dotted surface with predictions within 200 m of the mapped catalogue removed. The downloaded raster contains 37,654 positive pixels. The reported predecessor contained 40,199 dots. The public description does not train a new geological detector. Its 'live mirror' generates assumptions about unseen faults and is not an independent holdout of those faults.

Let T be the maximum distance-weighted coverage summed over truth pixels, S the total emitted prediction mass, M the emitted mass weighted by proximity to truth, and G the number of truth pixels. The official equation simplifies to:

```
DTI = T / [0.2 (T + S - M) + 0.8 G]
```

This follows from FP = S - M and FN = G - T. A removal that reduces S while preserving T and M increases DTI. Removing redundant dots can therefore outperform a denser raster without discovering a single additional fault. If a removal loses coverage a, removes mass b, and loses proximity mass c, it helps exactly when:

```
a (1 - 0.2 DTI) < 0.2 DTI (b - c)
```

The catalogue-flank pruning result is consistent with this mechanism: new-fault labels may be less common immediately beside previously mapped traces. It does **not** prove that no new faults occur within 200 m, nor that the entire catalogue is masked out by the server. Those claims appear in some prior sites but are not established by the published problem description. A newly mapped parallel strand could be close to an old fault.

Dotted emission can preserve much of 300 m kernel coverage while paying for fewer positive pixels. It does not turn a prediction into a physical fault trace. Indiscriminate pruning eventually loses enough T to hurt the score. Simply extending B from 2 to 3 or adjusting dot spacing is not a new discovery mechanism.

## Can we exceed 0.2778 or 0.3345?

Possibly, but neither a guarantee nor a numeric forecast is justified. A higher score requires a better coverage/false-positive trade-off on the unseen truth. The current 0.3345 leader is about 20.4% higher in DTI than 0.2778; this is **not** a claim that 20.4% more geological information suffices. We do not know either model's hidden TP, FP or FN.

The rational next experiment changes the detector rather than continually reallocating the same dots. H45-PHASE tests whether parity of differential responses distinguishes step-like structural boundaries from symmetric topographic features. At scales sigma=2 and 5 pixels, define:

```
O = sigma * sqrt(f_x^2 + f_y^2)
E = sigma^2 * sqrt(f_xx^2 + f_yy^2 + 2 f_xy^2)
F = O / (O + E + 1e-6)
Q = sigma^2 * (f_xx + f_yy) / [sqrt(2) (O + E + 1e-6)]
P = Q_sigma2 * Q_sigma5
```

Here f is a Gaussian-smoothed field; implementation applies sigma scaling to its derivatives exactly once. F measures differential odd-versus-even dominance; Q retains curvature sign, and P measures signed phase persistence. These are **differential parity descriptors**, not a Hilbert-transform phase-congruency estimator. Ratios are bounded; the classifier sees gradient and Hessian energy too, so a ratio in nearly flat terrain need not be treated as strong evidence. Hypotheses are frozen in `research/hypotheses.md` before evaluation.

We fit a small gradient-boosted classifier to the public catalogue. Unlabelled pixels may contain faults: this is noisy-label learning, **not** a mathematically validated positive-unlabelled probability estimator. The model output is a ranking, and deterministic spacing converts that ranking into a binary GeoTIFF. No old submission is blended in. The final raster targets all faults, as the official brief requests, rather than imposing a catalogue exclusion rule.

## Validation that does not manufacture a win

Four geographic quadrants are held out in rotation. A separate quadrant selects prediction budget, and the other two fit the classifier. A 3 km exclusion on both sides of split lines separates domains by at least 6 km. The same training samples, model settings, packing rule and budget choices are used in baseline and phase arms. The phase arm must beat the baseline in all four test quadrants to pass our local gate.

This is a new **catalogue-transfer proxy**, not the non-existent inherited GEMSDOE45 holdout, nor the organizers' private split. Spatial withholding prevents fitting labels in the test region, but does not make visible catalogue faults representative of newly discovered faults. Selecting a method after observing these four outcomes would also consume the holdout; future tuning needs fresh locked confirmation regions or independent fault labels.

The old H33 artifact is unsuitable as a fair labelled-catalogue baseline: it was intentionally constructed to avoid those labels. Beating it on them would be trivial and misleading. We therefore refuse to use that comparison as submission-slot approval. In particular, a local win does not establish that this file is better than 0.2778, let alone 0.3345. The slot gate stays closed without comparable incumbent holdout evidence and independent off-catalogue confirmation.

## Literature and physical limits

USGS describes derivative/filtering methods that delineate **both buried faults and contacts** at Northern Granite Springs Valley; that non-uniqueness is precisely why magnetic edges alone are insufficient. [1](https://www.usgs.gov/publications/a-geophysical-characterization-structure-and-geology-northern-granite-springs-valley)

DOE-hosted work describes complementary gravity, magnetic and MT responses and argues for geophysical integration with geological constraints to characterize geothermal systems. It does not validate H45 or a universal conductivity threshold. [2](https://www.osti.gov/servlets/purl/1724108)

USGS documents how lidar profiles at Hebgen Lake reveal both modern and prehistoric offsets, as well as scarp gaps and preservation effects. A missing topographic scarp does not imply a missing fault. [1](https://www.usgs.gov/observatories/yvo/news/using-modern-tools-look-past-earthquakes-how-lidar-data-help-better)

The official GeoDAWN release documents 400 m traverse spacing in Area 2 and variable terrain clearance. A 100 m raster is not evidence of independent measurements every 100 m: resampling cannot create physical resolution. https://doi.org/10.5066/P93LGLVQ

The experiment does not claim that a fault necessarily has hot fluid, a permeable reservoir, a vent, or commercially useful heat. Those need additional thermal, hydrological, geochemical and geological evidence. A regional 100 m detector cannot replace 1 m scarp mapping, expert interpretation or field verification.

## Irregularities and limitations

1. Empty starting repository: the supplied claim that an existing GPU train→inference pipeline was ready was inapplicable here. We built and ran a new CPU pipeline.
2. Authenticated organizer downloads still redirect to login. Owner mirrors are checksum-pinned, not organizer-authenticated. Dropbox and ordinary HTTPS requests failed in this sandbox; GitHub API restoration succeeded. No TLS validation was disabled.
3. The mirrored sample contains 60,988 ones coincident with public catalogue positives, unlike the official example description of all zeros. Its values are never used as predictions; only grid and footprint are used.
4. Mirror metadata are not authoritative scientific descriptions. In particular band 6 `tc` has an ambiguous description ('tilt angle or total curvature'); the official page also illustrates total radiometric counts. It is an uninterpreted raw input here. Band 15's mirror 'basement' wording conflicts with the competition's 'conductive base' language. No physical interpretation of those two bands is used by H45 transforms.
5. The portal error's actual server path cannot be diagnosed without a receipt or validator source. NaN, infinities, and negative sentinels are possible causes, not proven facts about your specific rejection. Our output has finite stored values throughout [0,1], no nodata sentinel, zeros outside footprint, and an internal GDAL mask marking that footprint invalid. This reconciles null-area semantics with a raw all-finite check, but portal mask handling is untested.
6. 42GEMSDOE/43GEMSDOE/44GEMSDOE resolve to GEMSDOE42/43/44. GEMSDOE43 and GEMSDOE44 only had README titles at review. The inventory preserves failed aliases as well as resolved names.
7. A source scan is not a proof that every line of every prior repository was reviewed. TIF uniqueness is explicitly scoped to the inventory. Original implementation plus changed pixel predictions is stronger than a renamed file, but cannot certify absence of an unobserved duplicate.
8. No weekly slot used. No competition authentication or personal eligibility certifications were supplied or requested. A local format receipt is not organizer acceptance.
9. External official data licensing matters. GDR1391 is CC BY 4.0; GeoDAWN release is marked CC0. Mirror access alone does not establish organizer redistribution rights. Raw data stay ignored; no claimed grant of rights beyond the sources.
10. Generative AI assisted research, software and documentation. Disclose this use in any finalist narrative as required by rules section 3.2. The classifier is a conventional numerical model, not generated imagery.

---

## Addendum (2026-10-06, H51): what the live ladder actually licenses

Everything in this addendum is first-hand measurement on rasters retrieved byte-exact from their
own repository blobs, or algebra on the published metric. Scores are owner-reported
(`user_report` in `research/score_claims.csv`) and have **no organizer receipt** (IR-45-004).

### 1. The three best members of the family are nested pixel sets

Verified by exact integer-coordinate set inclusion, not by trusting the descriptions:

```
h33-2-b2        37,654 px   live 0.2778   strict subset of h27-4
h27-4-solo-d28  40,199 px   live 0.2708   strict subset of dotted-d2-8
dotted-d2-8     44,090 px   live 0.2600
h25-1-d1-5      60,069 px   live 0.2477   superset in the same line of descent
```

So the leader of this family is the **most heavily pruned** member: 6,436 predicted pixels were
deleted from `d2.8` in two steps, and each deletion **raised** the live score
(+0.0108 for 3,891 removed, +0.0070 for 2,545 removed). Per removed dot the gain was
2.78 × 10⁻⁶ and 2.75 × 10⁻⁶ — the same to two significant figures across two independent removals
at different exclusion radii (≈1.4 px and ≈2.0 px).

### 2. Why deletion can raise the score, exactly

Write the metric as `DTI = T / (0.2T + 0.2F + 0.8K)`. Deleting a dot at x removes its contribution
`t` to `T` and reduces `F` by `1 − k(x)`, where `k(x)` is the kernel weight it was receiving from the
nearest truth pixel. Substituting and comparing:

```
delete  <=>  0.2 * (t + 1 - k)  >  DTI * t
```

For the family's files `t` is of order 0.1–0.3 and `k` of order 0.3, so the left side is ≈ 0.18
against a right side of ≈ 0.05: **deletion is profitable whenever the removed dot is not carrying
much unique coverage**. That is the entire mechanism. It is also why simply drawing more dots
cannot fix a weak detector: every extra dot adds `0.2` to the denominator and only `t` to the
numerator.

### 3. Count is not the explanation — placement is

Three controls from the same corpus rule count out:

| file | px | live | what it isolates |
|---|---|---|---|
| `h34-scatter-q50-arr-matched` | 37,654 | **0.0778** | the *same count* as the leader, displaced — only 644 px shared, mean nearest-neighbour displacement 9.74 px |
| `h35` | 39,530 | **0.0418** | nearly solid mass (mean NN 1.03 px, 99.9 % of dots within 3 px of another) — under a MAX-over-predictions TP_w, clustered dots self-cancel and are pure FP_w |
| `r13-lattice-s5` | 206,895 | **0.0904** | a content-blind uniform lattice: an unbiased but weak probe, worth about a third of the tuned family |

Same count, same marginal band statistics, same distance-to-catalogue profile — and a 20-point
DTI gap. **Which pixels survive is what is being scored.**

### 4. What this licences, and what it does not

The ladder licenses one statement: *do not emit more mass than the 37,654-dot member of this
family*. It does **not** license:

* a claim that continuing to prune at 2.8 × 10⁻⁶ per dot reaches 0.30 or 0.3345. The constant is
  measured over two adjacent removals totalling 6,436 dots; the productivity of the remaining dots
  must rise as the cheap deletions are exhausted, and where the curve turns is unmeasured;
* a rank for any candidate built on a **different** candidate set. The ladder is a within-family
  dose-response.

### 5. Why the emission emulator cannot rank candidates either

`evidence/h51_live_fit.json` fits the expected-coverage emulator to 46 unique live-scored rasters
over a grid of `a_tip ∈ {0, 0.5, 1, 2}` × `r_zero ∈ {1, 2}` × `K ∈ {2,000 … 34,000}`. Best point:
RMSE 0.1374, **Pearson +0.176**, R² +0.031. It ranks the live-best file (0.2778 → model 0.0863)
*below* a file that scores 0.0904 (model 0.4158), and it cannot separate `h33-2-b2` (0.2778) from
`h34-scatter` (0.0778) even though a multi-physics corroboration field was added specifically to
separate them. The field statistics at the two files' dots are nearly identical
(`det_elev_slope` percentile 0.855 vs 0.890). **Conclusion: the supplied 19-band stack does not
separate mapped from unmapped faults at the 300 m kernel's scale, and no blend of it that we tested
does so either.** That is a negative result about the data, and it is the honest reason H51 ships as
an unscored candidate rather than as a claimed improvement.

### 6. The one thing that would change the answer

An independent, off-catalogue inventory of *known-real* faults — USGS QFAULT
(<https://earthquake.usgs.gov/ws/qfault/>) or the USGS State Geologic Map Compilation
(<https://mrdata.usgs.gov/geology/state/>) — would supply the only available population resembling
the hidden truth. Both are free and official; both were checked for obtainability from this
environment and **failed** (HTTP 000 on plain HTTPS), alongside drivendata.org, gdr.openei.org and
dropbox.com. That remains the highest-value unblocked action, and it is a network problem, not an
analysis problem.
