# Three-pass implementation review

## Pass 1 — implement and measure

- Confirmed the starting checkout had only an 11-byte README: no inherited pipeline or holdout existed.
- Read official problem/metric/submission documentation, resources, leaderboard and every extracted chunk of the September 2026 rules PDF.
- Inspected prior landing sources, restored checksum-pinned feature/label/template mirrors and the GEMSDOE32 reference artifact. Resolved reversed aliases for repositories 42–44.
- Registered four geological hypotheses and a frozen four-quadrant experiment before implementation. Created an original differential phase-parity detector, not an edit of a prior submission.
- Ran the complete CPU prepare→fit→calibrate→holdout→full inference→write→re-read pipeline.
- Measured baseline mean DTI 0.1930835521 and candidate 0.2085207088; positive paired differences in all four quadrants. Generated 82,678 positive cells.
- Created a unique finite-value GeoTIFF, internal footprint mask, ZIP, note, sources, score register, reproducible pipeline and download-first website.

## Pass 2 — assumptions, bugs and edge cases

- Independent brute-force implementation agrees with the optimized DTI on random binary and probabilistic grids. Tests include empty truth/predictions, perfect predictions, exact kernel boundary, masking, NaN/infinity/out-of-range rejection, spacing and deterministic ties.
- Verified phase bounds, polarity reversal, and a synthetic step-versus-ridge sanity check. These establish mathematical behavior, not geological efficacy.
- Fixed the inventory-refresh list to retain resolved GEMSDOE42/43/44 on future refreshes, rather than losing those sources after the initial alias failure.
- Rejected NaN/sentinel output as the default; all stored cells are finite, outside footprint is zero, and the internal mask conveys invalid footprint cells. Portal acceptance remains untested.
- Refused to use a catalogue-pruned historic artifact as a misleading visible-label holdout baseline. The local proxy gate and competition-slot gate are separate.
- Fingerprinted 311 unique prior Git blobs in the inventory: 309 compatible single-band grids, 2 incompatible grids/bands, zero unavailable downloads. None of the compatible normalized prediction arrays is identical to H45. Outside-only NaN/zero changes are normalized and cannot fake uniqueness.
- Flagged ambiguous mirror metadata (`tc`, depth-to-base), positive example-template values, source authentication limits, stale leaderboard target, and selection of the highest budget in every calibration fold.
- Source refresh fails closed on login/schema/network errors. The sandbox HTTPS failures are logged, while dated research-tool verification remains separate. No TLS checks were disabled and no authentication secrets were requested.

## Pass 3 — complete request cross-check

- Checked released TIFF bytes against the stored SHA-256, grid and mask; checked ZIP contains exactly the same TIFF bytes.
- Recomputed all reported holdout DTI values from TP/FP/FN, checked calibration choices, paired deltas and protocol hash.
- Added automated checks for all generated local site links, accessibility basics, stale evidence copies and future source-feed parsing failures.
- The README preserves a normalized standing project prompt, all supplied score observations, official links and supplied data links. It explicitly says this is not a verbatim transcript. AGENTS.md directs future sessions to read it first.
- GitHub Actions tests the release and deploys Pages; a daily job refreshes source health without retraining or consuming submission slots. A separate reproducibility workflow can regenerate the research TIF on CPU.
- Repeated the entire frozen fit/calibration/holdout/inference stage independently: identical fold results, identical prediction hash and byte-identical TIFF. See `evidence/reproducibility.json`. A repeat on the same data is a reproducibility check, not fresh statistical confirmation.
- Local HTTP preview and TIFF download returned 200, and downloaded bytes matched the release hash. Browser visual automation could not run because the Chromium download was blocked by sandbox network access; no screenshot review is claimed.
- Existing GitHub Pages is configured for legacy main/root publication. The integration denied a Pages configuration update (403); root redirect pages and `.nojekyll` provide a compatible fallback without changing that setting. Workflow deployment is independently checked after merge, not assumed from configuration.
- Final PR/deployment details are recorded in `evidence/delivery.json` when available. Do not infer success merely from workflow configuration.

## Requirements deliberately not claimed as satisfied

- No proof of a competition score above 0.2778, 0.3195 or the now-observed 0.3345.
- No comparison to a nonexistent inherited GEMSDOE45 holdout best; no independent new-fault test labels; no actual competition upload or use of a weekly slot.
- No organizer authentication of mirrored data, no guarantee of portal acceptance, and no blanket line-by-line verification of every historical repository or every geological interpretation.
- No exhaustive uniqueness claim beyond the publicly inventoried raster blobs; no new 1 m DEM or external geothermal field study was performed.
- The README prompt is normalized, not a complete verbatim archival transcript. Repeated instructions and malformed duplicate links were consolidated transparently.

These limitations are visible on the website and in the executive summary rather than hidden behind a passing test count.
