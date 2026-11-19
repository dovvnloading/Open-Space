# Independent-referee-style review — ACT DR6-era depth-1 search

**Verdict:** **Accept only as a public-data inventory and preliminary software test, not as an astrophysical transient search, detection, null result or rate-limit paper.** The claim-first README is appropriately restricted. The original science objective remains unfulfilled.

**Scope reviewed:** frozen catalog manifest and SODACAN receipt; `freeze_mapcat.py`, `freeze_scoped.py`, `download_products.py`, `mf_core.py`, `run_pair.py`, `test_mf_core.py`; GitHub Actions evidence including success of eight synthetic tests. No full FITS observation data, independent map-minus-coadd products or real astrophysical candidates were submitted for review. Nothing reported here qualifies for `[REMARKABLE]`.

## Falsification and limitations

| Test / acceptance criterion | Independent assessment | Required before scientific escalation |
| --- | --- | --- |
| Source list is real and frozen | **Supported at catalog level.** Public MapCat SQLite 294,342,656 B with SHA-256 `c4e73e…8efd` was retrieved by Actions. 31,261 deterministic JSONL rows; selection of 27,512 rows at 90/150/220 GHz. | Re-run on pinned source; make compressed manifest and scoped SHA reviewable. Verify that each entry actually belongs to the intended DR6 release and no relevant products are omitted. |
| Full survey denominator | **Not measured.** Rows are map/array records, not independent spacetime exposure. | Validate FITS WCS, usable ivar pixel counts, masked area, duplicated frequencies, sky coverage, effective observation duration. |
| Correct map reference | **Not demonstrated.** The core rejects an unasserted reference and the pilot demands metadata, but the *data* of a leave-one-out coadd were not built or tested. | Provenance and checksums for an observation-excluded coadd. A contaminated reference can attenuate an astrophysical transient and bias error bars. |
| Matched-filter calibration | **Code only.** Diagonal inverse variance, identical pixel grids and a supplied PSF give the correct formal white-noise weighted least-squares amplitude on synthetic tests; stripe-correlated ACT noise and transfer functions can radically alter significance. | Measured beam (not the synthetic Gaussian), map units/normalization, correlated residual-noise spectrum, scan-direction response, multi-band response, finite-patch/edge calibration. |
| 8 synthetic unit tests | **Passed** in [Actions run 44442148269](https://github.com/dovvnloading/Open-Space/actions/runs/44442148269). | Still need real-data injection/recovery and null-distribution characterization, including false positives in split/jackknife and sign-flipped data. |
| Published transients | **None recovered from this data.** Naess 2021, Li 2023 and Biermann 2024 (published 2025) were cited only. | Freeze event coordinates/times and original catalog IDs, prove their sky pixels/epochs eligible, compare measured positions, light curves, fluxes and beam significance to source papers. |
| Moving/known-source vetoes | **Not run.** | JPL Horizons topocentric positions; ACT, SPT, Planck radio/variable source catalogs; Gaia, VLASS, RACS; Swift/GCN, TNS and contemporaneous literature for every surviving candidate. |
| Detection significance | **Not calibrated.** White-noise pixel SNR is neither a trials-corrected search significance nor an astrophysical event probability. | Exposure-matched empirical null, independence-aware peak clustering, injection completeness by flux/duration/position and correction for selection trials. |
| Rate limit | **Mathematically gated, but absent.** The software refuses to calculate a numerical zero-event limit unless passed a validated effective exposure. | Real null with (N=0), `E(F, duration)` from exposure/completeness, and a stated event model; clustered observations and varying source durations need separate treatment. |
| Novelty | **Not established.** Biermann already reported 34 ACT single-observation transient events at 27 sky positions; the 2017–2022 epoch overlap is direct. | Compare any measured new candidate at the same epoch/position with all prior ACT transient tables and cataloged variable-source outbursts, document mismatches. |

## Adversarial code review

1. A manifest file SHA identifies catalog rows, not the map/ivar bytes: **download staging and per-product hashes remain mandatory**. A public SQLite table may contain absent, relocated or different-version products.
2. `mf_core.matched_filter` assumes independent observation and reference noise; a coadd containing the tested observation violates the stated formula, even if the sky signal cancels visually. The wrapper's JSON receipt is a provenance declaration and is not an independent proof; require hashes of both reference FITS images and trace back their constituent scans.
3. `fftconvolve` is a correct *simple white-noise* reference implementation, not the usual scan-correlated ACT matched filter. Null stripes, atmospheric residuals, intrinsic source spectral differences, and local non-Gaussian noise can mimic or hide peaks. A non-detection of this uncalibrated filter cannot constrain astrophysical rates.
4. `test_synthetic_injection_recovery` injects nine bright sources into artificial white-noise maps. Nine successes are **not** a completeness curve; injected sky positions, fluxes and durations must sample the actual multi-epoch data and be run blind through all vetoes.
5. The 90, 150 and 220 GHz array catalogs contain partially synchronized scans. Do not count their signal/noise detections as independent evidence; cross-band agreements need covariance and accurate bandpasses.
6. Naive-string catalog times are not validated against original FITS time maps; UTC versus TAI/TT, per-pixel mean time and scan duration matter for Horizons or fast afterglows. Median/start/stop map times cannot replace each source pixel's observed time.
7. All candidate-source names and taxonomy (orphan GRB, TDE, magnetar, blazar) require multiwavelength and association calculations, not inference from mm brightness alone.
8. The user-requested scientific outcome (novel mm transient **or** calibrated null with quantitative rate limit) **cannot be awarded** for this scope. It would be misleading to label this an astrophysical null or log it in `docs/REMARKABLE_FINDINGS.md`.

## Review decision / disposition

**For a methodology/inventory-only contribution:** technically reproducible at the catalog level with checkable source and manifest digests, and explicitly labeled as incomplete; acceptable to merge as **inconclusive preparation**, if the scoped file and receipt are committed, the synthetic tests pass, and the root catalog wording does not imply a scientific null.

**For the user's discovery/rate-limit goal:** major revisions / more data required. No real FITS analysis, journal-quality negative control, published-event recovery or completeness denominator has been completed. The present work is neither a survey result nor an evidence-based upper limit.

**Reviewed 2026-11-18 (America/New_York)**. This is an adversarial methodological referee-style assessment prepared within the same task; it is **not** independent external peer review.
