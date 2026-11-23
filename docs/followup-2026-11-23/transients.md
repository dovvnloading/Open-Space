# Nuclear-transient follow-up: 23 November 2026

Base: `972cd5169eea2b83d7407df84eb60c167b3c084e`. This review advances the NEOWISE and JWST leads with numerical replays, an uncertainty check, a fresh archive query, and an additional public JWST epoch. **The seven infrared flares remain candidate dust-obscured transients; COSMOS-022 and GOODS-S 185 remain unclassified nuclear brightenings. None is a newly confirmed black hole or a classified tidal disruption event (TDE).**

## What changed

1. A public **4 February 2026** GOODS-S observation was absent from the existing historical search, which enumerated programs 1180, 3215 and 6541. Program **5324, observation 27** supplies valid small-aperture measurements in four filters. The nucleus was near its earlier state in February and brighter in October. This provisionally improves the timing constraint, with the imaging limitations below.
2. All **113 NEOWISE table rows** reproduce from deposited intermediate measurements, and the stated cuts select the same seven objects. All **seven selected JWST analysis scripts** reproduce their frozen outputs within numerical tolerance; three focused referee test methods also pass. These are computational checks on deposited evidence, not independent observations.
3. UGC 2653's possible radio rise is sensitive to relative calibration: the 2024-minus-preflare contrast is **0.48354 mJy/beam**, or 3.22 times the statistical error. A scenario with independent 10% epoch gain uncertainties reduces this to **2.40 times the propagated error**.
4. The NEOWISE result page now corrects an interpretation error: **coronal lines can occur in interacting supernovae**. Their presence alone does not identify a TDE. SN 2005ip is a direct counterexample with both coronal lines and infrared dust emission ([Smith et al. 2009, ApJ 695, 1334](https://doi.org/10.1088/0004-637X/695/2/1334); [author manuscript](https://arxiv.org/abs/0809.5079)).

## Additional GOODS-S 185 epoch

### Question, observing unit and frozen method

The question was whether program 5324 supplies a valid intermediate measurement of **site 185, J2000 53.188672°, −27.812554°**, not whether a new source can be discovered. A fresh public [MAST API](https://mast.stsci.edu/api/v0/invoke) footprint query returned 113 NIRCam rows, including this epoch. The four filters are one observing epoch. The [protocol](transients-goods5324-protocol.json) fixed coordinates, products, primary F277W measurement, corroborating F356W measurement, aperture and coverage checks before inspecting the new fluxes.

The extraction uses the existing `goods-s-oct2026/measure.py` prescription: a 0.15″ aperture, 0.27–0.42″ median-background annulus, and bilinear sampling on a 0.03″ sky grid. Quoted uncertainties are 1.4826 times the median absolute deviation of 28 displaced apertures. They exclude source Poisson uncertainty, differential PSF effects, relative calibration and registration uncertainty. Fluxes are **small-aperture measurements without a total-flux correction**.

| Filter | February 4, 2026 aperture flux, nJy | October 9–10, 2026 deposited aperture flux, nJy | October minus February, nJy |
| --- | ---: | ---: | ---: |
| F150W | 56.10 ± 5.93 | 81.29 ± 1.93 | +25.18 ± 6.24 |
| F200W | 52.76 ± 3.64 | 86.38 ± 2.37 | +33.62 ± 4.34 |
| **F277W** | **44.46 ± 3.37** | **70.07 ± 1.81** | **+25.60 ± 3.83** |
| F356W | 54.22 ± 2.13 | 67.94 ± 0.77 | +13.73 ± 2.26 |

Difference errors in this table are quadrature combinations of the reported local scatters. They are **not full uncertainty budgets or calibrated discovery significances**. The previously deposited F277W values were 47.88 ± 0.92 nJy on 6 October 2025 and 51.59 ± 1.65 nJy on 18 November 2025. The February measurement is consistent with the earlier low state at this precision. Common coordinate offsets of −0.02, 0 and +0.02″ in each axis give February F277W fluxes of **40.94–46.30 nJy**, still below the October high state. This is a coordinate sensitivity test, not an astrometric calibration.

The new observation therefore supports a **provisional last-low to first-high interval of 4 February–9 October 2026 (247.492 observer-frame days in F277W)**. It does not date the start, peak, or duration of the event, and intermittent variability inside the gap remains possible. The host still has a photometric redshift, not a transient spectrum. Source-specific PSF effects, calibrated relative astrometry and native-exposure checks for program 5324 remain to be completed before promoting this timing constraint.

**Local comparison-source check.** After inspecting the target, a separate diagnostic used the existing `motion_dithers.py` rule: compact October catalog sources with catalog S/N > 28 and `nn_dist > 5`, then the eight nearest within 60″, selected before examining their February fluxes. Seven have valid positive aperture measurements in both epochs; one lacks February coverage. Their February/October flux ratios have median **0.932657** and scaled MAD **0.042629**. The nucleus has ratio **0.634565**, or **0.680385 after division by the control median**. An illustrative common multiplicative adjustment would put February F277W at 47.67 nJy. Thus a common 7% scale difference does not explain the nucleus's full change, while the control offsets demonstrate why the raw local-scatter error alone is incomplete. This small, post-inspection comparison lacks variability vetting and a full source-specific PSF/registration solution. It supports the provisional timing interpretation without supplying a calibrated detection probability. [Selected controls](transients-goods5324-controls-selection.json) · [Measurements](transients-goods5324-controls.csv) · [Receipts](transients-goods5324-controls-receipts.json) · [Code](transients-goods5324-controls.py)

The frozen [photometry](transients-goods5324-photometry.csv), [range receipts](transients-goods5324-receipts.json) and [extraction code](transients-goods5324.py) identify the four public STScI products, each requested byte range, SHA-256, observation date and pipeline/CRDS versions. Only 11,773,196 bytes were fetched for the target extraction. A receipt-checked replay from those cached bytes reproduces the photometry byte-for-byte. This replay shares the original pixels.

The February and October F277W image headers both give `CAL_VER=3.0.0` and `CRDS_CTX=jwst_1584.pmap`. Matching those metadata removes a simple processing-version mismatch as an explanation; it does not make their PSFs, source sampling or relative zeropoints identical.

### Archive coverage and priority limits

Both fresh MAST requests completed on **2026-11-23 at about 23:21 UTC** with HTTP 200 and a single complete page: 79 NIRCam footprint rows for COSMOS-022 and 113 for GOODS-S 185. No listed public product has a later observation start than the already inspected May 2026 COSMOS epoch or October 2026 GOODS-S epoch. A catalog footprint result alone does not establish usable pixels; the February candidate passed a separate finite-pixel check.

The original responses and queries are preserved: [COSMOS response](transients-cosmos-mast.json), [COSMOS receipt](transients-cosmos-mast-receipt.json), [GOODS-S response](transients-goods-mast.json), [GOODS-S receipt](transients-goods-mast-receipt.json). Their response SHA-256 values are respectively `ed1fabbca968a8fa91ae156ade2300cbe1bce4aa17a734e5730d00ab3adca4b8` and `3718ef475af1f0ac0b0f48245ce37120cf9a013e2e0366213fb035f98c5c10f7`.

Web searches on 23 November used the exact names CGCG 539-025 and UGC 2653, the two JWST positions, host IDs, and candidate names, with a second search engine and a 30-day recency check for recent reports. No inspected result established a new classification or a matching transient report. These searches are not an official TNS cone search, a complete ADS/catalog audit, or a fresh check of all seven nearby hosts. The prior [TNS guest requests](../../results/jwst-repeat-imaging-transients/REFEREE_REVIEW.md) returned throttling or unvalidated HTML; this review does not turn those failures into negative matches. Discovery priority remains unestablished. The NEOWISE label “unreported” still refers to the **14–15 November** service snapshots.

## Target list and next observations

The ordering below is an observing recommendation based on the deposited evidence and likely information gain, not a probability ranking. Coordinates for the nearby sample are those actually used in [flares.csv](../../results/neowise-nuclear-flares/data/flares.csv); use a current image for acquisition because a galaxy catalog coordinate is not a precision transient position. In particular, CGCG 539-025's prose position and table coordinate differ by about an arcsecond.

| Priority within nearby sample | Target | RA, Dec (J2000 degrees) | Evidence that motivates follow-up | Next discriminating measurement |
| --- | --- | --- | --- | --- |
| 1 | **UGC 2653 / MCG+06-08-012** | 49.65928, +37.60806 | 2022 onset; SPHEREx red excess in the 2025–26 combined product; possible 2024 radio increase | Nuclear spectrum plus calibrated, spatially resolved radio spectrum; compare later time-resolved infrared observations |
| 2 | **CGCG 539-025 / Z 539-25** | 36.72150, +41.68951 | Largest integrated infrared event; 7.1-fold W2 peak in 2018; 2024 residual; SPHEREx excess is marginal | Nuclear and adjacent-host spectrum; current infrared state; classify through line profiles and multiwavelength evolution |
| 3 | ESO 161-23 | 101.01440, −55.90860 | 2022 event, peak two-band blackbody estimate 1.3 × 10⁴³ erg/s | Southern optical spectrum and current infrared measurement |
| 4 | LEDA 204150 | 291.71210, −11.78802 | About 49% of peak W2 excess remained at the final NEOWISE visit | Current infrared state, then nuclear spectrum if an excess persists |
| 5 | LEDA 140951 | 319.37120, −80.80237 | Long-lived 2017 flare, precise WISE centroid consistency | Southern nuclear spectrum and archived host comparison |
| 6 | 2MASX J01110461−4558435 | 17.76929, −45.97878 | 2019 flare, near its old level by the final NEOWISE visit | Host spectrum; archive search for observations during activity |
| 7 | LEDA 707155 | 177.92800, −31.26870 | Broad 2017–18 flare, no remaining NEOWISE excess | Host spectrum and historical classification; avoid using its inconsistent AllWISE/SPHEREx baseline as a fading test |

“Remaining excess” in this table refers to observations through **2024**, not to the flux on 23 November 2026. SPHEREx DR1 provides a combined 2025–26 product without a unique exposure time for each extracted color. The [IRSA mission page](https://irsa.ipac.caltech.edu/Missions/spherex.html), checked on this review date, lists DR1 in November 2026 and provides individual-image spectrophotometry tools. A temporal decay claim needs the contributing exposure times and consistent calibration.

For **UGC 2653 and CGCG 539-025**, obtain spectra containing the host continuum, Balmer lines, He II, [O III], [N II]/[S II], and available coronal diagnostics such as [Fe VII], [Fe X] and [Fe XIV]. Preserve instrumental resolution, slit/IFU extraction, adjacent-host light and flux-calibration information. Interpret line ratios, widths, evolution, radio spectrum and infrared history together: a single coronal-line detection is compatible with more than one physical mechanism. Neither a late nondetection nor a steady narrow-line host alone excludes a past TDE. Exposure times require current source brightness and an instrument sensitivity calculation; none is asserted here.

For the **UGC 2653 radio test**, the first goal is to determine whether the modest change survives a matched beam, fixed-position source fit, and local comparison-source flux calibration. A subsequent pair of radio bands can constrain the spectral shape of any compact variable component. Explicitly separate the nucleus from SN 2007D, 6.3″ away, and from extended star formation. The NVSS 1.4 GHz flux and VLASS 3 GHz peak flux have different epochs and resolutions; their ratio is not a simultaneous nuclear spectral index.

The two distant JWST leads require separate programs:

| Target | Coordinates (J2000 degrees) | Measurement to secure | Decision enabled |
| --- | --- | --- | --- |
| **COSMOS-022 / COSMOS2020 953994** | 150.122033, +2.323386 | Later F162M and F182M imaging with independent dithers, plus a spectrum while elevated | Fading versus recurring/persistent nuclear variation; classification at host spectroscopic z = 0.5067 |
| **GOODS-S 185 / JADES DR2 103607** | 53.188672, −27.812554 | Validate the new February epoch on native exposures; later F277W/F356W imaging; nucleus-specific spectrum | Refine rise history and test decline; obtain the correct host/transient redshift and classification |

Each future image comparison should keep the observing epoch, filter, detector/exposure IDs, coordinates, aperture/PSF prescription and masks fixed before inspecting candidate flux. Include nearby stationary sources, signed artificial sources on the real background, and native detector/DQ checks. Report sensitivity when a predicted signal is not recovered. A second implementation on the same exposures tests robustness; a separate observing epoch supplies new physical evidence. The adjacent z = 5.2509 GOODS-S galaxy is not site 185's host.

## UGC 2653: what the radio and SPHEREx statistics support

The radio calculation reads the exact deposited [radio.json](../../results/neowise-nuclear-flares/data/radio.json). The inverse-variance weighted 2019/2021 baseline is **0.754929 mJy/beam**. The 5 July 2024 peak is **1.238471 mJy/beam**. The baseline weights are fixed to the statistical errors across all scenarios.

| Assumed error model | Error on 2024 minus baseline, mJy/beam | Contrast divided by error |
| --- | ---: | ---: |
| Deposited statistical errors only | 0.149983 | 3.224 |
| Additional 10% gain error, independent between epochs | 0.201860 | 2.395 |
| Additional 10% gain error, common pairwise correlation 0.5 | 0.181081 | 2.670 |
| Additional 10% gain error, fully correlated between epochs | 0.157585 | 3.068 |

These gain-error models are explicit sensitivity assumptions, not empirically fitted VLASS calibration errors. With early weights `w_i`, early fluxes `F_i`, late flux `F_l`, contrast `d` and fractional gain uncertainty `f`, the added variance is `f²[(1−rho)(F_l² + sum(w_i² F_i²)) + rho*d²]`. The [NRAO quick-look guide](https://science.nrao.edu/vlass/data-access/vlass-epoch-1-quick-look-users-guide) documents campaign-dependent flux biases. All three UGC products here are campaigns 1.2, 2.2 and 3.2, so this review does not incorrectly apply the special 1.1 offset. The original extractor takes the maximum pixel within 2.5″ in each image. Neither that selection bias nor source morphology is calibrated by the calculation. “Possible brightening” is the appropriate wording.

SPHEREx's 0.477267-mag color residual divided by the 0.052004-mag comparison scatter equals **9.1775**. There are **0 of 65** calibration galaxies as red. This establishes a substantial color outlier under that calibration; it does not measure a Gaussian nine-sigma false-alarm probability or identify the emitter. Those same controls help fit the calibration and do not constitute an independent extreme-tail holdout. CGCG 539-025 is only 0.118931 mag above its color prediction, with 3 of 65 controls as red. The [replay code](transients-radio-sensitivity.py) and [output](transients-radio-sensitivity.json) preserve the exact source hashes and assumptions.

## Reproduction and completed checks

| Check completed in this review | Result | Scope |
| --- | --- | --- |
| NEOWISE original `summarize.py`, all inputs deposited | 113 rows match; class counts 55 active nuclei, 27 unreported, 13 published, 13 reported transients, 5 off-nuclear; strict subset 7 | No new 550 MB IRSA download or remeasurement of single exposures |
| Source hash validation | 19 of 19 available files agree with the deposited source manifest | Includes the radio, SPHEREx, 2MRS and main flare inputs used here |
| Focused `test_jwst_pr129_referee` methods | Three methods pass: catalog/COSMOS flux, native-exposure/spectrum tables, GOODS-S candidate/blind-depth tables | Assertions on stored measurements |
| Seven original JWST analysis scripts | All regenerated outputs agree within relative/absolute tolerance 3 × 10⁻⁵ | Scratch copies; no changes to frozen result tables |
| Fresh MAST inventory | Two complete successful responses; additional February product identified | Metadata query followed by finite-pixel target extraction |
| February extraction and cache replay | Four valid measurements; exact repeat of the CSV with source-byte hashes checked | Stage-3 mosaic checkpoint; native-exposure follow-up pending |
| Local F277W controls | Seven valid pairs from eight frozen comparison positions; median February/October ratio 0.9327 versus nucleus 0.6346 | Post-inspection robustness check; source-specific systematics remain |

The initial MUSE and blind-summary replay attempts lacked three binary inputs in the sparse checkout. They were rerun after retrieving the missing deposited files and passed; these were missing-input errors, not discrepant scientific results. The [NEOWISE replay receipt](transients-neowise-replay.json) and [JWST replay receipt](transients-jwst-replay.json) retain commands, comparisons and the earlier failure reasons.

From the repository root, with NumPy, pandas, SciPy and Astropy installed:

```sh
python docs/followup-2026-11-23/transients-radio-sensitivity.py
PYTHONPATH=tests python -m unittest \
  test_jwst_pr129_referee.JWSTPR129Referee.test_catalog_census_and_cosmos_flux \
  test_jwst_pr129_referee.JWSTPR129Referee.test_native_exposures_and_spectrum \
  test_jwst_pr129_referee.JWSTPR129Referee.test_goods_candidate_and_blind_depth -v
python docs/followup-2026-11-23/transients-goods5324.py \
  --cache /path/to/scratch/goods5324-cache --out /path/to/scratch/goods5324-replay
python docs/followup-2026-11-23/transients-goods5324-controls.py \
  --cache /path/to/scratch/goods5324-controls-cache --out /path/to/scratch/goods5324-controls
```

The last command downloads only the specified ranges if absent from cache and rejects bytes that disagree with the frozen receipt. No telescope request, scientist message, TNS report or paper submission was sent in this review. The existing [infrared follow-up message](../../papers/rnaas/messages/03-flare-followup-email.md) can use this packet after the observer-facing claims and current target state are reviewed.
