# Mutual-orbit families for (556068) 2014 JR80

The combined data favor a **103.03-day mutual orbit**, with model-dependent eccentricity **0.289**, relative semimajor axis **22,870 km** and total system mass **8.93 × 10¹⁹ kg**. This is an inference for a previously known binary. Sparse sampling permits aliases; the next fitted family has a 154.62-day period. Johnston's archive already listed a 150-day **estimate**, based on presumed physical properties; the 2026 announcement gives projected separations without an orbital solution. [CBET 5739](https://www.cbat.eps.harvard.edu/iau/cbet/005700/CBET005739.txt) · [Johnston's entry and estimate flag](https://johnstonsarchive.net/astro/astmoons/am-556068.html)

## Data and model

Thirteen independent observing nights span 2020–2026. Nine nights come from CFHT/MegaCam, one from HST/WFC3 and three from new Rubin measurements. The nine CFHT night IDs and contributing exposures are enumerated in `cfht_ids.json`; every accepted archived frame is labelled in the accompanying archive analysis. We use median east/north components per CFHT night after removing the 180° component-label ambiguity within each night, average the four HST exposures, and average Rubin bands within each night. The final inputs exclude background-contaminated 2501872p. The earlier frozen prediction used one preliminary exposure per CFHT night and included that frame; see `prediction_training.csv`. These training sets are explicitly separate. [Archived exposures](../archive/measurements.csv) · [Final 13 vectors](observations.csv)

Orbit fitting uses reception times at exposure midpoints (Rubin MJD TAI minus 37 seconds), removes differential light travel time and rotates a common three-dimensional Keplerian ellipse into each epoch's local east–north frame. Relative astrometry labels the components only up to an overall reflection through the origin. The grid fits all quadrant combinations between separated observing blocks. A separate run lets each Rubin night independently reverse its labels. Trial periods span 14–1,100 days at 0.000012 cycles/day spacing; trial eccentricity spans 0–0.9. Phase and four linear sky-plane coefficients are profiled at each grid point. The 200 best separated grid minima are then refined nonlinearly, testing both mirrored three-dimensional orientations and allowing eccentricity to 0.94. This finite search cannot prove global uniqueness. The CSV node angle is measured from the reference epoch's east direction toward north, not conventional astronomical position angle.

The HST per-coordinate error floor is 0.006″. Ground-based nights receive 0.04″, except 0.05″ for the most weakly resolved 2021 CFHT epoch. The best-fitting sum of squared normalized residuals is 2.26 with 19 nominal continuous degrees of freedom. Residuals this small are consistent with conservative floors and/or correlated systematics. The fitted mass follows Kepler's third law from the relative semimajor axis and period. The uncertainties below come from local least-squares curvature with the assigned floors held fixed; they include neither the discrete aliases nor absolute plate-calibration errors.

| Period family | χ² | Eccentricity | Semimajor axis | Total mass |
| --- | ---: | ---: | ---: | ---: |
| 103.026 d | 2.26 | 0.289 | 22,873 km | 8.93 × 10¹⁹ kg |
| 154.623 d | 17.11 | 0.432 | 27,736 km | 7.07 × 10¹⁹ kg |
| 309.451 d | 34.59 | 0.645 | 44,677 km | 7.38 × 10¹⁹ kg |
| 77.283 d | 56.32 | 0.531 | 27,604 km | 2.79 × 10²⁰ kg |

Within the 103-day family the nominal conditional standard errors are 0.029 d, 0.019, 608 km and 0.71 × 10¹⁹ kg respectively. Inflating all ground errors by 75% gives 0.050 d, 0.034, 1,060 km and 1.24 × 10¹⁹ kg. The 103-to-155-day Δχ² then decreases from 14.86 to 4.95. The leading family also survives independent component signs per Rubin night, equalization of the total weight per observing block, and replacement of the pipeline-PSF science fits by difference-image or empirical-stellar-PSF fits. [Alias search](fit_all.csv) · [Robustness and mass propagation](sensitivity.csv)

## What the Rubin data test

Holding all three Rubin nights out of a new period search leaves 103.025 d, 154.689 d and 309.659 d as the best three families. With no Rubin vector used to fit them, they give these observable errors, each the smaller of the residual lengths for the two possible component orderings. Thus Rubin tests an archive-only prediction separately on each night.

| Predicted from CFHT + HST | June 27 Rubin | July 7 Rubin | July 9 Rubin |
| --- | ---: | ---: | ---: |
| 103 d | 0.035″ | 0.017″ | 0.012″ |
| 155 d | 0.072″ | 0.056″ | 0.090″ |
| 310 d | 0.104″ | 0.106″ | 0.146″ |

The archive-only Δχ² between the first two is 10.11; adding Rubin raises it to 14.86. Data source selection for this investigation used Rubin morphology, so this leave-out calculation is a conditional cross-validation, not a prospective discovery sample. Withholding the HST night instead predicts its measured vector within 0.0085″. Removing MJD 59403, where target trailing could affect the most compact CFHT pair, predicts this held-out vector within 0.032″ for the 103-day family and 0.772″ for the 155-day one. [All Rubin nights held out](fit_no_rubin.csv) · [HST held out](fit_HST.csv) · [Late 2021 CFHT held out](fit_CFHT_59403.csv)

Whole-season tests illustrate an important ambiguity. Without 2020, periods 97.60, 92.72 and 103.03 d are almost tied in training χ² (2.11, 2.17, 2.21); the held-out 2020 residuals span 0.46–0.54″, 0.67–0.73″ and 0.004–0.009″ respectively. Without 2021, the training χ² is 0.86, 0.89 and 0.90 for 103.03, 108.17 and 98.35 d. Only the first of these reproduces all six held-out vectors within 0.042″. Neither season alone selects a unique cycle count. [2020 held out](fit_2020.csv) · [2021 held out](fit_2021.csv)

## Predictions and checks

The first forecast was frozen at [65252389](https://github.com/dovvnloading/Open-Space/commit/65252389d4c0f98cc20ca987f1d92e5d04816310) before the archival analyst supplied DECam measurements; the second batch at [3dec65cb](https://github.com/dovvnloading/Open-Space/commit/3dec65cb3dda035d80eee31faf834d70f3843ccc). Additional frame forecasts are preserved in `second_predictions.csv` and `last_predictions.csv`. They all use the preliminary CFHT fit. A model agreement in an off-detector footprint, noisy single-source image, or failed two-source fit supplies no confirmation. In the first batch the two on-chip April 25, 2024 frames and the March 2025 frame provided no admissible vector; the first three deep 2022–23 frames missed the usable detector area.

The August 6, 2022 DECam i-band exposure supplies a blind archival test outside the 13-night orbit fit. The preliminary 103-day model predicted **0.592″ at 146.56°**. An independent image reduction measured **0.660″ at 147.85°**, with a nonparametric moments angle of 148.2°. The measured vector is 0.069″ from the prediction. The 155-, 309- and 77-day predictions miss by 0.587″, 1.668″ and 0.197″ respectively. The secondary/primary ratio hit the fitting lower bound, so this test is stronger in orientation than separation. The target's pair-fit improvement is 40.40, versus at most 3.21 in six comparison stars. Predictions were hashed and communicated before the image analyst supplied this vector. The archive analyst knew the candidate period scale but had not received the predicted vectors. An existing CBET statement that binary structure was detected in 2019–2023 might encompass this image; CBET does not provide its measured vector. [Blind test residuals](blind_residuals.csv) · [Image and comparison-star measurements](../archive/README.md)

A second admitted frame, June 2, 2019 (MJD 58636.35649), independently gives **0.654″ at 30.92°**. The frozen 103-day prediction, **0.667″ at 29.05°**, is 0.026″ away; the 155- and 309-day predictions miss by 0.622″ and 2.164″. An aperture-moments measurement gives 30.71°. This second observation is a different year, sky field, and detector pointing. A third independently accepted frame from May 19, 2023 gives **0.660″ at 42.09°**, compared with the frozen 103-day **0.618″ at 46.19°**, a 0.062″ vector residual. Its alternative moments axis is 39.04°. The 155- and 309-day predictions miss by 1.241″ and 2.025″. All three DECam vectors stay outside the orbit fit. The release preserves all five forecast families at every admitted epoch.

Two future measurements around opposition can distinguish the final 103- and 155-day families. Dates were chosen to maximize sign-independent vector differences subject to two conditions at 05:00 UTC as seen from Rubin: solar elongation exceeds 160° and target elevation exceeds 60°. The two dates must be at least ten days apart. These are model-conditioned predictions; position angle is east of north, modulo 180°. [Full coordinates, ranges, and observability](forecast.csv) · [JPL Horizons response](opposition_horizons.txt)

| UTC date at 05:00 | 103 d predicted separation and PA | 155 d predicted separation and PA |
| --- | ---: | ---: |
| 2027 Aug 6 | 0.751″, 31.47° | 1.369″, 148.53° |
| 2027 Aug 16 | 0.612″, 64.66° | 1.332″, 156.04° |

## Reproduce

Run from this directory with Python 3.10+, NumPy and SciPy. The independent-quadrant scan uses more memory; a single BLAS thread reduces overhead.

```sh
OPENBLAS_NUM_THREADS=1 python fit.py
OPENBLAS_NUM_THREADS=1 python fit.py --independent-rubin
for omitted in 2020 2021 no_rubin HST CFHT_59403 2026-07-09; do
    OPENBLAS_NUM_THREADS=1 python fit.py --omit "$omitted"
done
python sensitivity.py
python historical_predictions.py
python validate_blind.py
python forecast.py
```

`opposition_horizons.url` gives the exact forward ephemeris request. The archived Horizons responses for CFHT and blind DECam frames contain the observer code, target and UTC epoch list. All angles in `fit_*.csv` parameter columns are radians, except where degrees are explicitly named. Repeated rows in these files can represent separate grid starts converging on the same local solution.
