# Can Gaia predict an omitted season of eclipse timing?

TIC 96100077 has 38 usable three-day eclipse-clock measurements in five TESS sectors. The last two sectors are consecutive in 2026, so we excluded them together. Each of the four folds withholds an entire season, fits all timing coefficients on the other three, and predicts the excluded season at its original timestamps. A separate forward check trains on the 2019, 2021 and 2023 seasons and forecasts 2026. These are retrospective tests of the existing measurements. No new TESS visit entered this validation.

## Absolute timing

Each cell is a root mean square timing error in seconds. The all-season score gives each withheld season equal weight. The comparison models use the same measured times throughout. The Gaia model uses the published outer period, eccentricity, orientation and periastron epoch, with one fitted signed amplitude. There is no test-season intercept in these absolute predictions.

| Model | Four excluded seasons, equal weight | 2026 forecast from earlier seasons |
| --- | ---: | ---: |
| Constant inner period | 171.9 | 265.6 |
| Constant inner period + Gaia orbit | 266.0 | 296.6 |
| Slowly changing inner period (quadratic ephemeris) | 496.4 | 828.6 |
| Quadratic ephemeris + Gaia orbit | 135.5 | 159.5 |
| Constant period + free annual sine and cosine | 150.9 | 205.7 |
| Constant period + free semiannual sine and cosine | 148.8 | 198.9 |
| Quadratic ephemeris + free sine and cosine at Gaia's period | 618.6 | 1127.6 |

Adding Gaia to the quadratic ephemeris improves each of the four omitted-season forecasts. That baseline alone extrapolates poorly. The gain over the simple constant-period model is smaller. Adding the orbit on top of a constant-period ephemeris makes the absolute forecasts worse. In the quadratic analysis, the individual primary and secondary clocks both improve (566 to 143 s and 444 to 202 s), while their difference gets worse (123 to 146 s). A common light-path signal would cancel in that difference. Local waveform changes, shared detrending, spot modulation and the limited number of observing seasons can still move these clocks. Their fitted amplitudes have no direct mass interpretation here.

## Timing changes within visits

The four common-clock slopes are −7.94, +3.38, +5.10 and +3.08 seconds per day. Gaia predicts +9.81, −2.98, −2.73 and −2.05 seconds per day for its RV primary. The observed sign reversal is consistent with the eclipse source lying on the other side of the outer barycenter. Three predicted slopes are similar. Two of the 23 other permutations of the slopes fit better than their original assignment. This sector-level comparison preserves correlations inside each season and allows a common slope and free signed scaling. It cannot estimate a population false-alarm probability.

A separate test predicts the *change* in timing after reading only the first timing window of an omitted season. Each training season has its own fitted epoch. On 2026 data, adding the Gaia term to a linear clock reduces the error from 95.1 to 24.7 s using this first-window anchor. The analogous four-season score changes from 83.4 to 75.9 s. An additional centered-shape diagnostic subtracts the whole held season's mean from its measurements and model. That diagnostic cannot assess absolute predictive performance.

## Phase, period and ephemeris uncertainty

Of the 99 alternative shifts of the outer phase, 16 have smaller absolute quadratic cross-validation error than the nominal ephemeris. A shift of −0.09 cycle gives 33.2 s. Choosing that shift from this grid uses all four seasons. It is an alternative fit, not an independent estimate of the true outer phase. On a separate 120–600 d grid at fixed periastron epoch and fixed shape, two of 241 trial periods (198 and 462 days) do better than the published period. The choices of phase origin and comparison grid matter. Neither count has a calibrated probability interpretation.

We drew 2,000 joint realizations of C, H, eccentricity, period and periastron epoch from the catalog's formal multivariate-normal approximation. The covariance is the corresponding marginal block of Gaia's 15-parameter AstroSpectroSB1 correlation matrix, multiplied by its parameter errors. Times use TCB in the orbital model and TDB in the TESS table. The middle 95% of four-season absolute RMS errors span 64–251 s for the quadratic-plus-Gaia model and 215–294 s for the linear-plus-Gaia model. Each draw refits all timing coefficients on training seasons only. These ranges propagate the published orbital covariance conditional on the measured eclipse times. They do not include a physical model of starspots, orbital model systematics or uncertainty in the choice of timing baseline.

The covariance gives approximate central 95% windows of October 26–December 4, 2026, and July 24–September 4, 2027 (UTC), for two predicted periastra. A spectrum obtained near an expected periastron can constrain the 273-day orbit. Establishing which velocity belongs to the eclipsing pair still requires additional component information.

## Replay

The inputs are the original investigation's committed `outputs/timings.csv` and `data/gaia-orbit.csv` at Open-Space commit 062d4c61aa5b41608a61cf6aed9f82331b548f30. The source directory is an explicit argument. From this directory, point the command below to that source or an equivalent checkout. Python 3.12, NumPy, pandas and Astropy were used. All result tables are regenerated locally, including every fold, phase trial, period trial and slope permutation.

```sh
OPENBLAS_NUM_THREADS=1 python validate.py /path/to/Open-Space/results/tic-96100077
sha256sum -c SHA256SUMS
```

[Input measurements and extraction](https://github.com/dovvnloading/Open-Space/tree/062d4c61aa5b41608a61cf6aed9f82331b548f30/results/tic-96100077) · [Gaia covariance order and orbital conventions](https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_non--single_stars_tables/ssec_dm_nss_two_body_orbit.html)
