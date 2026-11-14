# Predicting the next reported spin glitch in PSR J1341−6220

14 November 2026

A model fitted to 33 historical intervals predicted four later Parkes intervals with a mean absolute error of **82.52 days**. Always predicting the historical median, 263 days, gave **131.75 days**. Glitch size reduced the mean error by 49.23 days (37.36%). The gain was concentrated in the longest interval, for which the model still predicted the next glitch about 218 days too early.

Fuentes et al. had already reported a correlation between glitch size and the following wait for this pulsar, also called B1338−62. Shang et al. discovered the seven later glitches. This analysis measures how well a rule fitted to the older observations predicted their spacing. The target was selected because a later series was available. This was a retrospective test in disjoint observing periods, and some historical measurements also came from Parkes. [Fuentes et al., Table 3](https://doi.org/10.1051/0004-6361/201935939) · [Shang et al., Table 2](https://doi.org/10.1051/0004-6361/202661329)

## Four scored Parkes intervals

Each forecast uses only the size of the first glitch in the pair. Size is the fractional spin-frequency step in units of 10⁻⁹. Absolute errors below use the midpoint epochs from Table 2. A positive gain means the model had the lower absolute error.

| Events | Size | Observed wait (days) | Forecast (days) | Model error (days) | Median error (days) | Gain (days) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| G36 → G37 | 136.4 | 149 | 224.79 | 75.79 | 114 | 38.21 |
| G39 → G40 | 705.1 | 331 | 319.86 | 11.14 | 68 | 56.86 |
| G40 → G41 | 385.7 | 292 | 266.46 | 25.54 | 29 | 3.46 |
| G41 → G42 | 953.4 | 579 | 361.37 | 217.63 | 316 | 98.37 |

Parkes observed from MJD 58222 through 60589. Table 2 reports seven glitches in order. G37 and G38 lack a measured size, which removes their outgoing forecasts. G37 still closes the forecast made at G36. The observers could not maintain a continuous timing solution between G37 and G38 and discuss possible intervening small glitches. G42 has no successor in the series. The target quantity is days to the next **reported** glitch.

The catalogue series ends at MJD 58214, and its last eligible training interval ends at MJD 58178. The size model uses a Theil–Sen slope, the median of the slopes between historical points of unequal size. Constrain the slope to be nonnegative and take the median fitted residual as the intercept. The resulting forecast is 201.98469 + 0.16717695 × size, in days, with a zero floor. Both the coefficients and the median baseline use only eligible historical intervals.

Removing any one scored pair leaves an average gain of at least 32.85 days on the other three. Under exchangeable waits with the four predictions fixed, only one of the 24 possible assignments matches or improves the observed model error. The dates come from one pulsar and adjacent intervals can be dependent. This conditional permutation result gives limited evidence about performance across future events.

Allowing both endpoint dates to span their published observing brackets and the starting size to vary by one quoted standard error gives outer bounds of **30.02–49.27 days** on the mean gain for these four intervals. The winner on either middle pair can change. The bounds hold historical coefficients fixed and allow a shared epoch to move separately for adjacent pairs. They describe measurement sensitivity in this sample. Future-event uncertainty also depends on the training sample, pulsar behavior and observing coverage.

## Historical catalogue screen

The 4 November Jodrell Bank table, retrieved on 14 November 2026, contained 737 glitches in 226 sources. Of 511 adjacent catalogue pairs, 483 remained eligible after checking starting sizes, epoch ordering, intervals containing conflicting events and known observational breaks. The audit used the ATNF compilation and published timing studies. A second catalogue may cite the same measurements. For example, Xia et al. measured two glitches between Jodrell Bank's consecutive J2229+6114 epochs 59262 and 61071, so that interval was excluded. [Jodrell Bank](https://www.jb.man.ac.uk/pulsar/glitches.html) · [ATNF](https://www.atnf.csiro.au/people/pulsar/psrcat/glitchTbl.html) · [Xia et al., Table 1](https://arxiv.org/html/2605.13011v1#S3.T1)

Seven sources had at least eight training intervals and five chronological test intervals. For the known J0537−6910 control, train on RXTE and test on NICER, dropping the multi-year gap. For each exploratory source, test on the latest third of eligible intervals, rounded up, with a minimum of five, and fit on the earlier intervals. The model and the constant median were fixed for each test. The J1341−6220 catalogue test below uses 22 training intervals. The later Parkes test above uses all 33 eligible intervals in the older catalogue. [J0537−6910 timing study](https://doi.org/10.1093/mnras/staf2264)

| Pulsar | Train | Test | Model MAE (days) | Median MAE (days) | MAE reduction | Permutation p | Holm p |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| J0537−6910 (control) | 44 | 22 | 16.32 | 54.32 | 69.95% | 0.00005 | — |
| J0631+1036 | 10 | 6 | 256.76 | 548.83 | 53.22% | 0.10000 | 0.5000 |
| J1341−6220 | 22 | 11 | 96.70 | 124.05 | 22.04% | 0.00725 | 0.0435 |
| J1740−3015 | 22 | 12 | 348.15 | 374.75 | 7.10% | 0.14465 | 0.5786 |
| J0534+2200 | 16 | 9 | 607.80 | 607.80 | 0.00% | 1.00000 | 1.0000 |
| J0835−4510 | 16 | 9 | 323.69 | 313.51 | −3.25% | 0.49925 | 0.9985 |
| J1413−6141 | 8 | 5 | 77.05 | 66.63 | −15.64% | 0.30000 | 0.9000 |

The exploratory follow-up criteria required a positive fitted slope, at least 25% lower MAE, positive average gain after deleting any one test pair, and Holm-adjusted p ≤ 0.05. None of the six exploratory pulsars passed all four. The control is a previously published relation and sits outside the six-source multiplicity correction. The permutation statistic holds forecasts fixed and permutes waits. It enumerates permutations up to eight test points and uses 19,999 seeded permutations plus the observed assignment for larger test sets. Serial dependence and observing-cadence changes limit that test.

## Reproduction

Run from the repository root with Python 3.10 or later. The verifier needs only the standard library. It fetches the original Jodrell Bank HTML, checks its SHA-256 digest, looks up the selected events by row ID, reconstructs waits from consecutive epochs of each pulsar, and refits both estimators. It also checks the seven historical scores, six Holm results and external sensitivity bounds.

```sh
python results/pulsar-glitch-forecast/verify.py
```

The live catalogue can change. For exact replay after a revision, use a saved copy with SHA-256 `e811419f2037080bd2c99c1db30555dd4ee6f462f36b91796908f0fe90aabeb4`.

```sh
python results/pulsar-glitch-forecast/verify.py --catalog-html /path/to/original-gTable.html
```

The numeric IDs in the CSV select the starting event of each eligible interval. The successor is always the immediately following event at the same source, including events with an unavailable size. Together, J1341−6220's training and catalogue-test IDs define its 33-event-pair external training set. The JSON contains all seven later events, inclusion decisions, per-pair results, and the paper PDF fingerprint. [Historical screen and row IDs](screen.csv) · [Parkes measurements and predictions](parkes-intervals.json) · [Verifier](verify.py)

The Jodrell Bank compilation credits Basu et al. and Espinoza et al. External dates, sizes and uncertainties are a transcription of selected columns of Table 2 in L. H. Shang, Y. L. Liang, Y. Z. Qiao, P. Liu and W. J. Xiao, *Astronomy & Astrophysics* 715, A132 (2026). The authors published the article under CC BY 4.0. Waiting times, forecasts, error comparisons and sensitivity bounds were calculated for this analysis. [Basu et al.](https://doi.org/10.1093/mnras/stab3336) · [Espinoza et al.](https://doi.org/10.1111/j.1365-2966.2011.18503.x) · [Shang et al.](https://doi.org/10.1051/0004-6361/202661329) · [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
