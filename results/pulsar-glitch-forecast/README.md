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

## What could distinguish the two forecasts?

A post-hoc comparison checks whether the four gains can be explained by a single shift in typical waiting time. Even a constant chosen after seeing all four waits gives a minimum MAE of 117.25 days, compared with 82.52 days for the historical size-based rule. This only rules out a single constant as an equally good fit to these four midpoints. Secular changes in cadence or glitch behavior are still possible.

For G42, the fitted rule gives 257.04 days after MJD 60422, or MJD 60679.04 (4 January 2025). The historical median gives MJD 60685. The study's last exposure was MJD 60589, 90 days before the size forecast, and the NASA 3PC timing solution ends at MJD 57277.06. We did not verify a successor in the available follow-up. Shang et al. already estimated 222–520 days for this wait from the post-glitch recovery, a different input. The six-day separation between the size and median dates is smaller than the two-to-four-week Parkes cadence described in the paper. An event with a measured step around 950 × 10⁻⁹ would separate the frozen forecasts by roughly 98 days. A prospective test could lock both dates at detection and score the subsequent independently timed event. [Parkes timing study, Section 4.2](https://doi.org/10.1051/0004-6361/202661329) · [NASA 3PC ephemeris](https://fermi.gsfc.nasa.gov/ssc/data/access/lat/3rd_PSR_catalog/3PC_HTML/J1341-6220.html)

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

Run from the repository root using Python 3.10 or later. The full-population verifier fetches both catalogues and checks their SHA-256 hashes. It orders Jodrell Bank events before applying exclusions, rebuilds all eligible source populations and their chronological partitions, compares those partitions to the published event IDs, refits the estimators and recomputes the errors and permutation results. The primary screen uses the thresholds described above. The four-interval uncertainty bounds and best-constant comparison are additional, retrospective diagnostics. Numerical reproduction cannot establish when the model was chosen or the events were first examined.

```sh
python results/pulsar-glitch-forecast/verify.py
```

The live catalogues can change. For replay after a revision, supply copies matching these fingerprints. The original catalogue HTML is required to verify the full population. A later release can change event IDs or introduce a glitch between two historical neighbours.

| Snapshot | SHA-256 |
| --- | --- |
| Jodrell Bank | `e811419f2037080bd2c99c1db30555dd4ee6f462f36b91796908f0fe90aabeb4` |
| ATNF | `4e605281a11fb36457a891cf8b0952d7e59f3abc5c91544dbe364f0dbc748c51` |
| Shang et al., October 2026 proof | `dad1b82c08c2179da6221fceee8a8b0c946619282c3f96ac36555ea67dcd8a36` |

The default calculation uses the transcribed Parkes data. With a copy of the original proof and Poppler's `pdftotext` on PATH, the command below also checks every transcribed epoch, size, missing value and reported uncertainty against Table 2 on page 4. The public evidence includes the per-interval bounds and exclusions. [Original proof](https://pdfs.semanticscholar.org/0a94/a62fa71a5679b4a0cedcdafe1e5af9d4fce4.pdf) · [Screen and event IDs](screen.csv) · [Parkes transcription and results](parkes-intervals.json) · [Verifier](verify.py)

```sh
python results/pulsar-glitch-forecast/verify.py \
  --catalog-html /path/to/jodrell.html --atnf-html /path/to/atnf.html \
  --paper-pdf /path/to/shang2026.pdf
python -m unittest discover -s results/pulsar-glitch-forecast -p 'test_*.py'
```

The ATNF check matches J-names and accepts epoch differences within the larger of three days or twice the sum of reported uncertainties. The floor rises to 30 days when either epoch uncertainty is missing. When an ATNF event has several plausible matches, the verifier leaves the identity unresolved and retains those intervals unless another rule excludes them. Separate exclusions cover early low-cadence Crab monitoring, a disputed Crab event near MJD 50489.7, the RXTE–NICER gap and two intervening CHIME events. A retained pair may still contain missed glitches. The replay applies the cited literature decisions and cannot establish complete exposure or guarantee the physical event sequence. [Crab cadence](https://academic.oup.com/mnras/article/482/3/3736/5142868) · [Disputed Crab event](https://doi.org/10.1093/mnras/sty1294) · [RXTE/NICER observing spans](https://doi.org/10.1093/mnras/staf2264) · [CHIME timings](https://arxiv.org/html/2605.13011v1#S3.T1)

The Jodrell Bank compilation credits Basu et al. and Espinoza et al. External dates, sizes and uncertainties are a transcription of selected columns of Table 2 in L. H. Shang, Y. L. Liang, Y. Z. Qiao, P. Liu and W. J. Xiao, *Astronomy & Astrophysics* 715, A132 (2026). The authors published the article under CC BY 4.0. Waiting times, forecasts, error comparisons and sensitivity bounds were calculated for this analysis. [Basu et al.](https://doi.org/10.1093/mnras/stab3336) · [Espinoza et al.](https://doi.org/10.1111/j.1365-2966.2011.18503.x) · [Shang et al.](https://doi.org/10.1051/0004-6361/202661329) · [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
