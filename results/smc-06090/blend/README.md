# OGLE-SMC-ECL-06090: Can two unchanging eclipse clocks explain the timing reversal?

An otherwise stationary blended pair of eclipsing sources can alter the observed center of an overlapping dip. This test gives that explanation four independently adjustable eclipse profiles and compares its predictions with a model whose two eclipse centers move in opposite directions.

## Data and fitting

All 637 OGLE-III and 711 OGLE-IV I-band measurements pass the criteria of finite time, magnitude and uncertainty with 0 < uncertainty ≤ 0.10 mag. The catalog reference period is 1.1696719 days, primary epoch HJD 2460000.0369. OGLE-III dates have 2450000 added. Source URLs and SHA-256 checksums are enforced by compare.py. The source catalogs are [Pawlak et al. 2013](https://arxiv.org/abs/1310.3272) and [Głowacki et al. 2025](https://arxiv.org/abs/2503.15596).

Convert magnitudes to linear flux relative to I = 18.23. Propagate the magnitude errors with an additional 0.012 mag floor in quadrature. Every model includes two freely fitted survey baseline levels and the same weak cosine variation at twice the catalog orbital frequency. All eclipse deficits add in flux and have a wrapped Gaussian shape. Use the sum of twice the soft-L1 objective, with standardized residuals and scale 1.8, for training and evaluation. This is a robust prediction score. It is not a chi-squared statistic.

| Model | Parameters | Eclipse centers in the folded light curve |
| --- | ---: | --- |
| Common constant period | 10 | Two dips share a fitted linear phase drift, with distinct centers, depths and widths. |
| Two independent constant periods | 11 | Two dips have independent linear phase drifts. Each can belong to a different source with an undetected secondary eclipse. |
| Two blended binaries | 17 | Each of two independent periods carries its own two dips. Each of the four has separately fitted central phase, constant depth and constant width, so the sum can change shape. |
| Oppositely moving centers | 13 | The two dip centers share a linear phase drift plus equal and opposite sinusoidal displacements. Constant depth and width of each eclipse are fitted. The sinusoid describes the leading timing signature of apsidal motion. |

The contemporaneous zero point is HJD 2456000. A trial orbital phase drift of s cycles/year corresponds to a trial frequency 1/P_catalog − s/365.25. The baseline blend optimization allows up to 0.055 and 0.11 cycle/year of drift in its two clocks and permits depths to vanish. Their epoch centers may move within specified basins for the associated primary and secondary; see setup(). The optimization is nonlinear and cannot certify that every alternative period and every possible shape have been excluded. The blend's globally best model among the starts had periods 1.16966840 and 1.16969013 days. The extra model parameters do not rescue its predictive accuracy.

## Unequal-epoch, held-block predictions

Each run withheld one contiguous time block, fitted to all remaining observations afresh, and evaluated every withheld flux. A point never predicted its own held-block fit. Five of the six predictions can use data both before and/or after the block. The last block is a chronological forecast using only earlier photometry. The data were explored during target selection; these are retrospective validation runs.

| Held-out OGLE block, chronological | Measurements | Common clock | Separate clocks | Two binaries | Opposite motion |
| --- | ---: | ---: | ---: | ---: | ---: |
| 2001–early 2005 | 319 | 471.5 | 1278.3 | 1251.5 | 240.4 |
| Mid-2005–2009 | 318 | 533.5 | 399.9 | 539.8 | 339.7 |
| 2010–2014 | 382 | 1289.8 | 1306.0 | 1254.8 | 479.9 |
| Late 2014–mid-2018 | 147 | 328.1 | 198.9 | 305.7 | 172.1 |
| Mid-2018–early 2020 | 86 | 403.7 | 258.7 | 316.1 | 114.6 |
| 2022–early 2024 | 96 | 574.8 | 423.7 | 509.8 | 146.7 |
| Total robust score | 1348 | 3601.3 | 3865.5 | 4177.6 | 1493.5 |
| Pooled RMS flux error | | 0.05247 | 0.06263 | 0.06450 | 0.03085 |

The six date blocks use boundaries HJD 2453500, 2455200, 2457000, 2458300 and 2459700. All four models have a baseline fitted for each survey even when an individual block from that survey is held out. In the full fit the opposite-motion model has a timescale of 27.29 years. This is an empirical sinusoid fitted across about 23 years, rather than an independently observed full cycle.

Grouping the already predicted errors by observing seasons separated by >60-day gaps yields 20 groups. The opposite-motion model beats the two-binary blend in all 20. Its mean robust loss per test observation is 1.99 units lower, with a seasonal resampling range of 1.42–2.59 (2.5th–97.5th percentiles, 10,000 replicates). This describes the season-to-season spread conditional on the cross-validation fits. The overlapping training sets preclude an independent-trial significance interpretation.

## Alias and optimization checks

- For the separate-clock model, independently allow each of the two orbital frequencies to lie near the catalog frequency, its +/−1 per sidereal day alternatives, or its +/−1 per year alternatives. Refit all 25 frequency-pair choices with local period corrections and six starting points per choice on every training set. The catalog-frequency neighborhood won on training data for the full fit and every cross-validation fold. The best non-catalog pair in the full dataset scored 4,375 vs. 2,495 for the catalog-neighborhood pair. Negative formal frequency in the lower daily sideband is equivalent to a positive period with reversed phase for a symmetric dip.
- The 4-dip blend used 96 starting points in the full fit, 36 for each held block. Another 160 starts each expanded the phase/rate bounds and challenged the optima for the full fit and chronological forecast. All four independent batches of 40 recovered the same best full-data score within 0.01 and terminal-block prediction around 510.
- The score replicates across seasons with different cadences. It disfavors a two-source blend of the type tested. Evolution of a source's own eclipse shape or brightness, a more complex multiple system, instrumental issues unique to this source and photometrically unresolved sources with very different periods remain possible. Physical interpretation requires additional data.

## Replay

Run in this directory with Python 3, NumPy and SciPy. These commands use the supplied primary OGLE files and make no external requests.

```sh
OPENBLAS_NUM_THREADS=1 python compare.py --model basic --output basic.json
OPENBLAS_NUM_THREADS=1 python compare.py --model blend --output blend.json
OPENBLAS_NUM_THREADS=1 python aliases.py
OPENBLAS_NUM_THREADS=1 python stress.py
python summarize.py
```

compare.py and aliases.py write all training fits. summarize.py reconstructs held-block predictions from their parameters to calculate the reported scores and seasonal intervals. Photometric data credit the Optical Gravitational Lensing Experiment collaboration.
