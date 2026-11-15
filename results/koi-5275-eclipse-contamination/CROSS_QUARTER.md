# Predicting the eclipse depth in the other quarter

The binary's light curve predicts a dip of **77.0 parts per million** in the KOI-5275 target during quarter 14 and **77.7 parts per million** in quarter 15. The measured dips in the same windows are 77.6 and 79.0 ppm. We fit the amount of contaminating light in one quarter and carried it into the other. The experiment tests how much of the shallow aperture signal this neighboring source can explain.

| Quarter measured | Quarter used for dilution | Complete orbits | Predicted dip | Measured dip | Residual and 95% interval |
|---|---|---:|---:|---:|---:|
| 14 | 15 | 23 | 77.0 ppm | 77.6 ppm | +0.6 ppm (−45.1 to +44.8) |
| 15 | 14 | 23 | 77.7 ppm | 79.0 ppm | +1.3 ppm (−45.4 to +50.7) |

The point estimates account for the measured depth. The error intervals allow an additional unresolved signal of several tens of ppm. The dilutions come from simultaneous light curves extracted for the binary KIC 6367628 and the candidate target KIC 6367663. Source light reaches the target stamp from outside its boundary. The published DR25 source-probability model favored the target star with probability 0.9959. The independent pixel gradients and simultaneous flux variations support a specific off-target contribution.

We used the Gaia period of 3.77990376255 days. The Q14 *source* data alone set the eclipse phase by the faintest of 100 equal-width phase bins. We fixed that phase for both targets and quarters. We joined data by spacecraft cadence number and selected unflagged finite measurements. We retained orbits with at least two samples within 0.047 cycle of minimum and four baseline samples between 0.15 and 0.35 cycle away on either side. Within each quarter we regressed the orbit-demeaned target flux on the orbit-demeaned binary flux. We predicted the held-out target's depth by multiplying its simultaneous binary depth by the opposite quarter's fitted dilution. The intervals independently resample the 23 orbits in each quarter with replacement 3,000 times. They are conditional on the chosen orbital phase and extraction.

Run `python dilution_test.py` with NumPy and Astropy. The script downloads the four public Kepler light curves. Pass `--cache DIR` to use existing copies. The [results table](cross_quarter_test.csv) preserves unrounded values. [Kepler data products](https://archive.stsci.edu/missions-and-data/kepler)
