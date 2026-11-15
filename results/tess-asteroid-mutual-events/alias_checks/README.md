# Testing detector position, sky contamination and a cadence alias

All tests use the cleaned TSSYS-DR1 exposure rows, the event windows in the parent measurements table, and the Pál et al. catalog primary-spin periods as starting values. Both asteroids crossed Sector 3, Camera 1, CCD 3. The data columns follow the release documentation. X and Y are predicted CCD positions. The sky level and its scatter refer to the background on differential images.

| Test | (29489) | (3177) |
|---|---:|---:|
| Original median event depth (mag) | 0.1282 | 0.0477 |
| Depth predicted from fractional pixel position and sky level/RMS (mag) | −0.0008 | +0.0007 |
| Depth after subtracting that prediction (mag) | 0.1272 | 0.0437 |
| Positive events after subtraction | 50/50 | 31/32 |
| Depth after removing the highest quarter of sky RMS (mag) | 0.1282 | 0.0493 |
| Positive events after removing the highest quarter of sky RMS | 40/40 | 22/22 |
| Positive events after removing the most extreme quarter of differential sky levels | 43/44 | 26/26 |
| Depth after subtracting 6 / 14 / 28 primary-spin harmonics (mag) | 0.1285 / 0.1252 / 0.1263 | 0.0476 / 0.0262 / 0.0180 |
| Nearest spin-harmonic/cadence alias among harmonics 1–32 (cycles/day away from event frequency) | 0.48888 | 0.01451 |
| Frequency resolution, 1 / time span (cycles/day) | 0.03661 | 0.03670 |

For the detector-position check we regress the target residuals against a two-dimensional Fourier basis on the fractional X and Y pixel positions. The 28 terms include both single-coordinate harmonics and mixed directions. Sky level and RMS enter as fractional ranks, their squares and a product. We omit the central 54% of each event cycle from training. Each prediction uses six other interleaved groups of complete cycles, so the corrected event has contributed no photometry to its own model. The analysis measures the median difference between in-event and away-from-event medians for each qualifying cycle. The prediction and correction medians need not add to the original because medians are nonlinear. Both event trains span all eight 1/8-cadence bins. A localized stationary contaminant, an unmodeled finer pixel pattern or a problem in individual difference images could survive these tests.

**The primary-spin alias materially weakens the interpretation for (3177).** Its fitted primary frequency of 6.76581 cycles/day has a 14th harmonic at 94.7214 cycles/day. Sampled at 48 times per day, that harmonic is indistinguishable on a perfectly regular grid from 1.27860 cycles/day, close to the measured 1.29311 cycles/day recurrence. The 0.01451 cycles/day gap is smaller than this observing run's frequency resolution. On the actual time stamps, reoptimizing the primary frequency with 28 freely fitted spin harmonics absorbs about 62% of the original median event depth. Thirty of the 32 events still dim, at a median 0.0180 mag. A flexible spin model can absorb a genuine second signal near an alias as well. A new observation at a different sampling cadence would separate a persistent low-frequency clock from an alias moving with the sampling frequency.

For (29489), no alias through the 32nd spin harmonic is within 13 frequency-resolution elements of the observed recurrence. Adding up to 28 harmonics changes the measured depth by less than 0.004 mag. The high-order fits are deliberate stress tests of a single primary spin. They use unweighted least squares in magnitude, at the retained native time stamps, over a narrow ±0.0025 cycles/day primary frequency window at each harmonic order. The primary was refined on the original cleaned measurements, using distance correction and a quadratic time trend. These fits do not forward-model the finite-integration or onboard cosmic-ray-rejection response. [TESS description of prime-mission cadence and integration](https://heasarc.gsfc.nasa.gov/docs/tess/TESS-CosmicRayPrimer.html)

From the parent results directory, after running its measurement script, run the following. The results directory contains the parent measurements CSV.

```sh
PYTHONPATH=. OPENBLAS_NUM_THREADS=1 python alias_checks/detector_geometry.py /path/to/tssys-release . --out alias_checks
PYTHONPATH=. OPENBLAS_NUM_THREADS=1 python alias_checks/spin_alias.py /path/to/tssys-release . --out alias_checks
```

[Exposure and sky checks](detector_geometry.csv) · [Aliased spin models](spin_alias.csv)
