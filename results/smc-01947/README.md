# A 4.252-hour brightness variation toward OGLE-SMC-ECL-01947

A sinusoid with a **4.25168-hour period** remains after subtracting the known 2.40135-day eclipsing binary. OGLE-IV measures an I-band amplitude of **0.01823 mag peak to peak**. Photometry from EROS-2 recovers the predicted phase in two bands. In four TESS visits the strongest residual peak in a search from 0.1 to 15 cycles/day lands at the OGLE frequency. Two other visits rank it second. The evidence establishes a periodic brightness variation along this line of sight. It leaves the identity of the variable star open.

## Resolving a daily alias

The earlier Open-Space search stopped at 5 cycles/day. The strongest entry it retained for this target was 0.215422 days, corresponding to 4.6421 cycles/day. Extending the search to 12 cycles/day finds a stronger peak at **5.64483 cycles/day**. Their difference, 1.00274 cycles/day, is the frequency of Earth's rotation relative to the stars. Observations at similar local times can produce both peaks from a single oscillation. TESS's much more continuous sampling distinguishes them. [Earlier screening table](../smc-01293/data/population.csv.gz)

This sightline was one of five selected for exploratory short-period follow-up (01844, 01947, 02360, 03848 and 05439) after inspecting the earlier search. We fitted 11 orbital harmonics with robust inverse-variance weights to 2,375 OGLE-IV I measurements, searched the clipped residuals over 0.1–12 cycles/day, and refined the strongest peak by jointly fitting the orbital waveform and a sinusoid. All four solar- and sidereal-day alternatives tested give a higher robust loss. The smallest increase is 243. This loss measures model fit and has no assigned population false-alarm probability.

| Photometry | Measurements | Nights | Signed semiamplitude at the OGLE phase (mmag, ±1 SE) |
|---|---:|---:|---:|
| OGLE-IV I, training | 2,375 | 949 | 9.12 ± 0.28 |
| OGLE-IV V | 80 | 79 | 9.00 ± 1.68 |
| OGLE-III I | 1,464 | 720 | 7.32 ± 0.58 |
| OGLE-III V | 53 | 53 | 6.44 ± 1.90 |
| EROS-2 red | 878 | 752 | 11.96 ± 0.83 |
| EROS-2 blue | 959 | 830 | 10.12 ± 0.62 |

Each comparison fits its own orbital waveform and projects the fitted sine and cosine coefficients on the OGLE-IV phase. Errors group measurements by integer HJD and condition on the chosen ephemeris. OGLE-III contributed to the earlier selection. The EROS channels share observing dates. The comparison EROS sources have lower precision. [Individual fits, aliases and controls](result.json)

## Independent TESS search

We retrieved every TESSCut sector returned at (12.668879641°, −73.215620831°). For each visit we used the pixel nearest the binary, subtracted the median of the cutout's four corners, and retained quality-zero cadences. Each visit fit nine harmonics of the known orbit and separate offsets and slopes in 0.8-day intervals. We then ranked local maxima of its residual Lomb–Scargle periodogram. The interval for matching the OGLE peak was ±0.025 cycles/day. The interval and analysis choices are exploratory. Adjacent Fourier samples are correlated.

| TESS sector | 1 | 27 | 28 | 67 | 68 | 94 | 95 |
|---|---:|---:|---:|---:|---:|---:|---:|
| Rank of the matching peak | 25 | 1 | 2 | 1 | 1 | 2 | 1 |

TESS pixels subtend about 21 arcseconds. Exploratory centroid fits in individual sectors were inconsistent, spanning tens of arcseconds, so we leave the variable unlocalized within the crowded field. A short-period pulsation is possible. A different variable inside the photometric aperture also fits the evidence. Resolved high-cadence photometry can distinguish these explanations.

## Do the brightness variations share the proposed third-body delay?

Hong et al. fitted a 6.9-year sinusoid to eclipse arrival times after accounting for apsidal motion. If the 4.252-hour source belongs to the inner binary and remains a stable clock, light-travel time should shift its maxima by the same delay. We fixed their published central phase and angular rate, then fitted the delay amplitude of the new clock. In 1,269 out-of-eclipse OGLE-IV measurements on 516 nights, the signed delay was **−1.0 ± 2.5 minutes**, versus **+11.5 ± 7.2 minutes** from their eclipse-timing fit. Each series was free to choose its constant oscillation period and phase. Testing nine combinations of published phase and rate offsets gave amplitudes from −2.7 to +1.7 minutes in OGLE-IV. OGLE-III and EROS estimates are less precise. These data give no positive evidence that the two clocks share the central published delay. The original amplitude is uncertain, phase and period can drift, and the brightness variation could come from another star. The test leaves the proposed third body open. [Fixed-phase test and sensitivity grid](delay.json)

## Prior coverage and reproduction

This binary is also SMC-SC5-185154, MACHO 212.16020.84 and Gaia DR3 4685932607040745856. Hong et al. reported apsidal motion and a possible third body with a 6.9-year modulation of its **eclipse arrival times**. The four-hour signal here is a direct modulation of **brightness**. We checked the published paper, the OGLE catalog, the VSX cone, and Zhou's catalog of pulsating eclipsing binaries. Its nearest coordinate-readable entry is 73.2 arcminutes away. We found no earlier report assigning the 4.252-hour signal to this sightline. Catalog coverage and alternate identifiers limit statements about priority. [Hong et al. 2016](https://academic.oup.com/mnras/article/460/1/650/2608975) · [Zhou catalog check](prior.json) · [OGLE release](https://www.astrouw.edu.pl/ogle/ogle4/OCVS/smc/ecl/)

Run the following from this directory with Python, NumPy, SciPy and Astropy. The first command uses the frozen numerical data. The last downloads the cutouts and compares newly extracted time and flux samples. Files regenerated by TESSCut can have different headers and compressed bytes.

```sh
OPENBLAS_NUM_THREADS=1 python measure.py
OPENBLAS_NUM_THREADS=1 python delay_test.py
sha256sum -c SHA256SUMS
python extract_tess.py refreshed
```

The frozen OGLE files contain original HJD, magnitude and uncertainty columns. EROS timestamps are HJD minus 2450000. TESS times are BTJD and use an approximate −0.0008-day conversion to HJD(UTC) in the analysis. The small timing approximation is immaterial to the sector periodograms. Credit for the observations belongs to OGLE, EROS-2, TESS and MAST. CDS/VizieR and SIMBAD supplied catalog data. [Input URLs and source digests](sources.json)
