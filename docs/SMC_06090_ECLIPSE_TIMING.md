# The changing eclipse interval of OGLE-SMC-ECL-06090

The time from the primary eclipse to the following secondary eclipse in this 28.07-hour binary first shortened, then lengthened. OGLE measured an increase of **2.21 hours** between its 2010–2014 and 2022–early 2024 observations. TESS and Swift independently detected a later increase in optical and ultraviolet light. An eccentric orbit that rotates slowly predicted the two latest TESS eclipse intervals within three minutes when fitted only to pre-2022 observations. A single dip in the interval followed by a plateau also fits the available record. Whether the interval will begin shortening again is the open question. [Evidence and reproduction instructions](../results/smc-06090/README.md) · [Illustrated report (PDF)](../results/smc-06090/orbit-test.pdf)

## Which star and which interval?

OGLE-SMC-ECL-06090 is a cataloged eclipsing binary at J2000 **01ʰ 24ᵐ 15.06ˢ, −73° 19′ 33.4″**. Its earlier name is OGLE-SMC-ECL-6090. Other identifiers are SMC121.3.2046, SMC738.16.542 and Gaia DR3 4686449235706841984. Throughout this analysis, *primary* means the dip near phase zero of the fixed OGLE-IV ephemeris, even when the two depths are similar. The reference period is 1.1696719 days and the primary epoch is HJD 2460000.0369. The observable is the time **within one orbit** between the fitted centers of the two dips. It can change while the duration of a full orbit stays nearly constant. [OGLE-III target photometry](https://www.astrouw.edu.pl/ogle/ogle3/OIII-CVS/smc/ecl/phot/I/OGLE-SMC-ECL-6090.dat) · [OGLE-IV target photometry](https://www.astrouw.edu.pl/ogle/ogle4/OCVS/smc/ecl/phot/I/OGLE-SMC-ECL-06090.dat)

The original Open-Space screen had 5,213 adequately sampled binaries. It ranked changes of the same sign in two fixed eclipse windows. Of 2,275 systems with opposite signs, five had apparent changes of at least 0.10 mag in each window. This target ranked highest among those five by the weaker of its two nominal signal-to-error ratios. The target selection was exploratory. Fitting the centers showed that the eclipses themselves moved relative to the fixed windows. [Original survey and selection](SMC_ECLIPSE_DEPTHS.md)

## Optical timing and source identification

The 637 OGLE-III and 711 OGLE-IV I-band measurements span June 2001 to February 2024. Each point passed the finite-value and positive-error cut of at most 0.10 mag. Two phase-wrapped Gaussian profiles were fitted per epoch, with free centers, widths, depths and baseline, and a 0.012-mag error floor. Resampling observing nights 800 times gave the intervals below. Each row estimates centers from many folded orbits.

| OGLE epoch | Measurements | Primary-to-secondary interval (orbital cycles) | 95% interval |
| --- | ---: | ---: | --- |
| 2001–early 2005 | 319 | 0.4982 | 0.4936–0.5032 |
| Mid-2005–2009 | 318 | 0.4657 | 0.4621–0.4694 |
| 2010–2014 | 382 | 0.4610 | 0.4578–0.4644 |
| Late 2014–mid-2018 | 147 | 0.4938 | 0.4899–0.4991 |
| Mid-2018–early 2020 | 86 | 0.5233 | 0.5126–0.5298 |
| 2022–early 2024 | 96 | 0.5398 | 0.5348–0.5474 |

During OGLE-III the interval shortened by 0.0325 cycles (95% interval for the change −0.0391 to −0.0266). Across OGLE-IV the primary moved earlier by 0.0411 cycles and the secondary later by 0.0377. Their separation grew by 0.0788 cycles, or 2.21 hours (95% interval 2.04–2.43 hours). Triangular eclipse profiles and refits using the independent Gaia catalog period recovered the change in the same direction.

TESS covered the target in sectors 1, 2, 13, 27, 28, 67, 68, 94, 95 and 96. Stacking the covered orbits produced intervals of **14.43 hours** in 2018–2019, **14.73** in 2020, **15.21** in 2023 and **15.14** in 2025. The earliest-to-latest change was 42 minutes (95% orbit-block bootstrap interval 24–63 minutes). A second background subtraction gave 29 minutes of increase and the same late value to two decimal places. A 2023-to-2025 change of −3.8 minutes had a 95% interval of −23.4 to +17.5 minutes, leaving the direction of the most recent change unresolved. [NASA TESS observations through MAST](https://mast.stsci.edu/tesscut/)

TESS and Swift blend a star 1.50 arcseconds from the target. DELVE-MC resolves the neighbor and measures both eclipses at the target position. Its red-band primary and secondary dips were 0.425 and 0.393 mag on two separate nights. Even complete loss of the neighbor's catalog flux would reduce their combined red light by only 0.087 mag. Gaia independently recorded primary and secondary dimming in its 2014–2017 G-band observations. [NOIRLab DELVE data](https://datalab.noirlab.edu/data/delve) · [ESA Gaia archive](https://gea.esac.esa.int/archive/)

## Ultraviolet observations

Swift's UVW1 images provided 275 accepted measurements on 267 nights between September 2012 and October 2026. At the OGLE-predicted phases, the target-to-reference flux ratio fell **41.2% ± 3.1%** during primary eclipse and **31.8% ± 5.4%** during secondary eclipse. These are differential aperture-flux depths with one-standard-error uncertainties grouped by night. The nearby comparison star showed no comparable dip at either phase. [Swift archive](https://heasarc.gsfc.nasa.gov/W3Browse/swift/swiftmastr.html)

With each eclipse center and depth fitted freely in each of two date ranges, Swift gave **14.66 hours in 2016–2020** (124 exposures) and **15.49 hours in 2021–2026** (142 exposures). The increase was **50.3 minutes**, with a 95% interval of **1.6–115.9 minutes** from 999 year-stratified resamplings of observing nights. The two groups exclude the isolated 2012 visit and use the optical folding period. Three aperture sizes gave increases of 41–50 minutes. Sparse sampling of the late secondary eclipse limits its precision. All three accepted exposures in 2026 missed both eclipse cores.

XMM-Newton's Optical Monitor caught a primary eclipse in observation **0741450101** on March 24, 2015. Three UVM2 exposures measured Vega magnitudes of 16.442, 16.866 and 16.435. The middle exposure was **0.428 ± 0.056 mag** fainter than the average of its neighbors. None of 75 similar-brightness field controls dimmed by as much. The predicted primary minimum fell within that 66.7-minute exposure, which establishes dimming but leaves the exact center uncertain. [XMM-Newton Optical Monitor archive](https://heasarc.gsfc.nasa.gov/W3Browse/xmm-newton/xmmomsuss.html)

## Does the eccentric orbit precess?

If the long axis of an eccentric orbit gradually turns, the stars pass through conjunction at changing orbital speeds and phases. To test the timing consequence, seven models were fitted to the same pre-2022 record of 1,252 OGLE I measurements, 43 Gaia G transits, and the 2018–2019 and 2020 TESS groups. All later observations were withheld. The two subsequent TESS estimates used 38 and 51 covered orbital cycles. The table shows three alternatives to the uniformly precessing Kepler orbit. Lower loss indicates more accurate predictions for the last 96 OGLE points. The test was retrospective.

| Pre-2022 model | Predicted TESS interval, 2023 | Predicted TESS interval, 2025 | Absolute errors, 2023 / 2025 | Held-out OGLE robust loss |
| --- | ---: | ---: | ---: | ---: |
| Precessing eccentric orbit | 15.246 h | 15.192 h | 2.4 / 3.0 min | 133.44 |
| Single trough returning to a plateau | 15.255 h | 15.406 h | 2.9 / 15.8 min | 135.91 |
| Sinusoid with a free mean | 14.956 h | 14.737 h | 15.0 / 24.3 min | 140.34 |
| Cubic trend | 15.661 h | 16.076 h | 27.3 / 56.0 min | 202.85 |
| **TESS measurement** | **15.206 h** | **15.142 h** | | |

With the entire OGLE, Gaia and TESS timing record included, the edge-on Kepler model gave **28.25 years** for one rotation (conditional 95% interval 27.59–28.87 years) and eccentricity **0.0672**. Directly calculating the overlap of two limb-darkened stellar disks in the OGLE I exposures gave **28.45 years** and **0.0669** for the forward-precessing fit. Both calculations use the same OGLE photons. Assumed extra light changes the inferred inclination from 83.9° to nearly 90°. Outside eclipse, the brightness modulation at twice the orbital frequency is weaker than a simple equilibrium-tide estimate predicts.

The optical record spans less than one fitted cycle. A single nonrepeating trough in the eclipse spacing, followed by a return to a plateau, came close to the precessing model on both the training and withheld data. In **102 of 300 simulations generated from the fitted trough**, the precessing model gained at least as much held-out predictive advantage as in the observations. The empirical finding is changing separation with an observed reversal. Periodicity, the sense of precession and its physical driver remain unresolved.

## How unusual is the reversal?

A comparison selected 676 non-contact SMC binaries of similar catalog magnitude, period and eclipse depth. Of 442 with adequate data in both early and late windows in **both** OGLE campaigns, the target and four others had opposite-direction changes at least as large as the target's changes. The comparison recovered OGLE-SMC-ECL-03825, already flagged for visible apsidal motion. This is a count of measured behavior among eligible systems. Selection of the original target after inspecting the survey prevents treating it as a discovery probability.

A November 13, 2026 search of target identifiers, OGLE remarks, SIMBAD references and published SMC apsidal-motion target lists found no earlier report of this particular reversal. Hong et al.'s study of 90 systems required I brighter than 18.0 mag. The target has OGLE-III I = 18.232 mag. One of those 90 systems had a shorter orbital period. Other aliases, surveys and unindexed work limit any priority claim. [OGLE-III catalog](https://arxiv.org/abs/1310.3272) · [OGLE-IV catalog](https://arxiv.org/abs/2503.15596) · [Hong et al. 2016](https://doi.org/10.1093/mnras/stw955) · [Zasche & Wolf 2019](https://arxiv.org/abs/1902.04961)

## An observation that can separate the explanations

On **July 1, 2027** the full-data precession and plateau fits predict gaps of **14.92 and 15.32 hours**, respectively. By July 2028 the difference grows to 37 minutes. A second observable is the eclipse width. For July 2027, the forward-precession models predict the secondary lasts roughly 9–12 minutes longer than the primary at half depth. Reversing the precession gives a secondary about 8–9 minutes shorter. Present TESS width estimates depend too strongly on the aperture and noise model to settle this test.

One of 17 feasible two-site pairs is **October 17–18, 2027**. Siding Spring, Australia, can cover the primary on October 17 from **11:06 to 15:25 UTC**. Cerro Tololo, Chile, can cover the following secondary on October 18 from **02:07 to 06:26 UTC**. For that orbit the two models predict **891.7 and 919.5 minutes**, a 27.8-minute difference. Repeated, spatially resolved measurements covering ingress, egress and baseline in the same filter would test both the spacing and width ordering. The windows describe modeled target and solar altitude. Telescope time has not been allocated. [Observing windows, fitted results and replay instructions](../results/smc-06090/README.md)

The original measurements belong to the OGLE collaboration, NASA TESS and STScI/MAST, ESA Gaia/DPAC, the Neil Gehrels Swift Observatory and S-CUBED investigators, ESA XMM-Newton and its Optical Monitor team, NASA HEASARC, and the DELVE collaboration and NOIRLab. This repository records an archival reanalysis, with status as of November 14, 2026.
