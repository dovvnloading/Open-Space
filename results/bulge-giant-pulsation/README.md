# Three Galactic-bulge giants with changing pulsation periods

**OGLE-BLG-LPV-046289 changed from a predominantly 137–139-day brightness cycle to a 264–271-day cycle.** EROS-2 and early OGLE independently detected the shorter period. Later OGLE and Gaia recover the longer period. In two other Miras, the dominant optical period shortened by 17.5% and 13.0% across successive disjoint intervals.

## OGLE-BLG-LPV-046289

Each row is a fresh, floating-mean sinusoid search over 70–1,200 days, with one measurement per populated day. EROS red and blue share their observing nights. Dropping the overlap between EROS and early OGLE still leaves a 138.0-day cycle in 357 earlier EROS red days.

| Observations | Occupied days | Best period (days) | Variance explained |
|---|---:|---:|---:|
| EROS-2 red, 1997–2002 | 444 | 137.1 | 51% |
| EROS-2 blue, 1997–2002 | 430 | 137.0 | 30% |
| OGLE-III I, 2001–2006 | 310 | 139.1 | 64% |
| OGLE-IV I, 2010–2016 | 602 | 264.7 | 89% |
| OGLE-IV V, 2010–2016 | 24 | 263.8 | 85% |
| Gaia DR3 G, 2014–2017 | 24 | 271.0 | 86% |

The earlier OGLE catalog lists 138.87 and 274.4 days. The later catalog gives 263.3 days. The independently measured change in **which cycle dominates**, and by how much, is the result. Simultaneously fitting the early and late I-band frequencies gives short/long peak-to-peak amplitudes of **0.391/0.211 mag** early and **0.171/1.138 mag** late. The gap between early 2006 and early 2010 leaves the transition date unresolved.

We tested a single repeating waveform with a fundamental and an exact second harmonic, allowing their amplitudes and phases to differ between surveys. The common period that best predicts omitted observing seasons was 273.5 days. Separate optima for EROS red, OGLE-III and OGLE-IV were 275.5, 279.25 and 264.75 days. In late OGLE the RMS error fell from 0.282 to 0.144 mag. Period selection used those same folds. The test disfavors an exactly stationary clock as a description of all three spans. A change in harmonic balance or waveform, with some period evolution, remains viable. Spectroscopy through another change could help distinguish it from a switch of pulsation mode.

Gaia DR3 4068994807280968448 lies 0.11 arcsec from the OGLE position. The next Gaia source is 2.81 arcsec away and 4.74 mag fainter in mean G. The Gaia-corrected EROS position is within 0.11 arcsec of OGLE. The EROS neighbor 3.31 arcsec away has a best-period variance fraction of only 3.4%. The data still allow an unresolved blend.

Two other candidates illustrate the effect of the observation window. At -020182, successive EROS red, OGLE-III I and OGLE-IV I searches give 254.2, 130.7 and 248.8 days. At -151705 they give 296.3, 141.6 and 290.0 days. In both cases the long cycle appears before and after the short-cycle observations. Different passbands and changing harmonic content complicate their interpretation.

## Two rapidly shortening Miras

Within each observing era, we split the time span into equal contiguous intervals, taking the integer part of span/950 days with at least one interval. Each interval received a free frequency search.

| Target | Date range | Successive dominant I-band periods (days) | Independent checks |
|---|---|---|---|
| OGLE-BLG-LPV-225280 | 2001–2016 | **257.7 → 228.6 → 212.6** | Gaia gives 218.4 days over its 2014–2017 window. OGLE-IV V gives 221.7 days across its longer baseline. |
| OGLE-BLG-LPV-178146 | 2001–2019 | **273.5 → 265.8 → 261.0 → 242.7 → 237.9** | Three-day I-band bins give 274.3 → 267.5 → 261.2 → 244.0 → 238.3 days. OGLE-IV V gives 243.3 days. Gaia's strongest daily-binned peak is near the half-period, 118.5 days. |

At -225280 the three fits explain 93%, 93% and 95% of the variance. Dropping one season at a time gives period ranges of 255.7–259.6, 226.7–229.3 and 212.0–213.5 days. At -178146 the first and last ranges are 272.2–275.2 and 235.5–242.2 days. Its last fit explains 60% of the variance, down from 94% initially. The late signal has become less regular.

We selected up to five angularly nearest earlier-catalog Miras in each target's OGLE-IV field, within 1.2 mag and 0.18 in absolute natural log period. None of the nine comparison stars changes by more than 4% between its first and last block. Selection was exploratory and these controls measure local behavior. Later optical observations can distinguish persistence below 220 and 245 days, respectively, from a return toward the early periods. A helium-shell flash is one possible physical explanation among several. [Prior work on changing-period Miras](https://doi.org/10.1051/0004-6361/201116463)

## Sample, data and priority

Exact OGLE identifiers link 6,539 earlier Mira or semiregular stars to the later bulge Mira catalog. **52** have an absolute log period ratio above 0.07. Both era-specific public OGLE-IV-release I-band files existed for **45** of them. Detailed selection used temporal coverage, the size of the discrepancy and its behavior within an era. The largest catalog discrepancy, 319.8 versus 104.8 days for -001386, has only 374 days of later coverage and a strong competing period. The downloadable catalogs have SHA-256 checks in the replay script. The full screen records each missing file. [Soszyński et al. OGLE-III](https://ogle.astrouw.edu.pl/cont/4_main/var/ogleiii/lpv_blg/) · [Iwanek et al. OGLE-IV](https://arxiv.org/abs/2203.16552)

All these stars, and their catalog periods, were previously published. Iwanek et al. explicitly discuss Mira and semiregular classifications changing with amplitude. In the SIMBAD object bibliographies checked on 14 November 2026, each of the five stars has entries for the OGLE-III, OGLE-IV, OGLE Galactic map and Gaia DR3 long-period-variable papers. We searched the OGLE designations, Gaia IDs, cg1183l11257, 2MASS J18120215-2812421 and J18011508-2649316 together with terms for period change, mode switching and thermal pulse. No earlier object-specific report of the independently checked amplitude reversal or the two within-era shortening histories emerged. Historical priority remains open.

NEOWISE provides only 14–21 half-year visits per target after flag and position cuts. Multiple aliases appear, including approximately 93 and 176 days in W1 at the two shortening targets, and we excluded the infrared estimates from the shortening claim. The nearest VIRAC2 main-catalog sources for the three headline stars are 2.25–2.99 arcsec from the optical positions. ZTF returned zero detections within the tested cones.

## Reproduction

Run the code from this directory with Python, NumPy, pandas and SciPy 1.15 or later. The first invocation fetches the target, comparison-star and external photometry. Subsequent invocations reuse the saved bytes. The optional catalog screen adds the original three catalog files and all 52 follow-ups.

```sh
python replay.py /tmp/bulge-replay --screen
```

The command writes the target measurements, a complete screen CSV and a source receipt with URLs, sizes and hashes. The published block-level measurements are available in [observations.csv](observations.csv). The light curves are the OGLE and EROS heliocentric times relative to JD 2450000, and Gaia's barycentric TCB times relative to JD 2455197.5. Only periods are compared across missions. OGLE photometric errors must be positive and below 0.35 mag. EROS requires 8–23 mag and positive errors below 0.18 mag. Gaia measurements rejected by photometry or by the G-band variability pipeline are removed. Digit precision in machine-readable output denotes numerical precision.

Data credit goes to the OGLE and EROS collaborations, ESA/Gaia and DPAC. VizieR is operated by CDS. [EROS-2 archive release](https://doi.org/10.1051/0004-6361/202557321) · [Gaia epoch photometry](https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_photometry/ssec_dm_epoch_photometry.html)
