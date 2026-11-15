# Pulsed brightenings on two “missing” orbits of EXO 2030+375

EXO 2030+375 brightened on **2022 September 28** and **2023 September 26**, during two of the eight orbits marked as missing outbursts in Figure 2 of Huang, Ji & Yan (2025). MAXI/GSC resolved both brightenings in the 2–20 keV band. Fermi/GBM recorded a simultaneous rise and fall in the pulsed flux of the 41-second neutron star. Both dates pass comparisons with Cyg X-3 and Cyg X-1. On the second date an independent fold of MAXI photons recovers the spin at the frequency supplied by GBM.

Earlier work established that bright outbursts occurred about every two 46-day orbits after the 2021 giant outburst. The narrower result here is that some intervening orbits also hosted weaker, resolved activity at the pulsar. Huang et al. defined “missing” by the significance of their fitted outburst component. Our local flux contrast and frequency-selected tests measure different quantities. They do not establish the amplitude their fitting procedure would have recovered or the mechanism behind the modulation. [Huang et al.](https://doi.org/10.3847/1538-4357/adc573) · [Published figure](https://content.cld.iop.org/journals/0004-637X/984/1/66/revision1/apjadc573f2_lr.jpg) · [Nakajima et al. 2023](https://www.asj.or.jp/nenkai/archive/2023a/pdf/W39a.pdf)

| Measurement | 2022 September 28 | 2023 September 26 |
|---|---:|---:|
| Figure 2 arrow, counted from left | 3 of 8 | 7 of 8 |
| GBM selection epoch, MJD | 59850.52 | 60213.01 |
| Offset from nominal periastron | +7.07 days | +1.39 days |
| GBM pulse period | 41.1267 s | 41.0913 s |
| MAXI 2–20 keV peak minus local baseline, photons cm⁻² s⁻¹ | 0.0555 ± 0.0059 | 0.0695 ± 0.0057 |
| MAXI distinct days, left / central / right windows | 5 / 5 / 5 | 5 / 5 / 5 |
| Central GBM flux above the preceding / following window, in formal errors | 8.7 / 7.1 | 4.7 / 3.6 |
| MAXI 3–15 keV photons / scans for independent timing | 4,531 / 93 | 4,373 / 93 |
| Two-harmonic scan-rotation tail fraction | 0.0138 | 0.00030 |

## What the three instruments measured

We started with six Swift/BAT 15–50 keV daily histories. We retained quality-zero rows with positive errors. The exploratory screen used five-day inverse-variance means, peak spacing of at least 0.24 of each approximate orbital period, prominence above 0.003, rate above 0.005 count cm⁻² s⁻¹ and nominal signal/error above seven. The 137 EXO peaks include secondary peaks within long episodes. The per-source denominators and dispositions are in [review.csv](review.csv). The six-source list was committed before the photometry was inspected.

The eight-arrow sample came from the published figure after the exploratory screen. For the GBM ephemeris, periastron is MJD 52756.17 + 46.0213 × orbit. The arrow positions are offset slightly from exact periastron in the drawing. We use the corresponding orbit numbers. On each orbit we select the brightest detected GBM bin centered between −1 and +12 days relative to periastron. Three duplicate time stamps in the full GBM history are collapsed by retaining the first occurrence. GBM flux is RMS pulsed energy flux in keV cm⁻² s⁻¹. Its selected epoch is a multi-day measurement, so the rounded MJD does not imply a precisely known maximum. [GBM source, ephemeris and data](https://gammaray.msfc.nasa.gov/gbm/science/pulsars/lightcurves/exo2030.html)

For MAXI and BAT we compare a ±2.5-day window centered on the GBM time against equal-width windows centered six days before and after. The reference line interpolates between the inverse-variance weighted times and fluxes on both sides. The uncertainty propagates the three bins' quoted errors, assuming independent measurements. Missing coverage remains blank. With alternate (half-width, spacing) pairs of (1.5,4), (1.5,8), (2,7), (2.5,7.5) and (3.5,8) days, the two MAXI excesses range from 5.5 to 10.5 and 8.0 to 15.7 formal errors. At the nominal windows the 4–10 keV excesses alone are 10.0 and 11.1 formal errors. The per-event GBM rise/fall check compares the central ±1.8 days with measurements 1.8–8.5 days before and after it. GBM selected the centers, so its peak selection favors local maxima. The independent MAXI and BAT flux measurements use those centers without optimizing their own peak times. The analysis choices were exploratory. The formal errors supply neither survey-wide significance nor noise estimates for all intrinsic variability. [MAXI source light curves](https://maxi.riken.jp/star_data/J2032+376/J2032+376.html) · [Swift/BAT source](https://swift.gsfc.nasa.gov/results/transients/EXO2030p375/)

Cyg X-3 is 3.32° from EXO. Its simultaneous MAXI 2–20 keV contrast is +0.2 and +0.5 formal errors on the two dates. The corresponding values for Cyg X-1 are both −2.0. Both comparison sources do brighten on 2023 June 26 (+6.5 and +41.6 formal errors). We therefore withhold a clean-source interpretation of the June MAXI amplitude despite GBM activity at the EXO frequency. A comparison-source coincidence by itself does not quantify leakage.

## Independent timing and the other six passages

We selected 3–15 keV MAXI/GSC low-bit-rate events within 1.25° of EXO in the six days centered on each GBM epoch rounded to the nearest half MJD for event extraction. The flux contrasts use the unrounded epochs. We applied the position, energy and detector quality cuts in [replay.py](replay.py). Sky tiles 139 and 171 cover the EXO aperture. The same tiles and EXO timing corrections define a nearby control aperture centered on Cyg X-3. That control is a check on field contamination at the EXO frequency. We used the archive TT reference, the Earth light-time correction and the GBM eccentric binary orbit including apsidal advance. Omitting the ISS geocentric position changes arrival times by at most about 0.023 s. We used a constant spin frequency within each six-day test. Both timing instruments share the orbital model. MAXI's daily flux and its photons also share an instrument.

The test evaluates the fundamental and second harmonic at the one GBM frequency. No frequency grid is searched. The 2023 September event has powers 114.1 and 5.0 in the EXO aperture versus 1.6 and 3.6 in the nearby control. For a sampling-aware null, we define scans by gaps over 400 s, rotate each scan's entire photon phase pattern at random, and recompute the sum of harmonic powers. Two of 10,000 draws reach or exceed the observed power for September 2023. Including the observed case gives 3/10,001. A Bonferroni bound across eight flagged orbits is 0.0024 for that fixed two-harmonic diagnostic. It does not cover the selection of this source, energy cut, aperture, duration or statistic. Nonuniform sampling within each scan is preserved. Recurring phase-locked instrumental effects across scans are a limitation of this null. [MAXI archive and sky tiling](https://heasarc.gsfc.nasa.gov/W3Browse/maxi/maximaster.html)

| Arrow | Approximate peak MJD | Disposition |
|---:|---:|---|
| 1 | 59658.50 | A declining bright tail dominates. No local rise in the pulsed flux. |
| 2 | 59757.52 | Weak local structure. The MAXI 2–20 keV contrast is 1.4 formal errors. |
| 3 | 59850.52 | Resolved MAXI brightening and GBM rise/fall. The BAT local excess is 2.9 formal errors. |
| 4 | 59943.00 | A pulsed rise/fall and a 3.9-formal-error BAT excess. No central MAXI flux coverage. |
| 5 | 60029.01 | A 2.8-formal-error BAT excess with partial coverage. No central MAXI flux coverage. |
| 6 | 60121.01 | A strong MAXI excess coincides with activity in both comparison sources. |
| 7 | 60213.01 | Resolved MAXI brightening, GBM rise/fall, and independent MAXI pulse recovery. BAT lacks a complete sideband. |
| 8 | 60401.04 | Sparse MAXI data and no compelling pulsed rise/fall. |

## Replay and prior reports

[observations.csv](observations.csv) contains the rounded measurements for all eight passages. [flagged-orbits.csv](flagged-orbits.csv) fixes the Figure 2 selection. Run `python replay.py ./cache` from this directory with Python 3, NumPy, SciPy and Astropy, or add `--monitors-only` for the flux analysis. Event files are available from HEASARC. The original screen examined six Swift snapshots generated on 2026 October 21; the source data were retrieved on 2026 November 15. BAT's integer daily dates are assigned to bin centers by adding 0.5 d. MAXI uses the published spacecraft-frame daily centers and GBM uses BARYTIME. The clock differences are negligible at the daily flux resolution. The GBM file and the daily monitor endpoints roll forward. The replay writes per-file SHA-256 values including explicit 404 gaps and unrounded outputs for comparison. A changed live file can change the answer. The archives, rather than synthetic or generated photons, supply every measurement.

Searches on 2026 November 15 used the source name, each of the four principal MJDs and dates, and combinations of “missing,” “weak outburst,” “alternating,” “pulsed flux,” and “periastron” in journal, arXiv, ADS-indexed and conference records. Nakajima et al. (2023) reported the alternating luminosities. Huang et al. (2025) analyzed the total-flux histories and marked the eight passages. Wilson et al. recovered faint outbursts in the *1990s*. Thalhammer et al. describe subsequent NICER monitoring but restrict their published spectral/pulse analysis to the earlier giant outburst. We located no report measuring the two quoted 2022–23 brightenings. A particularly useful falsification would be a joint scan-by-scan spatial response fit for EXO and the Cygnus sources, followed by a blind pulse search in the separated photons. That could test source leakage and independence from the GBM frequency selection. [Huang et al.](https://arxiv.org/abs/2505.05102) · [Nakajima et al.](https://www.asj.or.jp/nenkai/archive/2023a/pdf/W39a.pdf) · [Wilson et al.](https://ntrs.nasa.gov/citations/20020022487) · [Thalhammer et al.](https://doi.org/10.1051/0004-6361/202348594)
