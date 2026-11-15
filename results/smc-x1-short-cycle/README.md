# Four short X-ray brightness cycles in SMC X-1

Swift recorded five consecutive broad minima between January 25 and June 25, 2025. The intervals were **38, 36, 37 and 40 days**. Four cycles took 151 days, about 69 days fewer than four at the preceding 55-day rate. An independent fit to MAXI measured a 38.6-day cycle. The modulation subsequently lengthened to about 49 days in both instruments and to 58.5 days in MAXI in 2026. These are measurements of the recurring bright and faint states of this known neutron-star binary.

| Observing window (MJD, upper bound exclusive) | Swift/BAT, 15–50 keV | MAXI/GSC, 2–20 keV | Days observed, Swift / MAXI |
| --- | ---: | ---: | ---: |
| 60490–60695, preceding interval | 55.0 d | 54.8 d | 170 / 178 |
| 60695–60860, short cycles | 38.9 d | 38.6 d | 152 / 116 |
| 60865–61080, subsequent interval | 49.4 d | 49.0 d | 173 / 157 |
| 61136–61349, April–November 2026 | Insufficient coverage | 58.5 d | 34 / 187 |

![Swift and MAXI light curves over the four short cycles and the fitted cycle lengths in the four observing windows](cycle-timing.png)

The upper panel shows four-day weighted averages, with each instrument centered on its median and scaled by its own 10th-to-90th percentile range. Lines show separate sinusoidal fits. Dotted lines mark the five Swift minima. The lower panel shows the duration of the windows used for the fits. The horizontal spans indicate dates, and the values describe a changing clock.

## What was known

Mihara reported shortening to about 45 days in September 2023 and an early turn-on in May 2024. Hu et al. measured a 47-day cycle ending January 29, 2025 (MJD 60704.5), and displayed monitoring through MJD 60730 in their spin analysis. The beginning of our sequence overlaps that paper. The four-cycle measurement, its historical rank under the stated cuts, and the subsequent recovery are the additions. A much earlier study identified four excursions through 2022. These dates and studies leave open whether the 2025 episode merits a separate numbered excursion. [MAXI notice](https://maxi.riken.jp/news/en/MAXI_News_s5.html) · [NinjaSat campaign, Sections 3.1 and 3.4](https://arxiv.org/html/2511.05016v1) · [MOOSE monitoring through MJD 59810](https://doi.org/10.1093/mnras/stad384)

Niwano et al. studied optical and X-ray modulation through TESS Sector 68. A separate 2026 spectroscopy study used observations through 2022. Searches on November 14, 2026 covered web-indexed arXiv, journal articles, ATels, and MAXI notices, including “SMC X-1 excursion 2025 2026”, “SMC X-1 38 days”, “SMC X-1 60700 60800”, and Japanese terms for superorbital shortening. We found the earlier observations above. Priority remains open to unindexed material. [Optical study](https://doi.org/10.1093/pasj/psag069) · [Spectroscopy](https://doi.org/10.1093/mnras/stag168)

## Cycle timing and competing explanations

We inspected monitoring products for SMC X-1, Her X-1, LMC X-4, EXO 2030+375 and 4U 1820–30 while choosing a question. SMC X-1 and the reported windows were selected after inspection. The fitted periods use native daily averages, inverse-variance weights, a free offset and a single sinusoid. We searched 24–76 days in 0.25-day steps and refined the best grid interval. The independent minimum count uses a three-day Gaussian smoothing width, inverse-variance normalization across missing days, minimum separation of 29 days and minimum prominence of 0.0023 count cm⁻² s⁻¹ for BAT. The cycle count is particularly useful because the shape changes and the individual periods differ.

- **Historical comparison.** Across 6,899 eligible Swift days beginning in 2005, we ranked four consecutive intervals between five detected minima. Each minimum needed a real measurement within five days. Each run needed consecutive minima no more than 75 days apart and no observation gap longer than 12 days. The 2025 sequence ranked first of 112 eligible overlapping runs. The next nonoverlapping sequence averaged 41.75 days in 2014–2015. Two- and four-day smoothing put the 2025 sequence first as well. Five- and six-day kernels merge a shallow minimum and miss the complete run. This rank depends on detectability, the filters and the sampling. The full list of eligible runs is in the results.
- **Sampling and aliases.** Swift observed all five trough dates with a maximum three-day gap during the short window. MAXI independently fits 38.6 days but has gaps up to 17 days, and three smoothed trough dates fall 5–7 days from an actual sample. Omitting one of four time blocks pushes one MAXI fit to 63 days. The corresponding Swift tests span 38.2–40.2 days. Shifting both edges of the short window by up to ten days gives 38.5–38.9 days in Swift and 38.1–38.6 in MAXI. A two-harmonic fit to eclipse-excised Swift pointings can prefer a doubled period at the upper grid boundary. Daily photometry and individual troughs resolve the repeated bright/faint cycle.
- **The 3.89-day binary eclipse.** We used the NinjaSat orbital ephemeris to remove pointings within 0.17 orbital cycles of mid-eclipse, then averaged the remaining measurements by day. A fundamental search on these data gives 38.5 days in Swift and 38.7 in MAXI. Holding the Swift daily waveform's period and phase fixed and fitting only a scale and offset lowers MAXI's out-of-eclipse weighted RMS from 0.05411 to 0.03639 photons cm⁻² s⁻¹. In the eclipse core (within 0.055 cycles) the RMS is 0.02525 before and 0.02522 after. This phase test ties the modulation to the eclipsing source and disfavors a varying neighbor or common background as the dominant signal. The two reductions of each mission share photons.
- **A third telescope gave an inconclusive result.** Fermi/GBM Earth-occultation data include 26 four-day bins in the short window. Their best free period is 26.1 days. The Swift waveform reduces the 12–25 keV weighted RMS only from 45.98 to 45.73 mCrab. We count Swift and MAXI as the two supporting instruments. Fermi GBM Earth occultation results provided by the Fermi GBM Earth occultation Guest Investigation teams at NASA/MSFC and LSU. The team labels the product preliminary and flags 50–55-day sampling artifacts. [Fermi/GBM monitoring notes](https://gammaray.nsstc.nasa.gov/gbm/science/earth_occ.html)
- **A second sky position.** Applying the same window and fundamental search to LMC X-4 yields 30.0 days in Swift and 30.4 days in MAXI during the short interval. The publicly measured cycle of this independent binary is about 30 days. [LMC X-4 timing study](https://doi.org/10.3847/1538-4357/ae4691)

The two missions measure the duration of the apparent brightness cycle. A changing disk shape can move the trough, even if the underlying precession rate changes less. A direct measurement of disk orientation or phase-resolved polarimetry could distinguish these explanations. We assign no population false-alarm probability or exact physical precession period.

## Energy dependence during the four short cycles

MAXI's Gas Slit Camera measured a brighter middle third than the two outer sixths combined in each of the four Swift-defined cycles, in all three energy bands. The difference is smallest at 10–20 keV. All 1,509 MAXI orbit bins between the five Swift minima (MJD 60700.5 to 60851.5, upper bound exclusive) are released with their original subband fluxes and errors. One bin had a zero subband error. The other 1,508 pass the existing broad-band quality cut and have positive errors in every band. [Orbit-bin subset](data/maxi-orbit-bands.csv.gz) · [Numerical results by cycle](energy-bands.json) · [MAXI column definitions and units](https://maxi.riken.jp/top/lc_readme.txt)

| MAXI band | Middle third | Outer sixths combined | Middle minus outer | 95% within-cycle day bootstrap |
| --- | ---: | ---: | ---: | ---: |
| 2–4 keV | 0.04638 | 0.01729 | +0.02909 | +0.02192 to +0.03602 |
| 4–10 keV | 0.05707 | 0.01825 | +0.03883 | +0.03289 to +0.04486 |
| 10–20 keV | 0.02625 | 0.01880 | +0.00745 | +0.00052 to +0.01419 |

Values are band-integrated photon fluxes in photons cm⁻² s⁻¹. We interpolated phase separately between consecutive *Swift* troughs, removed MAXI pointings within 0.17 orbital cycles of binary mid-eclipse, and included negative fluxes. We averaged orbit bins within each sampled integer MJD day, days within each phase region, and then the four cycles with equal weights. The middle and outer intervals contribute 434 and 247 orbit bins, respectively, across 44 and 27 day-by-cycle cells. Flux cuts could favor positive noise in the faint sample. We applied quality cuts only to flux uncertainties and finite values.

The 10–20 keV fraction of the summed three-band flux rose from 0.202 in the middle to 0.346 near the troughs. The middle-minus-outer difference is −0.144, with a stratified day-bootstrap interval of −0.204 to −0.085. Using windows a quarter-cycle wide gives −0.135. Omitting each cycle in turn gives −0.178 to −0.125 with the original windows. The bootstrap resamples days within the eight cycle and phase cells and preserves all three bands together. It conditions on these four cycles and the Swift trough times. Serial dependence, diffuse background, the instrument response and missing days can affect the estimate. In particular, the eclipse-core control also has a positive middle-minus-outer contrast of 0.0109 and 0.0104 in the two softer bands. Its 10–20 keV contrast is −0.0039. The core has only 13 middle and 10 outer day-by-cycle cells. Its near-zero outer summed flux makes a hardness ratio unstable.

Subtracting the mean eclipse-core flux in each band and phase interval reduces the middle-minus-outer hard-fraction difference to **−0.031** (95% day-bootstrap interval −0.166 to +0.121). All 10,000 resamples had positive summed fluxes in both states. The middle interval remains brighter in each band, by 0.0182, 0.0284 and 0.0114 photons cm⁻² s⁻¹. This control can remove real scattered emission from the binary as well as field background, and the sparse eclipse sample widens the interval. The aperture-level color change therefore does not establish an intrinsic change in spectral shape.

Hu et al. found stable hardness in an earlier NinjaSat cycle using different bands and phase selection. Here the measured change concerns the relative band fluxes in the MAXI aperture during the four short cycles. Source hardening, a hard persistent contribution at low flux and energy-dependent obscuration remain possible. An absorption column or emission-line strength would require a response-folded spectrum. [NinjaSat timing and spectroscopy](https://arxiv.org/abs/2511.05016)

### Pointed Swift/XRT photon energies

We retrieved the cleaned photon-counting event file and exposure map for each of 15 Swift survey pointings in the original short-cycle window. They provide 666.9 s of XRT exposure near the optical position. Four exposures fall in the middle third outside binary eclipse. Two fall in the outer sixths outside eclipse. In the 6–20-pixel annulus around SMC X-1, the middle observations give 349 events in 167.35 s of live time and the outer observations give 14 in 82.43 s. A 30–60-pixel local annulus contains 115 and 3 events, respectively. The numbers below use exposure-map area for background subtraction and mean annular exposure for the rates. [Observation IDs, event files, exposure maps and SHA-256 digests](data/xrt-pointings.csv) · [Compact photon list](data/xrt-events.csv) · [Counts for all 15 visits and both annuli](xrt-counts.json)

| XRT pulse-invariant energy | Middle events | Outer events | Middle net rate (counts/s) | Outer net rate (counts/s) |
| --- | ---: | ---: | ---: | ---: |
| 0.3–1 keV | 47 | 1 | 0.286 | 0.010 |
| 1–2 keV | 98 | 4 | 0.601 | 0.053 |
| 2–10 keV | 204 | 9 | 1.242 | 0.117 |

The observations sample a small subset of each phase interval. XRT independently localizes the bright signal and shows a larger rate in all three bands during the sampled middle visits. Its faint-state counts are too few to establish a hardness change. Several full-aperture rates cross the photon-counting pile-up threshold. We discard the central 6 pixels and repeat at 10 pixels, which still gives a larger middle rate in every band. The annular rates omit the PSF-core and wing corrections, and the local annulus can include source wings. Exposure-area weighting allows for bad pixels and unequal background coverage under a uniform-surface-brightness approximation. These measurements describe detected counts. Physical continuum or line fits would require per-visit response files, an encircled-energy correction, stronger pile-up checks and more exposure. [UK Swift pile-up guidance](https://www.swift.ac.uk/analysis/xrt/pileup.php) · [UK Swift spectral analysis](https://www.swift.ac.uk/analysis/xrt/spectra.php)

We selected Swift pointings with start MJD in [60695, 60860] within 0.25° of (19.2714°, −73.4433°) in the HEASARC swiftmastr table. All 15 are centered on the survey tile SMCJ0117.2−7326. For extraction we projected the SMC X-1 optical coordinate (19.27144054°, −73.44333744°) through each EVENTS sky WCS. The compact event table contains every cleaned event within 60 sky pixels, before grade and energy cuts. We require photon-counting grades 0–12 and PI channels 30–1000 inclusive. A channel width is 0.01 keV. Times are exposure-weighted GTI midpoints in the mission MJD system. The core excision, wings, background and image aperture use sky-pixel radii with an outer-exclusive boundary. The XRT pointings and MAXI bins are separate observations. Swift/BAT supplies only the superorbital phase boundaries for these comparisons. [HEASARC Swift master catalog](https://heasarc.gsfc.nasa.gov/W3Browse/swift/swiftmastr.html) · [SMC X-1 optical coordinates](https://simbad.harvard.edu/simbad/sim-id?Ident=X+SMC+X-1) · [Swift Small Magellanic Cloud Survey](https://www.swift.psu.edu/S-CUBED/)

Offline replay requires NumPy and SciPy. The XRT script also uses Astropy. Both JSON files can be regenerated from the pinned tables. To independently extract from the mission products, supply the original MAXI orbit file from the source manifest, or use XRT's extraction option with the linked files online or in a local directory.

```sh
python energy_bands.py --output /tmp/energy-bands.json
diff energy-bands.json /tmp/energy-bands.json
python xrt_pointings.py --output /tmp/xrt-counts.json
diff xrt-counts.json /tmp/xrt-counts.json
python energy_bands.py --raw-orbit /path/to/J0117-734_g_lc_1orb_all.dat --data /tmp/maxi-orbit-bands.csv.gz --output /tmp/maxi-rebuilt.json
python xrt_pointings.py --extract --events /tmp/xrt-events.csv --output /tmp/xrt-rebuilt.json
```

## A dated test

Using MAXI data ending November 3, 2026 (last retained sample MJD 61347.5), the full late window predicts the next broad flux minimum at MJD 61363.95, approximately November 20. Starting the fit later and allowing a second harmonic moves the forecast as late as November 23. These are model predictions with the same data cutoff. A new well-covered trough outside November 16–29 would miss our deliberately wider observing window for a continuing 58–60-day rhythm. SMC X-1 has already demonstrated that it can change that rhythm.

## Inputs and replay

The compact CSV contains the eligible daily rates and their published errors, the four-day 12–25 keV GBM points at their bin midpoints, and phase-selected daily aggregates from individual Swift pointings and MAXI orbit bins. Time is Modified Julian Date. BAT daily labels are placed at day center. For the broad orbital phase cut we approximate Swift MET as days since MJD 51910, ignoring second-scale clock corrections. MAXI times are ISS bin centers. Neither approximation affects a timing measurement on this scale. We require finite rates, positive errors below four times the median valid error of that product, and good Swift flags. Negative rates are retained. Each product's full SHA-256 digest, row counts, cutoff and source URL are in the manifest. The linked monitoring streams can acquire new data or revised values. [Pinned inputs](data/measurements.csv) · [Source manifest](data/sources.json) · [All fits, failures and ranks](measurements.json)

From this directory, with Python, NumPy and SciPy installed, run the first two commands for an offline numerical replay. Matplotlib is needed only for the figure. To request new mission files or rebuild from saved originals, use the extraction commands with a separate output directory.

```sh
python replay.py --output /tmp/smc-x1-replay.json
diff measurements.json /tmp/smc-x1-replay.json
python plot.py
python extract.py --output /tmp/smc-x1-new-download
python extract.py --raw-dir /path/to/originals --output /tmp/smc-x1-from-originals
```

Swift/BAT transient-monitor measurements are provided by the Swift/BAT team. MAXI data are provided by RIKEN, JAXA and the MAXI team. The signal directions are J0117−734 (SMC X-1) and J0532−663 (LMC X-4). [Swift/BAT monitor documentation](https://swift.gsfc.nasa.gov/results/transients/Transient_synopsis.html) · [MAXI data columns and time system](https://maxi.riken.jp/top/lc_readme.txt)
