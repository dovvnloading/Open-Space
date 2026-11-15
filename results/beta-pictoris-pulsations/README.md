# A reversal in the pulsation amplitudes of β Pictoris

TESS measured a reversal in the relative amplitudes of two resolved frequency bands in β Pictoris (HD 39060, TIC 270577175). The group around 47.28 cycles/day remained weaker than the 47.44 cycles/day peak in sector 87 (December 2024–January 2025). It was stronger in sector 95 (July–August 2025) and reached 2.36 times the 47.44 amplitude in sector 98 (November 2025–January 2026). The cadence coverage brackets the change in ordering. The data cannot date a unique crossing inside the gap.

| TESS sector | 47.28 group, ppm | 47.44 peak, ppm | Amplitude ratio | Ratio with the published close pair fitted separately |
|---|---:|---:|---:|---:|
| 32 | 279 | 1,071 | 0.26 | 0.27 |
| 34 | 419 | 1,024 | 0.41 | 0.41 |
| 87 | 465 | 553 | 0.84 | 0.84 |
| 95 | 452 | 357 | 1.27 | 1.26 |
| 98 | 554 | 234 | 2.36 | 2.36 |

The first three numeric columns use sinusoidal semi-amplitudes measured in PDCSAP. The last compares the sampled RMS of the reconstructed 47.2686 + 47.2853 cycles/day signal to that of the 47.44 peak. The close pair is narrower than the reciprocal of a typical sector duration. A changing mix of unresolved components, including beating, is a viable explanation for changes within either band. The measurements establish a change in the observed photometric spectrum. They do not measure energy transfer between individual eigenmodes.

Zieba et al. reported a developing 47.3 cycles/day feature in 2018–21 and tabulated its close components. Lecavelier des Etangs et al. described the dominant 47.44 cycles/day signal in the earlier data. The 2025 reversal and its persistence through the sector ending in 2026 are the new measurements. A targeted search for β Pic, HD 39060, the frequencies, sectors 87/95/98 and an amplitude reversal retrieved the earlier reports. Priority remains open. [Zieba et al. 2024, §5.5, Fig. 15 and Table A.1](https://www.aanda.org/articles/aa/pdf/2024/07/aa47754-23.pdf) · [Lecavelier des Etangs et al. 2022](https://doi.org/10.1038/s41598-022-09021-2) · [Zwintz et al. 2019](https://doi.org/10.1051/0004-6361/201834744)

## Checks on the amplitude reversal

- Each sector is fit at the native two-minute cadences with quality zero, after subtracting a linear baseline in each half-day interval and fitting neighboring pulsations jointly. Fits with each frequency free within ±0.025 cycles/day agree. The 47.28 and 47.44 peaks are separated by about 0.155 cycles/day.
- Both halves of sector 87 have ratios below one (0.79, 0.88). Both halves of sectors 95 (1.13, 1.46) and 98 (2.40, 2.26) have ratios above one. Raw-aperture SAP and corrected PDCSAP give the same ordering. These reductions share flux measurements. The separated sectors are the temporal replication.
- In two target-pixel stamps the two bands agree in their signed flux-weighted detector-column centroid to 0.022 pixel (sector 87) and 0.002 pixel (sector 95). Their profiles follow the separately known 50.49 cycles/day pulsation. This test sums the charge bleed along detector columns on this saturated star.
- The MAST target inventory and the TESScut footprint query are recorded in [coverage.json](coverage.json). The ten analyzed sectors exhaust the archived two-minute target light curves returned by the inventory. TESScut also reports full-frame coverage in sector 68. The original exocomet screen recovered the previously published sector 6 dip at BTJD 1486.4; the result reported here concerns stellar variability.

## A 15.28-hour signal in the same aperture

There is a second, weaker periodicity. The maximum of a joint 1.48–1.66 cycles/day search is 1.570639 cycles/day (15.2804 h, 29 ppm). At this frequency, target PDCSAP semi-amplitudes are 50.8, 41.6 and 38.2 ppm in sectors 87, 95 and 98. The corresponding SAP amplitudes are 54.1, 43.1 and 39.1 ppm. At a common BTJD 3000 epoch, the ten individual fitted phases span 1.16 radians. The three separated late sectors agree in phase to 0.36 radian. These sectors were examined during discovery, so these comparisons are retrospective replication.

The spectral window admits competing aliases. The next joint peaks at 1.571825 and 1.567886 cycles/day reach 24.0 and 23.6 ppm. The quoted digits specify the grid maximum. A local 1.57 cycles/day peak also recurs in each of the three widely separated observing seasons. The simple difference between the prominent 50.492 and 48.919 cycles/day pulsations is 1.573 cycles/day and drifts several phase cycles across the data relative to the preferred periodicity. No unique physical identification follows from the period.

The target's fixed-frequency amplitude divided by its median amplitude in the 0.85–2.3 cycles/day comparison band (excluding ±0.12 around the target frequency) is 4.6–5.7 in the three late sectors. A comparison on the same camera and CCD at each epoch produced ratios at most 2.2 in 18 control light curves. The target background and both spacecraft pointing-correction components have ratios below 1.7 at all three epochs. Controls were selected from the nearest fourteen other numeric TIC identifiers beginning 270 in each sector's SPOC target list, then filtered to the same detector and a separation above 0.04 degrees and ranked by sky separation. The six closest eligible products in each sector appear in [control_manifest.csv](control_manifest.csv). This convenience sample tests common detector behavior. The sideband ratio is descriptive and gives no search-wide false-alarm probability.

**Source attribution is unresolved.** In sector 87, the slow signal has its signed column centroid at 7.49 versus 6.02 for the 50.49 stellar pulsation. In sector 95 the comparison is 4.10 versus 5.32. Several bright neighboring detector columns contribute in opposite phase. These signed centroids are diagnostics of flux redistribution and cannot serve as sky positions. A blended variable, charge motion in the saturated image, and a source-dependent instrumental effect remain viable. We do not assign a stellar rotation period, a new mode, a planet, or an exocomet to the 15.28-hour modulation. Prior studies discuss the rotation and low-frequency behavior of this star. [Zwintz et al. 2019](https://doi.org/10.1051/0004-6361/201834744) · [Zieba et al. 2024](https://doi.org/10.1051/0004-6361/202347754)

## Reproduction

Run from this directory with Python 3, NumPy, SciPy and Astropy. The programs verify the manifest SHA-256 values before using cached files. The two optional pixel products together occupy approximately 356 MiB. All photometry comes from TESS/SPOC via MAST. No derived photometry is fabricated or injected.

```sh
OPENBLAS_NUM_THREADS=1 python measure.py --cache /tmp/beta-pic-tess
OPENBLAS_NUM_THREADS=1 python slow_signal.py --cache /tmp/beta-pic-tess
OPENBLAS_NUM_THREADS=1 python pixel_check.py --cache /tmp/beta-pic-tess
```

[Sector measurements](measurements.csv) · [Frequency aliases](frequency_aliases.csv) · [Target and field checks](slow_checks.csv) · [Signed pixel-column profiles](pixel_checks.csv) · [MAST TESS science products](https://archive.stsci.edu/missions-and-data/tess)
