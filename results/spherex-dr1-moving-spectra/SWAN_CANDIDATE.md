# SWAN C/2025 R2: preliminary repeat spectral-excess lead

This is an unvalidated candidate for further analysis, **not a discovery, a calibrated spectrum, or a new molecular identification**. Published 2026-11-16.

## Selection and provenance
A MOST query returned 348 intersections. After the faint Echeclus pilot, select the first nine D5 records in manifest order as a brighter-target feasibility screen, before examining their flux. `most-swan.tbl` freezes all returned records; `pilot-swan.json` freezes the first pass (one HTTP503), and `pilot-swan-retry.json` preserves a second attempt that retrieved the missing product. The script's cache now includes the complete URL with target center/size, avoiding cross-target cache collisions. All retrieved pilot images report quality Pass.

## What was measured
Raw 2.5-pixel aperture sums of MJy/sr image values, OUTLIER-inclusive diagnostic mask:

| UTC date | wavelength, microns | raw sum | formal error |
|---|---:|---:|---:|
| 2025-12-22 16:51 | 4.29284 | 13.637 | 0.304 |
| 2025-12-22 16:53 | 4.32943 | 2.762 | 0.222 |
| 2025-12-22 16:55 | 4.36619 | 1.931 | 0.232 |
| 2025-12-25 19:46 | 4.21545 | 5.701 | 0.251 |
| 2025-12-25 19:48 | 4.25171 | 17.941 | 0.316 |
| 2025-12-25 19:50 | 4.28817 | 13.440 | 0.296 |

The pattern repeats near the expected CO2 feature across two dates, each with adjacent observations separated by roughly two minutes. It merits testing, but every December 2.5-pixel strict-mask aperture fails because source pixels are flagged. OUTLIER inclusion is diagnostic only. The recovered Dec22 4.25660-micron aperture fails even the inclusive gate because its core also contains bit28 (and other flags); it is excluded, not a nondetection. The October continuum source has a bright off-target control and one bad core, so it is not used as a clean epoch comparison.

## Required falsification before any science claim
- Independently confirm ephemeris/centroid tracking and observation midpoint treatment.
- Inspect images and neighboring-source contamination, not just scalar sums; remove stationary-source contamination using independent track positions/epochs.
- Integrate the full spectral response, PSF/aperture and solid-angle corrections; verify the broad-wing response near 4.29 microns before assigning a molecule.
- Demonstrate feature stability to background/aperture choices and outlier handling with independent exposures.
- Account for correlated noise, color corrections, extended coma, and selection/trials. Formal errors above are not discovery significance.
- Complete literature search under C/2025 R2, SWAN, SWAN25B and instrument aliases. Absence of a search result does not prove novelty.

The separate Wierzchos C/2024 E1 pilot also shows a near4.26-micron excess, but CO2 emission is already published with JWST (Snodgrass et al., DOI10.1093/mnrasl/slaf046). It cannot be promoted as a new discovery. The SWAN finding is currently an exploratory lead only.

Reproduce the deterministic pilot:

```sh
python examples/spherex_moving_screen.py results/spherex-dr1-moving-spectra/most-swan.tbl --detectors D5 --limit 9 --output /tmp/swan.json
```
