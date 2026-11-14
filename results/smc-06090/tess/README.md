# TESS measured the timing after the final OGLE observation

The public [MAST TESScut service](https://mast.stsci.edu/tesscut/) lists ten sectors at RA 21.06275°, Dec −73.32594444° in the SMC field. NASA’s TESS mission provided these full-frame images. See [Brasseur et al. (2019)](https://mast.stsci.edu/tesscut/) for Astrocut and the [MAST data description](https://outerspace.stsci.edu/spaces/TESS/pages/14563420/2.0%2B-%2BData%2BProduct%2BOverview) for time stamps and background conventions. The archive and data credits belong to NASA, TESS and STScI/MAST. The measurements here come from a new reduction. Exact URLs and numerical input hashes are in tess_inputs.json.

## Data and extraction

Every inventory sector was requested as an 11 × 11 pixel cutout. There were 66,116 input images; the numerical receipt gives the individual counts. We required quality zero, finite times and complete finite pixel stamps. For the main extraction, six common modes of pixels more than 3.1 pixels from the catalog position, plus an intercept and slope in each 0.65-day block, model variations unrelated to the target. The combined background residual across outer pixels supplies an additional robust cadence cut. The primary aperture is the 3 × 3 pixels around the nearest catalog pixel. The other apertures use the nearest four pixels or one pixel. These apertures share data. A separate reduction in alternate_reduction.py uses the fainter outer pixels as a sky estimate and a segmentwise running median; it fits no detector modes. In these residual fluxes the measured depths depend on the reduction and aperture and are not intrinsic stellar depths.

We fold on the published OGLE-IV period and epoch (1.1696719 d and HJD 2460000.0369). An approximate 70-second translation from BJD_TDB to HJD_UTC is used. TESS times are corrected to CCD center. This approximation mainly affects the common phase of a given sector. For the free-center timing comparison, we used 120 phase bins, a Gaussian width of 0.025 cycle per eclipse, free depths and centers, and a sinusoidal baseline at the orbital frequency. Changing both widths independently is a sensitivity test. The bins average over real timestamps; finite exposures and binning can affect measured shapes and have little effect on the relative timing of symmetric dips. The spatial fits use axis-aligned local elliptical Gaussians plus a constant, and serve only as approximate centroid checks in crowded pixels.

## Timing and controls

selection.md records the fixed templates and the selection of sectors 94–96 as the external time check before the cutout flux was inspected. The free-center timing fit and its data-reduction sensitivities are exploratory. The file tess_timing.json records the outcomes for all four time groups. A 600-replicate bootstrap resamples orbital cycles independently within each sector. A second bootstrap resamples sets of up to three consecutive accepted cycles within sectors for the earliest-to-latest difference. Surrogate and comparison results apply to this target and do not estimate a sky-wide false alarm rate.

For the fixed-depth test in 2025, there were 49 well-covered individual cycles in the per-orbit regressions, with positive primary and secondary depth estimates. For the free-center phase stack 51 cycles passed its distinct coverage cut. The mean of the primary and secondary fixed-template depths exceeded all 63 tested common phase displacements in the three aperture sizes. These shifts overlap and are not independent statistical trials. The largest raw fixed-template dimming pixel occurred away from the source in sector 1; the bounded local fit and its target-region result are recorded separately. In sectors 94–96 the strongest combined dimming pixel was the pixel nearest the target every time.

## Replay

Python 3.12, NumPy, SciPy and Astropy are used by the full reduction. To refit the four TESS intervals offline, run tess_timing.py from this directory. For spatial controls, fetch and reduce the original pixel data first so that folded_pixels.npz exists. The following commands then replay all three analyses. Orbit-only predictions in tess_timing.py use the sibling orbit/orbit.json.

```sh
OPENBLAS_NUM_THREADS=1 python photometry.py --from-folded
OPENBLAS_NUM_THREADS=1 python tess_timing.py
OPENBLAS_NUM_THREADS=1 python independent_checks.py
```

For the full pixel replay, run fetch.py first (about 158 MiB of decompressed FITS). Re-run photometry.py without --from-folded, then orbit_folds.py, tess_timing.py, independent_checks.py, alternate_reduction.py and injection.py. Compare the numerical hashes in the resulting tess_inputs.json against the original copy. The numerical hash covers the full original TIME, QUALITY and FLUX arrays, including rejected cadences; packaging metadata may change. The compact per-orbit aperture sums and counts are in ../joint/inputs/orbital_folds.npz. The pixel-level commands recreate the sector-folded images locally. Compare the recomputed per-orbit arrays with those distributed in ../joint/inputs before fitting. The synthetic injection is explicitly labeled in injection.json and is used only as a cadence and recovery check.

A synthetic signal was added to the actual calibrated pixels at two comparison phases with spacing fixed to 0.5300 cycle. The identical extraction recovered 0.5317 in 2018–2019 and 0.5313 in 2025. This single injected shape checks timing bias for this cadence and extraction. Its recovered depth was attenuated, especially for the earlier sectors. Across real TESS extraction choices, the late 2025 timing was stable; the earlier group shifted by about 14 minutes when the sky/trend reduction was changed. The quoted measurement intervals are conditional on the chosen extraction.
