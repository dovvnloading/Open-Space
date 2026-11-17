# C/2025 V1 on October 27 with a trailed point-source template

The two clean ZTF difference images do not establish an earlier detection. The joint template reaches 3.15 background σ at (+2, −2) arcsec from the predicted position. Nine of 168 offset tracks produce as large a maximum, giving a rank probability of 10/169. The exposures individually reach 3.27σ and 0.96σ at that common offset. The [same-night subtraction](../science-pairs/README.md) yields a broad maximum at (+3, −2) arcsec with rank probability 0.24. The [aperture search](../stack/README.md) prefers (−1, +3) arcsec. Each sampled the Horizons three-sigma region.

| Palomar I41 exposure | Midpoint UTC, October 27 | Conditional Gaia-corrected RA, Dec in degrees | East, north residual in arcsec | Formal flux S/N | WCS tie |
| --- | --- | --- | --- | ---: | ---: |
| 20251027542164 | 13:01:05.003 | 170.77861829, −10.09786955 | +2.51, −2.61 | 2.96 | 31 stars |
| 20251027547164 | 13:08:17.003 | 170.78560267, −10.09372623 | −0.50, +0.16 | 2.98 | 21 stars |

The conditional coordinate errors are approximately (0.67, 0.61) and (0.61, 0.55) arcsec, using the formal pixel-coordinate errors, near-unity plate scale and measured single-star dispersion in quadrature. Local shifts to Gaia DR3 are (−0.021, −0.094) and (−0.095, −0.079) arcsec. Robust single-star dispersions are (0.072, 0.081) and (0.065, 0.073) arcsec. The residuals change by 4.1 arcsec in seven minutes. Peak selection and the chance of fitting sky structure dominate the conditional errors. Neither coordinate entered the main 117-position solution.

The known October 29 control peaks at zero common offset with 5.93σ. Its noisier second exposure contributes 0.88σ there, while the first gives 6.42σ. The October 27 local five-sigma flux scale is 158 DN at zero point 25 (r approximately 19.50) for this trailed point-source model. Measuring completeness at that magnitude would require artificial-source recovery.

The [script](../../centroid_check.py) convolves a Gaussian PSF with the apparent trail during one exposure and fits a local background plane under each 17 × 17 pixel stamp. The [measurements and receipts](measurements.json) preserve midpoints, calibration, controls, FITS hashes and URLs. Both compressed Gaia queries retain catalog coordinates and proper motions.
