# MACS J0416.1−2403, January 2026 against February 2023

A blind, same-filter subtraction around the Warhol and Spock arcs found **zero same-sign ≥5σ F150W/F444W peak pairs within 0.15″** in the **0.190136 arcmin²** common usable area. The screen detected two F150W peaks and six F444W peaks separately. The historically variable Warhol arc shows a repeatable F444W excess at a position confused with published microlensed stars.

## Search depth

The 36.06″ square at (64.0375°, −24.0686°) uses public NIRCam mosaics from GO/GTO 1176 and 6882. The older data were obtained on 2023 February 10 and the newer data on 2026 January 10. Each pixel must have ≥2 contributing exposures per epoch, positive weight and finite science and error values. The pre-existing El Gordo analysis fits affine alignment, measures compact-source PSFs, cross-convolves the images, removes smooth structure and evaluates the signed matched-filter response with local robust noise. Bright structure, edges and low exposure coverage reduce the usable area. F150W and F444W registration RMS were 0.0092″ (54 anchors) and 0.0065″ (62 anchors), respectively.

| Injected AB magnitude in each band | Recovered / injected, brightening | Recovered / injected, fading |
|---|---:|---:|
| 26.0 | 90/90 | 90/90 |
| 26.5 | 59/90 | 62/90 |
| 27.0 | 11/90 | 8/90 |
| 27.5 | 1/90 | 0/90 |

Linear interpolation of the pooled recovery fractions places 50% recovery at **AB 26.65**. The 2,700 signed trials cover AB 23–30 with spatially distributed equal-flux point sources. Completeness is conditional on the registered image, fixed PSF model and unmasked area. It excludes sources lost to bright-host masking, variables with different changes in the two bands, and variability that occurs between the sampled epochs. [All injection bins](injections.csv) · [Area and alignment](screen.json) · [Single-band peaks](single_band_peaks.csv) · [Two-band detections](joint_detections.csv)

## The Warhol residual

The F444W peak at **(64.03667859°, −24.06730000°)** is 0.072″ from the published position of D31-W5 and 0.118″ from D21-W3, in the spectroscopically identified z = 0.9397 Warhol arc. Both separations are smaller than the F444W PSF FWHM. Relative astrometric calibration of the historical catalogs adds further uncertainty. Distinct stellar components are unresolved here. [Yan et al., Table 2](https://doi.org/10.3847/1538-4365/ad0298) · [Warhol arc redshift](https://doi.org/10.1051/0004-6361/202349119)

| Filter | Earlier UTC epoch | January 2026 minus earlier flux (nJy) | Empirical S/N |
|---|---|---:|---:|
| F444W | 2022 October 7 | +80.7 ± 9.1 | 8.84 |
| F444W | 2022 December 29 | +63.4 ± 10.6 | 5.99 |
| F444W | 2023 January 12 | +43.0 ± 9.2 | 4.69 |
| F444W | 2023 February 10 | +77.5 ± 11.6 | 6.65 |
| F150W | 2023 February 10 | +13.2 ± 22.1 | 0.60 |
| F300M | 2024 January 17 | +47.2 ± 15.1 | 3.12 |

These are forced difference fluxes at the fixed F444W peak, including unresolved arc light and subtraction systematics. The four F444W differences share one late image. The peak also coincides with caustic-sensitive structure. The nearest zero of the tangential eigenvalue in the public CATS v4.1 model is **0.410″** away, rescaling its convergence and shear to zₛ = 0.9397 and zₗ = 0.396 with Ωₘ = 0.3 in flat ΛCDM. This is one macro-model and the number is model dependent. Source confusion precludes assigning an individual-star light curve, absolute magnitude, or discovery classification. A direct TNS cone request returned HTTP 403. The positional overlap with published sources independently rules out claiming a previously unreported source here. [All measured differences](forced_differences.csv) · [Model inputs and calculation](warhol_lens.json) · [HFF model archive](https://archive.stsci.edu/pub/hlsp/frontier/macs0416/models/cats/v4.1/)

Independent F444W dither splits reproduce the brightening. Four single-exposure I2D products per epoch cover the position on detector NRCBLONG. Each was aligned to its own epoch using surrounding structure with the central 0.66″ excluded from the fit. Comparing odd dithers from both epochs gives **+93.1 ± 13.5 nJy (6.90σ)**. The disjoint even comparison gives **+63.7 ± 10.8 nJy (5.92σ)**. All four individual dither pairs have positive differences. The error scale is the spatial MAD of the filtered difference in a masked local annulus. Shared mosaic PSFs and coherent residuals on the underlying structured arc can affect both tests. A single bad exposure cannot account for this pattern. [Split measurements](warhol_split.csv) · [Eight exposure products and byte-range receipts](dither_ranges.json)

The 48.9σ F150W peak comes from a compact negative pixel defect in the 2023 mosaic. In a 21×21 native-pixel box the 2023 minimum was −6.80 MJy/sr against a local median of +0.223, while the 2026 minimum was +0.151 against a median of +0.201. The ±47σ F444W pair is an opposite-sign positional dipole around one bright point source, also visible in F150W. The remaining low-significance peaks lie on an extended source, bright-source wings or masked structure, or lack the second-band detection. No pixel-level classification of an isolated supernova or individual star follows from them.

## Reproduce

With NumPy ≥2, SciPy, Astropy, and contourpy, run the following in this directory. The first two commands run the blind analysis without any follow-up epochs. Download receipts record public science-product byte ranges, ETags, and SHA-256 hashes. No external source catalog is used for blind selection.

```sh
python macs0416.py fetch cache
python macs0416.py measure cache
python epoch_context.py fetch cache
python epoch_context.py measure cache
python split_dithers.py cache
python lens_curve.py
```

[Primary FITS receipts](input_ranges.json) · [Additional epoch receipts](extra_ranges.json)
