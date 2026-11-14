# Where does SkyMapper see the eclipses?

Independent Southern Sky Survey images place the variable light close to the catalogued source, SMSS J173806.43−492057.1 (SkyMapper DR4 object 363533777). The known eclipse waveform is recovered in aperture-corrected magnitudes measured in a **2-arcsecond diameter** aperture. A separate astrometric measurement finds no significant motion of the image photocentre during eclipse. The precision is hundreds of milliarcseconds on the location of the variable light. The proposed outer relative orbit has a semimajor axis of only about 2.2 mas for the fixed TESS timing scale. [SkyMapper DR4](https://skymapper.anu.edu.au/data-release/dr4/) · [Data release paper](https://doi.org/10.1017/pasa.2024.53) · [Camera and archive fields](https://skymapper.anu.edu.au/table-browser-standalone/dr4/) · [Gaia orbital interpretation](../validation/spectroscopy/README.md)

## Measurement

SkyMapper recorded 58 target detections. We kept 28 exposures on seven nights in 2015–2021 across all six filters, requiring zero source-extraction flags, zero image-mask flags, PSF magnitude error below 0.03 mag, and at least three same-image, same-chip reference stars. Eligible frames had mean seeing FWHM 1.49–4.97 arcsec. We used the previously published Gaia G-band two-Gaussian-plus-ellipsoidal shape, period and epoch, with UTC exposure midpoints converted to barycentric JD(TCB). The fits allow an independent zero point for each filter. No outer-orbit prediction enters this localization.

| Brightness estimator | Fitted eclipse waveform amplitude relative to Gaia G |
| --- | ---: |
| Survey PSF magnitude | 1.05 |
| 2-arcsec aperture diameter | 0.90 |
| 5-arcsec aperture diameter | 1.01 |
| 8-arcsec aperture diameter | 1.04 |

The mean image position of the bright target agrees with Gaia to 0.02 arcsec in SkyMapper's master crossmatch. For the epoch test, we used a 0.15° local reference field. Reference stars have catalogued g between 10.5 and 15.8, no master extraction flags, nearest SkyMapper neighbor farther than 4 arcsec, Gaia match within 0.4 arcsec, Gaia RUWE below 1.4, multi-peak fraction at most 1%, and clean frame photometry with error below 0.08 mag. After propagating Gaia DR3 positions for proper motion and parallax, we subtracted each exposure's median astrometric residual from 6–29 comparison stars on the target CCD. [Observation log and exact query digests](observations-manifest.csv) · [Epoch results](localized-epochs.csv)

For variable light lost at a fixed position, the change in centre of light is its displacement from the bright photocentre with reversed sign, multiplied by the fractional light loss divided by the remaining flux. We regress both coordinates against the inverse relative flux minus one from the Gaia template and normalize the slope to the measured PSF depth. For all filters together, the inferred position of the removed light relative to the bright photocentre is **−0.28 ± 0.38 arcsec east and −0.31 ± 0.42 arcsec north**. The errors are seven-night jackknife standard errors, recomputing both the centroid slope and its depth normalization on each omission. Each coordinate is consistent with zero displacement. This describes the source of the short-period modulation, including any associated ellipsoidal variation. [Fitted quantities](localization.json)

A comparison set of 27 field stars on at least 15 of the same accepted frames and at least five nights gives PSF magnitude slopes from −0.31 to +1.21 versus +1.08 for the target. One comparison star exceeds the target's positive slope. These sparse dates alone would be insufficient to identify an unknown binary. The pooled photocentre fit also assumes similar fractional modulation across filters. Chromatic refraction, seeing, catalog centroid weighting, possible intrinsic reference-star variability and covariance of season with eclipse phase limit the spatial interpretation. The aperture is centered on each detection and has seeing-dependent throughput. It supplies no hard 1-arcsec exclusion radius or completeness curve. [Comparison controls](comparison-controls.csv)

The ExoFOP target page has zero imaging and spectroscopy entries. The Gemini 60-arcsec pointing query and the SOAR 2023 catalog query yielded no rows. We also checked all 5,316 resolved measurements and 2,376 unresolved observations in the SOAR 2024–25 machine-readable tables. The closest listed observation is 0.72° away. Within the hour-17 USNO Fourth Interferometric Catalog mirror we found no complete position header within 60 arcsec. These are checks for archived observations. They supply no diffraction-limited upper limit for this particular star. [Archive receipts](imaging-search.json) · [SOAR 2024–25](https://doi.org/10.3847/1538-3881/ae6cf5)

## Replay

Python, NumPy, pandas and Astropy suffice. From this directory run the following, supplying the path to the original target dossier. The compact input contains the analysis columns for target and field stars at the target's exposure IDs. The manifest records checksums of the original untrimmed archive responses. Positions and measured fluxes come from SkyMapper DR4 and Gaia DR3. The positional fits and their interpretation are derived here.

```sh
OPENBLAS_NUM_THREADS=1 python localize.py /path/to/Open-Space/results/tic-96100077
```
