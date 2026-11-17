# ATLAS photometry of the two possible giant-planet hosts

Independent ATLAS images cover three nights inside TIC 91529455's extrapolated 2.364306-day boxes. Twelve in-window exposures have a median deficit of 0.072% against the other 187. A bootstrap resampling independent observing nights gives a descriptive 95% interval of −0.34% to +1.33%. TIC 220578510's strongest 6.771955-day ephemeris has no in-window points. Both objects remain possible hosts in the TESS limit.

The tests cover all nine periodic triggers on TIC 220578510 and both on TIC 91529455. For TIC 220578510 the 0.966602-day fold gives a 0.048% deficit across 23 measurements on five nights. At 6.753547 days four in-window frames from one night give a median 5.14% deficit. The frame zero points range from 17.20 to 20.11 mag and their individual zero-point RMS values reach 0.196 mag. All four fail the secondary stricter calibration cut. At 6.820728 days six in-window points on two nights give a median 1.00% deficit. Two nights and a selected TESS alias supply limited evidence of recurrence. The 120-minute signals have TESS-epoch uncertainties that accumulate across the intervening years. No ATLAS fold establishes that a compact body crossed either WD.

## Measurements

1. Cone coordinates come from the archived Gentile Fusillo extension. CATCH returned 625 Sutherland footprints for TIC 220578510 and 1,030 for TIC 91529455 between 2022 January 1 and 2025 March 1. Select each star's latest 200 exposures by timestamp, independently of measured flux and TESS phase. Retrieve 240-arcsec calibrated FITS cutouts from the Small-Bodies Node Survey Image Service. The footprint search produced the input roster, and the FITS pixel values produced the measurements.
2. Use the primary TPV distortion solution and the proper-motion-propagated Gaia coordinates. Sum the signed sky-subtracted flux in a circular aperture of radius max(3 pixels, 1.55 FWHM) with a local annulus from max(9 pixels, 2.8 FWHM) to max(14 pixels, 4 FWHM). Estimated FWHM must lie in 1–6 pixels. The zero point must lie in 15–28 mag with at least 50 calibration stars and RMS below 0.2 mag. Reject edge, saturated and invalid apertures. Retain finite low and negative fluxes. Scale by exposure time and the frame zero point and normalize to each star's in-band median. The 193 and 199 accepted fluxes have median centroid offsets of 0.080 and 0.076 pixels. The archive products can include spatially varying calibration errors or unresolved neighbouring light.
3. Convert the shutter midpoint UTC MJD to barycentric TDB using the approximate Sutherland observatory position. Fold on each TESS period, epoch and full box duration without optimization. Report the difference of the out- and in-window median relative fluxes, the number of distinct nights, and a descriptive percentile interval from 1,000 night-bootstrap draws where two or more independent in-window nights exist. The per-fold intervals exclude period uncertainty and multiple testing. A secondary diagnostic applies zero-point RMS below 0.05 mag and zero point above 20.5 mag. Zero returned exposures in a box give no depth test. This is a finite archive sample with sparse phase coverage.

[Epochs, cutout URLs, calibration and FITS hashes](extension/atlas_sutherland/epochs.jsonl.gz) · [Query receipts, phase counts and all 11 folds](extension/atlas_sutherland/folds.json) · [CATCH survey holdings](https://catch.astro.umd.edu/about) · [SBN Survey Image Service](https://sbn-survey-image-service.readthedocs.io/en/latest/) · [ATLAS data products](https://pds-smallbodies.astro.umd.edu/data_sb/missions/atlas/ATLAS_PDS_Bundle_Overview.pdf)

Run from the search directory:

```bash
python atlas_epoch_pixels.py atlas-cache extension/atlas_sutherland
```
