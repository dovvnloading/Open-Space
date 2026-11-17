# Native 20-second pixels do not validate a short occultation

All 124 fast light-curve triggers have usable central and shoulder exposures in the raw target-pixel products. These 22 star-sector visits cover 16 white dwarfs. In every 1.5-pixel-radius aperture, the mean brightening of the shoulders above the local baseline exceeded the loss in the central box. Eleven windows lacked a complete central box or shoulder in the 120-second pixel remeasurement. At 20 seconds, those eleven yield central changes from a 0.87% brightening to a 1.95% loss in a stellar-template fit. Their flanks are 13.05–104.43% brighter than the local reference in the same fit. This explains the large box contrasts without requiring a deep occultation. [FITS and Gaia input hashes, central exposures, flank counts, aperture depths and matched controls](short_cadence/fast_trigger_pixel_checks.json.gz)

The largest apparent central loss in any stellar-template fit was 6.63% on TIC 14802607 at BTJD 2544.66849 in Sector 45. A positive source traverses the stamp during that visit. Twenty-one bright, off-edge residual positions follow a linear track at 39.43 pixels/day with 0.029-pixel radial RMS. The source stood 2.07 pixels from the WD at the trigger epoch and subsequently crossed its aperture. Relative to the quiet stamp 0.28–0.55 days before and after the box, the central aperture was **0.54% brighter**. Simultaneous fits of the moving source, the stationary WD image and the local sky yield a WD loss of 0.86% at the nominal Gaussian width of 0.85 pixels. Across widths 0.65–1.05 pixels the fitted loss ranges from −0.35% to +2.85%. The surrounding hours contain extra moving light and do not supply an appropriate undisturbed baseline. No orbital identity was assigned to the moving source. [Time-resolved track, pixel fit and archival product hash](short_cadence/fast_moving_source.json)

TIC 8988748 (GD 1212) had central stellar-template losses of 2.31% in Sector 42 and 2.24% in Sector 70, with flanks 6.01% and 6.47% above their local baselines. The windows occur in independently measured pulsating light curves. Neither event supplies an independently recurring transit ephemeris. Smaller intrinsic or blended signals can coexist with the variability. [Measured pulsation periods](short_cadence/pulsation_checks.json)

## Measurements and limits

For each fast trigger, the pixel coordinates follow the Gaia proper motion to its sector epoch. Quality-zero calibrated FLUX values support two extractions. The first sums a circle of radius 1.5 pixels. The second fits the measured, out-of-window WD image in a disk of radius 2.5 pixels together with a varying sky level and two slopes. Both retain the original central box and two distinct flanks. The context reference is the median outside twice the duration and within the greater of 0.12 days and five durations. For 113 boxes, disjoint off-time comparisons translate the same surviving spacecraft cadence-number pattern. For 11, the controls match counts separately within each box and flank. The counts of equally deep control windows are descriptive after trigger selection. These data share photons with the original SPOC light curves.

Run from the result directory, using the directory of 120-second output records from `roster_pixels.py` (or their JSON list compressed with gzip). The fast-product run needs approximately 6 GB of cache space.

```sh
python fast_trigger_pixels.py two_minute_records cache short_cadence/fast_trigger_pixel_checks.json.gz
python fast_moving_source.py short_cadence/fast_trigger_pixel_checks.json.gz cache short_cadence/fast_moving_source.json
```

This follow-up tests only the original isolated boxes. No sensitivity calibration for a time-variable debris population is available.
