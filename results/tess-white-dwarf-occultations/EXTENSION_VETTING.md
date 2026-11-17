# None of the 104 extension triggers establishes a white-dwarf occultation

The original aperture measurements reproduce all 19 periodic and 85 isolated trigger scores on 20 white dwarfs. No periodic ephemeris met the original detection criteria in either simultaneous SPOC light curve. Of 85 isolated boxes, 83 had zero quality-zero finite 20-second flux measurements in the central window. The two with usable fast data do not show a white-dwarf-localized deficit of the reported depth in the 120-second pixels. The rate analysis retains all compatible original periodic triggers in its numerator.

## Completed checks

| Fixed-ephemeris reduction | Periodic covered / above threshold | Isolated covered / above threshold |
|---|---:|---:|
| Published FFI aperture | 19 / 19 | 85 / 85 |
| Same FFI pixels and aperture, with fast quality mask | 15 / 0 | 2 / 2 |
| FFI target profile and quadratic sky, with fast quality mask | 15 / 0 | 2 / 0 |
| Simultaneous 20-second SPOC SAP | 16 / 0 | 2 / 0 |
| Simultaneous 120-second SPOC SAP | 16 / 0 | 1 / 0 |

“Covered” requires the original cadence-fraction requirements in the central box and each flank, in at least three cycles for a periodic signal. The table uses the saved ephemerides, including their durations, without refitting at another cadence. The fast quality mask requires at least eight quality-zero, finite 20-second measurements centered in each 200-second exposure. This retains 196,795 of 234,992 quality-zero FFIs on these 20 trigger-host visits. Among the 85 isolated boxes, all time-tagged in-window 20-second exposures carry bit 4096 in 74 boxes; eight other boxes carry other flags throughout, and one mixes exposures carrying 4096 and another exclusion flag. Seventy-five windows contain at least one exposure carrying bit 4096. Loss of coverage leaves the physical origin of individual events open.

The alternative FFI extraction fits the target with a pixel-integrated Gaussian of width 0.85 pixels, fits a quadratic spatial background, and includes sufficiently bright, separated Gaia neighbours on the stamp. All reductions use the original six-hour normalization, two local flanks, duration-specific noise scale and recurrence test. The products at different cadences share the same sky and spacecraft. The simultaneous comparisons use individual spacecraft exposure times and the FFI barycentric correction. The archived checks include relative flux, local sky, centroids, nearest comparison stars and individual periodic cycle coverage.

| Isolated event | Native FFI depth | SPOC 20-second depth | 120-second empirical WD-profile depth | Matched-control scatter |
|---|---:|---:|---:|---:|
| TIC 241392886, S108, BTJD 4310.27510, 4 min | 4.41% | −0.28% | −0.935% | 0.550% (1,326 controls) |
| TIC 288144896, S84, BTJD 3607.04999, 64 min | 1.33% | 0.73% | 0.133% | 0.157% (118 controls) |

The Sector 108 box has only one quality-zero central two-minute exposure. In the better sampled Sector 84 box, 33 central two-minute exposures show a 0.717% loss summed over the 1.5-pixel-radius aperture. The contrast is spread across the stamp. A fit to 19 pixels within 2.5 pixels of the propagated WD, using its measured out-of-window spatial profile and a local sky plane, assigns 0.661 electrons/s/pixel to the constant sky contrast. The fitted WD loss lies inside the central 95% range (−0.302% to +0.339%) of the disjoint off-time control contrasts. The profile and controls are computed by the same method for Sector 108. Controls preserve the number and spacecraft-cadence pattern of retained central and flank samples. The scatter describes empirical background variation and is not an occurrence limit or an independent detection significance. [120-second input hashes, difference images, source templates and controls](extension/extension_supported_pixels.json)

All 166 SPOC comparison light curves were obtained and identified from FITS headers. Forty follow the targets in the discovery visit (two cadences per target), 60 measure three other two-minute targets on the same CCD per field, and 66 sample other visits of the periodic trigger hosts. Comparison stars more than 120 arcsec away test shared detector behaviour. Gaia neighbours inside each science aperture are recorded separately in the frozen FFI measurements. Period and phase errors can accumulate between seasons; distant fixed-phase checks have no calibrated exclusion power here. [SPOC product URLs, hashes, quality counts and measurements](extension/extension_replications.json.gz) · [Every fixed-ephemeris score](extension/extension_fixed_ephemerides.json.gz) · [Original 9-by-9 difference images and centroid shifts](extension/extension_image_checks.json.gz) · [Gaia neighbour coordinates and dilution](extension/extension_measurements.json)

## Replay

Start in this result directory with NumPy, SciPy, Astropy and tess-point installed. The first command checks the frozen TIME, FLUX and QUALITY SHA-256 values for each FFI, acquires the comparison products and packs the evidence. Reuse an empty output directory for a fresh measurement. The remaining commands use that same cache.

```sh
python replication.py extension/extension_measurements.json extension/extension_parent.csv cache target_pixel_products.csv cache replication_records --receipts extension/extension_receipts.json --bulk bulk_receipts.json --pack extension/extension_replications.json.gz
python vet_extension_pixels.py extension/extension_measurements.json target_pixel_products.csv cache extension/extension_fixed_ephemerides.json.gz
python localize_supported.py extension/extension_fixed_ephemerides.json.gz extension/extension_parent.csv target_pixel_products.csv cache extension/extension_supported_pixels.json
```

## Remaining limits

The source of an individual quality-flagged event remains unresolved. An astrophysical occultation could have occurred during a spacecraft exclusion. This fixed-ephemeris comparison has no completeness calibration for a drifting period, a different-duration feature or aperiodic debris; the separate injection experiments quantify their stated periodic models. External catalogue and ground-based coverage are documented separately in the survey record.
