# Complete control screen and independent-target sensitivity pilot

Status: exploratory screening, not a discovery. 2026-11-16.

All 90 D4–D6 29P rows in the frozen MOST manifest were attempted. 82 FITS cutouts were retrieved and analyzed; eight failed with HTTP 503. Of the 82 retrieved cutouts, 15 pass the strict 2.5-pixel aperture gate and 75 pass the explicitly diagnostic OUTLIER-inclusive gate. These are quality-gate counts, not detections. `screen-29p.csv` preserves every attempted row, exact times, wavelengths, cutout hashes, raw sums, formal errors, and failures. The earlier 51-row JSON remains a frozen intermediate checkpoint rather than a full final sample.

The repeat localized 4.26-micron feature described in CONTROL_RESULT.md survives the strict mask in the strongest exposure. It is a known-comet positive control, not novel chemistry. The May 2026 observations do not provide an equivalent near-4.26-micron comparison. Do not interpret unmatched seasons as gas variability.

## Echeclus sensitivity pilot

An independent 174P/Echeclus MOST query returned 257 intersections. Select the first eight D5 rows in manifest order, before examining image signal. Five downloads succeeded and three failed with HTTP503. The successful wavelengths are approximately 3.82–4.06 microns, not a CO2-band test. Strict valid 2.5-pixel aperture measurements have negative, low-significance signal; some apertures fail quality masking. This pilot does not establish gas absence or a calibrated upper limit, and it is not evidence of an astrophysical discovery. Do not use the failed downloads as nondetections.

Reproduce:

```sh
python examples/spherex_moving_screen.py results/spherex-dr1-moving-spectra/most-174p.tbl --detectors D5 --limit 8 --output /tmp/echeclus.json
python -m unittest discover -s tests -p test_spherex_moving_screen.py
```

Three focused tests pass. No browser suite or whole-repository test run is needed. Actual photometry still needs PSF/aperture and pixel-solid-angle corrections, spectral response integration, contaminant rejection, and empirical systematics. These are screening artifacts, not calibrated spectra.
