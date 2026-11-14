# Published SPOC aperture check

The independent TESS-SPOC FFI reduction covers sectors 12 and 39 of TIC 96100077. Fifteen three-day timing windows pass the same phase and cadence requirements used for the custom pixel extraction. With the native SPOC timestamps and corrected aperture flux, the common eclipse clock has an RMS of 55.56 s about a linear ephemeris. Adding the unshifted Gaia delay template with one free coefficient reduces the RMS to 18.91 s and yields a scale of −0.965. Allowing a different constant per sector gives 50.34 s and 18.55 s respectively, with a scale of −0.937.

Both extractions share the same spacecraft exposures, and this check reuses the same clock and fitting code. It tests aperture and pipeline sensitivity. With two short observing intervals, the sector slopes can constrain a common scale even for some alternative outer phases. These data alone cannot establish a unique 272.72 d orbit, measure its phase independently, determine a mass or establish that a third star is bound.

The products are the original MAST TESS-SPOC high-level FFI light curves. The script keeps zero quality flags, finite TIME and finite PDCSAP_FLUX, and normalizes to each sector median. TIME is native BJD(TDB) − 2,457,000. Matching to the reconstructed target-corrected timestamps by nearest cadence gives median differences of +0.0009 s and +0.0011 s in sectors 12 and 39. The manifest records original FITS SHA-256 and the SHA-256 of the TIME, PDCSAP_FLUX and QUALITY columns in that order, using little-endian float64 or int64 and normalized floating-point NaNs. [MAST TESS-SPOC data collection](https://archive.stsci.edu/hlsp/tess-spoc)

From this directory, use the dependencies of the parent investigation. The default analysis directory is the parent directory containing measure.py.

```sh
python check.py /path/to/fits-cache /path/to/replay
diff timings.csv /path/to/replay/timings.csv
diff summary.csv /path/to/replay/summary.csv
```
