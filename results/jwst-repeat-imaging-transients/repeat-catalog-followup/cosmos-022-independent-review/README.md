# COSMOS-022: independent native-aperture diagnostic

**Noteworthy follow-up lead:** the measured near-infrared brightening remains positive under these post-inspection changes to the pixel calculation. This describes a candidate for follow-up, not a classification or a claim of a new transient. This analysis was designed after the candidate and original results were inspected; it is not a blind test.

The 34 CAL exposures, coordinates and 204 byte-range hashes are the frozen inputs of the [original validation](../cosmos-022-validation/README.md), which retains the observing context and mask discussion. The ranges total 41,814,016 bytes. The review independently computes aperture weights on a 5×5 subpixel grid, excludes every flagged DQ pixel, takes a median local sky at 0.40–0.58 arcsec, and multiplies the summed MJy/sr by FITS PIXAR_SR and 10¹⁵ to obtain nJy. It calculates radii 0.12, 0.16, 0.24 and 0.32 arcsec, at the reference coordinate and common ±0.02 arcsec east and north shifts. The same sky annulus is used at all radii. No missing-area correction is applied. The resulting [680 rows](diagnostics.csv) retain weighted clean coverage and sky-pixel counts.

For rows with at least 90% clean aperture coverage, centered early-to-late median changes at those four radii are respectively +206.6, +213.9, +236.6 and +202.9 nJy in F162M and +178.2, +192.3, +205.5 and +201.9 nJy in F182M. At 0.16 arcsec the four common shifts give changes of 202.9–226.9 nJy and 185.6–211.5 nJy. Coverage changes the set of contributing exposures by radius and shift, so these are not fixed-cohort growth curves. The common shifts do not test differential registration between visits; the calculation does not independently calibrate PSF changes, detector gains or host gradients. The measurements alone do not distinguish an AGN fluctuation from an unresolved transient.

Use Python with NumPy and Astropy. Supply separate writable cache and output directories on D:. `--fetch` obtains missing ranges from the frozen public URLs; without it replay is offline and missing bytes fail. Every existing or downloaded range is checked against the fixed SHA-256 and size before analysis; a mismatch fails rather than changing the receipt. The output path must be outside this package and cache. For example, from the repository root:

```sh
python results/jwst-repeat-imaging-transients/repeat-catalog-followup/cosmos-022-independent-review/diagnostic.py --cache D:/Open-Space/campaigns/cosmos022-review-ranges --out D:/Open-Space/campaigns/cosmos022-review-output --fetch
```

The recorded output was reproduced byte-for-byte on the test host using Python, NumPy and Astropy. Platform-dependent floating-point differences can affect the final digits. No original data products or results were modified.
