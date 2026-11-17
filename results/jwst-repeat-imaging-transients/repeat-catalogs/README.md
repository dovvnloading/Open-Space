# Matched-catalog search in four fields

This screen uses 48 distinct public STScI NIRCam Level 3 source catalogs across COSMOS, Abell 370, WHL0137 and MACS J0138. The dates, products and filters are in [pairs.csv](pairs.csv). The six pairs yield 25 filter comparisons. Both images must contain a cataloged source for this test to see a change. Transients visible in only one epoch, hostless objects and unresolved components on extended hosts require a different selection function.

The initial selector produced 356 single-filter triggers and 35 same-sign, multi-filter groups, including groups with opposite signs at the same position. These are catalog-level triggers under review. Cross-filter agreement alone can also arise from blending, PSF changes or catalog backgrounds. The score uses the larger of the formal uncertainty, the empirical scatter at comparable flux, and 2.5% of mean flux. The threshold is 5 times this scale and 0.30 mag in at least two filters. It is a ranking statistic. Object-level fluxes use the archive aperture-corrected total, with the reference multiplied by the empirical per-filter calibration factor.

Run from this directory using Python, NumPy, pandas and SciPy.

```sh
python fetch_catalogs.py /path/to/catalog-cache
python compare.py /path/to/catalog-cache
```

[catalog_inputs.csv](catalog_inputs.csv) identifies each input by SHA-256. [screen_counts.csv](screen_counts.csv) gives match counts by filter. [coincidences.csv](coincidences.csv) gives trigger positions and [matched_photometry.csv](matched_photometry.csv) records all matched-band fluxes at them. Product start MJDs represent the start of a mosaic's observation interval, which can span several days.
