# Independent measurement of the Hubble reference epoch

The known binary (556068) 2014 JR80 was resolved in four consecutive 350-second HST/WFC3 UVIS F350LP exposures on June 13, 2026. The primary-to-secondary displacement averages **+0.233″ east, +0.813″ north**, corresponding to **0.8461″ separation at 15.979° east of north**. The four exposure measurements have sample standard deviations of 0.0006″ and 0.133°. They are repeat measurements within one visit and share camera calibration. These scatter values describe repeatability only.

The original discovery telegram reports 0.846″, PA 16°, and a secondary 0.13 ± 0.01 mag fainter on the same visit. The independent reduction agrees with all three published quantities. The HST vector anchors the Rubin comparison. CBET 5739 also records six earlier binary epochs from 2019 to 2023 with 0.5–1.0″ separations. [Discovery telegram](https://www.cbat.eps.harvard.edu/iau/cbet/005700/CBET005739.txt)

## Data and measurement

The MAST HST target query for `2014JR80` returned four exposures from proposal 18010, all June 13. A second query for proposal 18010 returned the same four records under that name. The target-name query response is preserved in `mast_target_query.json`. Other designations or programs with different target names could exist. [MAST search](https://mast.stsci.edu/search/ui/) · [HST proposal 18010](https://www.stsci.edu/hst-program-info/program/?program=18010)

Each cutout contains a 73 × 73-pixel section of the charge-transfer-efficiency-corrected drizzled DRC SCI extension, in electrons per second. Coordinates in the manifest use zero-based pixel indices. The resolved pair repeats around the commanded moving-target position in all four frames. Larger apertures include nearby contamination in individual exposures. The primary measurement fits an axis-aligned elliptical Gaussian and local constant to each source core within 3.5 pixels of its seed. The independent centroid uses iterated first moments in a 3.5-pixel aperture, subtracting the median background beyond 7 pixels. The secondary is consistently fainter. The mean aperture magnitude difference is 0.131 mag with exposure-to-exposure scatter of 0.016 mag.

The DRC science header's CD matrix converts the centroid displacement to tangent-plane east and north; projection effects over the subarcsecond pair are negligible at the precision quoted. Mean aperture results are 0.8459″ at 15.945°. Two centroid estimates from the same image share photons. The HST exposures and the Rubin exposures are independent observations.

## Reproduce

Run `python measure.py` in this directory with NumPy and SciPy. It verifies each science cutout's pixel digest and emits `measurements.csv`. Run `python measure.py --archive` with Astropy to download every original DRC file, verify its SHA-256, crop the SCI array using the recorded bounds and remeasure. This transfers about 870 MB. Raw-image URLs and both sets of hashes are in `manifest.json`.
