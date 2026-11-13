# Testing a fixed-wavelength artifact against the Ca II K interpretation

Three absorption features in different galaxies were paired after inspecting the pilot output. All three spectra come from SDSS DR17, plate 6497, MJD 56329. A common sky or instrumental line predicts approximately the same observed wavelength on this plate. A transition in each galaxy predicts wavelengths proportional to its own `1 + z`. The DR17 first-pass table contains emission-line wavelengths. A separate SDSS line list tabulates Ca II K at 3934.777 Å and H at 3969.588 Å in vacuum. [SDSS absorption-line reference](https://classic.sdss.org/dr7/algorithms/linestable.php)

| Fiber | Catalog redshift | Observed centroid (Å) | Source-frame centroid (Å) | Apparent K offset (km/s) | Flagged OR pixels in core / full window |
| --- | ---: | ---: | ---: | ---: | --- |
| 0088 | 0.04228695 | 4100.179 | 3933.829 | −72.2 | 0/3 and 2/25 |
| 0096 | 0.03487289 | 4071.925 | 3934.710 | −5.1 | 2/3 and 2/25 |
| 0100 | 0.07697336 | 4238.359 | 3935.435 | +50.1 | 0/3 and 8/25 |

The observed centroids span 166.434 Å. Their source-frame span is 1.606 Å. The unweighted RMS observed residual around their best common wavelength (the mean) is 72.719 Å. Using 3934.777 Å times the catalog `1 + z` as the prediction for each galaxy gives 0.703 Å RMS. Formal redshift errors correspond to 1.6–2.4 km/s, but each recorded centroid uses a three-pixel sampling estimator. The largest residual is about one pixel. A velocity offset or probability of association would require centroid uncertainties and a calibrated line-profile model. The pipeline estimated systemic redshift from these spectra, so this comparison also shares information with the original data.

As a location control, inspect each galaxy at the other two fibers' observed centroids. All six comparison centers have clean 25-pixel windows under the default AND-mask rule. Within ±2 Å, one position also has an above-threshold absorption, in fiber 0088 at 1.12 Å below the fiber 0096 centroid. The other five have no such detection. The broad wavelength shift weighs against one fixed plate or terrestrial line explaining all three peaks. Calibration errors at different wavelengths and selection from many possible features remain plausible contributors. The comparison uses three selected features and supplies descriptive residuals only.

## Sensitivity to individual exposure flags

Every nominal K window contains a nonzero OR-mask pixel. Setting every flagged pixel to unusable breaks all three detector windows. The strict rerun has no coverage of their original centers. To inspect the surviving flux, apply a second, fixed-position index with no requirement for 25 consecutive samples. In each rest frame:

1. Take two sidebands 8–15 Å from the reference wavelength. Within each sideband take the median flux and median wavelength of usable samples. Require at least two pixels in each sideband.
2. Interpolate a straight line between the median anchors. Within ±3 Å of the reference wavelength, calculate `1 − sum(ivar × flux) / sum(ivar × continuum)`. Keep positive inverse variance throughout.
3. Recalculate after requiring every retained pixel to have zero OR mask. Positions and sidebands stay fixed for both reductions.

| Fiber | K flux deficit, default | K flux deficit, strict OR=0 | Strict K pixel counts, blue / core / red |
| --- | ---: | ---: | --- |
| 0088 | 57.7% | 57.7% | 8 / 7 / 8 |
| 0096 | 39.4% | 36.4% | 8 / 5 / 8 |
| 0100 | 52.4% | Insufficient coverage | 4 / 6 / 1 |

The dip persists in unflagged samples of two galaxies. The third lacks a usable red anchor under this stricter rule. The fixed-window index depends on the continuum model and is descriptive. Full individual exposures, additional nights and independent instruments would test the remaining concerns. [SDSS DR17 pixel quality](https://www.sdss4.org/dr17/spectro/quality/)

## The companion test

The baseline screen returns troughs at 3972.029 Å (0088) and 3972.760 Å (0100). Both centroids lie closer to H epsilon at 3971.195 Å than to Ca II H at 3969.588 Å. Fiber 0096 has no above-threshold feature within 4 Å of Ca II H. The fixed index at Ca II H gives positive absorption in all three galaxies, but its ±3 Å window includes hydrogen H epsilon. This is a blended region. Establishing a Ca-specific companion requires a joint stellar/line fit or more informative spectra. [MaNGA DR17 hydrogen reference wavelengths](https://www.sdss4.org/dr17/manga/manga-analysis-pipeline/) · [Published Ca II H and K absorption measurements](https://doi.org/10.1093/mnras/stv1277)

**Disposition.** An ordinary Ca II K feature at roughly the host systemic redshift is a plausible explanation worth recording for future agent triage. The present data leave its exact line formation and the Ca II H contribution unsettled. The three fibers are different galaxies observed with one telescope on one plate and night. The review used no independent instrument or repeat exposure.
