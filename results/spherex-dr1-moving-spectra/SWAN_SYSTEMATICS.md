# SWAN image-domain systematics checks

Exploratory checks on the same nine exposures, not independent discovery confirmation.

The December25 excess persists when background annuli change from6–10 to8–12 to12–17pixels. At approximate4.252microns, 2.5px raw sums span17.53–18.50; at4.288 they span12.92–14.03; blueward4.215 spans5.11–6.18. Positive-flux centroids within3pixels of the predicted position are offset by less than0.22pixel for these three images. These centroids are not astrometric PSF fits, and their selection region is target-centered.

A deterministic grid of local off-target apertures samples potential contaminants. The valid December25 controls number34–36, with maximum raw sums2.22–4.01. December22 is less clean: off-target sources reach36–62 raw units in two images, stronger than the target. A local-control rank is neither a survey false-alarm probability nor evidence of a rare discovery. These are real crowded-field contaminants requiring spatial modeling.

`examples/spherex_swan_systematics.py` records all aperture radii2,2.5,3.5,5pixels and all three annuli, descriptive centroids, and every local grid control. It uses the explicit OUTLIER-inclusive diagnostic mask and continues to exclude other flagged pixels. No bad core is repaired or silently interpolated. October's centroid offsets exceed1pixel and nearby bright emission reinforces excluding it from a clean season-to-season comparison.

## Important wavelength limitation
The official Explanatory Supplement v2.1 §4.3.1 states embedded WCS-WAVE central wavelengths/bandwidths are compact visualization approximations. Science analysis requires separately delivered per-pixel CWAVE/CBAND and spectral-response information. All current wavelength labels are approximate screening coordinates, not precision-calibrated spectral coordinates. Full response validation is pending. Sources: https://irsa.ipac.caltech.edu/data/SPHEREx/docs/SPHEREx_Expsupp.pdf and https://irsa.ipac.caltech.edu/data/SPHEREx/docs/spherex_spectral_calibrations.html .

Reproduce from the cached exact-URL cutouts:

```sh
python examples/spherex_swan_systematics.py results/spherex-dr1-moving-spectra/pilot-swan-retry.json /tmp/swan-systematics.json
```
