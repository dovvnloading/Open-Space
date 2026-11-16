# SWAN descriptive spatial falsification

Independent code, same observations. This is not independent observational confirmation.

## Reproduce

```bash
python swan_broadening_review.py --manifest /tmp/swan-pilot-retry.json --cache /tmp/spherex-cache --output swan_broadening_review.json
```

CLI requires numpy, scipy, astropy. Input manifest holds exact URLs; cache names are SHA256(URL). Output includes FITS hashes. No network needed after caching.

## Method and result

Matched Dec25 radius2.5-native-pixel cores; nearest detector-grid EPSF, constant fitted background, amplitude and centroid ±0.5px. Additional isotropic Gaussian sigma constrained0–0.75nativepixel. EPSF padded24 oversampled pixels on each side before 4sigma Gaussian convolution; original support33pixels with oversampling5. Unknown wings outside original EPSF are not reconstructed. Mask permits SOURCE and OUTLIER only, and rejects any other invalidcore pixel. All these exposures fail the strict OUTLIER gate.

| Approx wavelength µm | Point χ² | Broadened χ² | Gaussian sigma native px | Δχ² |
|---|---:|---:|---:|---:|
|4.21545|24.78|22.15|0.510|2.63|
|4.25171|100.51|26.20|0.638|74.31|
|4.28817|54.31|16.43|0.578|37.88|

There are21,19,19pixels respectively; point and broadened models have4 and5 fitted parameters. The extra shape freedom fits the two brighter channels appreciably better. The continuum's best-fit width is ALSO0.51px, but weakly constrained: do not claim gas is resolved relative to continuum merely because its improvement is greater. Brighter signals constrain deviations from a point model more strongly. This test does not yet compare a common broadening width across channels against distinct widths. Profile grids are in JSON.

No p-value or detection significance follows from these χ² differences: noise correlations, EPSF uncertainty, original stamp truncation, mask selection and background/shape degeneracy are omitted. Gaussian broadening is descriptive, not a physical coma, gas classification or discovery. Fits use same pixels previously screened. A physical trail model, local stars at matched signal level, and shared-width vs separate-width residual comparisons are next falsification tests. Low-level alternate WCS A is used only for EPSF-grid lookup; celestial WCS determines targetpixel. Astropy emits inherited SIP warnings for alternate WCS; these should be independently validated before a publication-grade morphology result.

## AllWISE neighbor scale, not exclusion of contamination

Using W2 zero point171.787Jy and Fν=Fν0×10^(-0.4m), nearest Dec25 catalogsource W2=17.193 at8.916arcsec corresponds to0.0228mJy; nearest Dec22source W2=16.490 at12.028arcsec corresponds to0.0436mJy. Dec22brightsource W2=8.658 at97.45arcsec corresponds to59.1mJy. These are historical broadband catalog flux estimates, NOT SPHEREx narrowband predictions. Primary conversion guidance: https://irsa.ipac.caltech.edu/data/WISE/docs/release/All-Sky/expsup/sec4_4h.html . Catalog input CSVs were supplied by engineering; no independent catalog re-query performed.

Even assuming100% capture, nearest Dec25source supplies only0.0228mJy at its historical W2 flux. Conditional on a verified~10mJy SPHEREx target scale and similar color/no variability, this is~0.23%. That conditional argument cannot exclude a narrow spectralfeature, changed brightness, uncatalogedsource, or catalogue upper-limit/quality issue. Need photometric quality/SNR flags before treating each W2 magnitude as a measurement. No catalog-neighbor subtraction performed. Bright-star far wings cannot be validated using a33/5pixel EPSF stamp: setting out-of-stamp intensityzero is a model limitation, not evidence that physical wings vanish.
