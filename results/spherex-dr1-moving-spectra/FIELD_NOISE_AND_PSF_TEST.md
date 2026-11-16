# Empirical field-noise and compact-source aperture tests

Question: could the held-out January 2 spectral structure be explained by ordinary background noise or compact-source aperture sampling? These tests do not identify a molecule or establish novelty.

## Method
Freeze the December Gaussian center/width. Translate the comet track over a deterministic east/north grid from -120 to +120 arcsec in 20-arcsec steps, retaining radii 60–140 arcsec. Apply identical 2.5-pixel apertures, 8–12-pixel background annuli and the explicit OUTLIER-inclusive mask. Require valid extraction in every exposure; do not reject controls based on flux or catalog membership. Evaluate the true local calibrated wavelength and response at each control's detector position. Fit a linear continuum plus the frozen feature shape and divide the feature amplitude by its formal error. Repeat on September 11.

Controls overlap and are not independent survey trials. Finite ranks are descriptive, not global false-alarm probabilities. Blank-field scatter does not estimate source-dependent calibration, PSF or mask systematics.

## Results
- January 2 target statistic: 14.347. Among 65 valid controls, median=0.0705, robust MAD scale=1.1048, maximum=2.9023. No control exceeds the target; descriptive finite rank=1/66. Do not call this a 14-sigma discovery.
- September 11 target statistic: 0.5227. Among 44 valid controls, median=0.0854, robust scale=1.1971, maximum=5.5455. Sixteen controls exceed the target; rank=17/45. No strong feature is established that epoch.
- January's frozen-feature fit still has chi-square 43.23 for one degree of freedom with a linear baseline. Field-noise controls do not rescue the mismatched spectral shape or justify an assumed 10% source-flux floor.

## Compact-source alternative
Evaluate each exposure's ePSF at the actual detector position and subpixel phase, optionally convolved with Gaussian widths 0.5 or 1 native pixel. Apply the same aperture/background geometry to a unit noiseless source. Zero padding avoids convolution edge loss, but cannot restore unmeasured PSF wings; these are fractions of a truncated model, not calibrated total-flux corrections.

Across January 2, radius-2.5 fractions span 0.9768–0.9793 for a point source, 0.9639–0.9727 with 0.5-pixel broadening, and 0.8577–0.8846 with 1-pixel broadening. These compact-source differences cannot alone generate the measured rise from raw sums 0.95/2.08 to 11.56/11.90. Extended wavelength-dependent coma and calibration effects remain unresolved.

The tested products report pipeline 7.0.5. The documented earlier QR1/QR2 PSF-header error is not a reason to impose a generic 5% floor on these products.

## Reproduction
Scripts: `examples/swan_empirical_field_test.py` and `examples/swan_psf_aperture_control.py`. Both accept input/output and calibration/cache paths through their CLI. Use frozen `swan-full-d5.json`, `empirical_feature_fit.json`, the official hashed calibration products and exact-URL FITS cache. Complete controls/rejections and aperture fractions are preserved in `empirical-field-test.json` and `psf-aperture-control.json`.

## Conclusion
The candidate passes stronger local empirical-noise and compact-source sampling checks. It remains a reproducible target-associated spectral-excess lead, not a validated CO2 abundance, temporal activity change, or important new discovery. A physical band/extended-source model and independent observational corroboration remain necessary.
