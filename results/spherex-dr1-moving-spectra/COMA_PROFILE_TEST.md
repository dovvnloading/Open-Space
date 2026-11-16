# Conditional steady-coma spatial test

## Hypothesis and limits
For optically thin stationary isotropic constant-speed outflow, with photodestruction negligible over the fitted aperture, a volume density proportional to r^-2 projects to surface brightness proportional to rho^-1. This is the inner-coma limit, not a full Haser or fluorescence retrieval. Dust may also follow this profile. NASA's [PSG atmosphere documentation](https://psg.gsfc.nasa.gov/helpatm.php) and [isocoma implementation](https://github.com/nasapsg/isocoma) describe the more general assumptions and destruction model.

We independently implemented softened projected profiles (rho² + core²)^(-p/2), convolved with each image's embedded ePSF on a five-times sampled 401-square grid. Grid p=0.5,1,1.5 and core=0.1,0.5 native pixels tests sensitivity, not physical parameter estimation. Each image is fitted with point+constant, coma+constant, and nonnegative point+coma+constant at radii 3.5,5,8 pixels. Centers are fixed to archive ephemerides. Formal variance, OUTLIER-inclusive masking, and finite ePSF wings impose the earlier limitations.

## Fixed p=1, core=0.1, radius=5 results
| Epoch / approximate µm | Pixels | Point χ² | Coma χ² | Mixture χ² |
|---|---:|---:|---:|---:|
| Dec25 / 4.2155 |79|142.97|65.58|65.54|
| Dec25 / 4.2517 |79|820.04|146.88|146.43|
| Dec25 / 4.2882 |78|580.62|99.89|99.89|
| Jan2 / 4.1677 |77|109.29|101.52|101.52|
| Jan2 / 4.2037 |79|108.57|94.79|94.79|
| Jan2 / 4.2398 |78|477.94|132.25|131.96|
| Jan2 / 4.2762 |78|464.37|157.50|147.51|

Two parameters are fitted in point and coma models; three in mixture. These are descriptive goodness-of-fit comparisons, not discovery significances. No multiple-grid-selection correction or correlated/source-dependent covariance has been applied.

## What this changes
The bright channels are poorly described by an unresolved source, and a steady-coma-like profile materially improves them. However the December continuum improves as well. Greater improvement in a brighter feature does not prove gas is more extended than dust. Bright-channel residuals still exceed formal noise, and slope/background/core choice changes component amplitudes. A simple spatial coma fit therefore does NOT repair or validate the previously failing frozen spectral shape.

## Reproduce
```sh
python examples/swan_coma_profile_test.py \
 --input results/spherex-dr1-moving-spectra/swan-full-d5.json \
 --cache /tmp/spherex-cache \
 --output results/spherex-dr1-moving-spectra/coma-profile-test.json
```
Cache entries use SHA256 of exact manifest URL; rerun the moving-screen acquisition to retrieve missing files. Model normalization inside five pixels is arbitrary and not a production rate. Full grid, input file hashes and masks are in JSON.

Next discriminating step is spatially weighted spectral response with an actual fluorescent molecular template plus a dust component, including field neighbors and calibrated systematic covariance. No chemical identification or validated new discovery follows from this test.
