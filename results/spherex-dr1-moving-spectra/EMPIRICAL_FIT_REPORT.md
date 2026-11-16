# SWAN empirical broadened-feature validation

## Verdict
A structured excess is inconsistent with the simplest smooth continua at nominal pixel-noise errors. But a common empirical feature shape fitted to December data does NOT transfer to January2 within those errors. Neither result identifies a molecule, proves an astrophysical anomaly or establishes a discovery.

## Model
Use exact official Band5 shifted middle-pixel responses from `swan_response.py`, which locates the target on full-resolution band-center maps. Convolve a normalized Gaussian in F_lambda with each response and ideal photon-count response weighting. This is a phenomenological shape, NOT a CO2 rotational/fluorescence model. Different epochs have separate amplitude and constant continuum. Training center permitted4.23–4.29um, Gaussian sigma0.0005–0.04um. Positive emission amplitude required. Background continuum may be signed.

Fit December22 and25: three accepted samples per epoch, six measurements total. Two shared shape parameters plus two amplitudes plus two constants = SIX parameters. Training has ZERO residual degrees of freedom. Exact interpolation must not be sold as detection-model validation. Allowing a linear continuum gives eight parameters for six observations: the shape is underidentified. A best-fit optimizer output in that model has no unique physical meaning.

Hold out January2 four-point sequence; freeze December shape but refit January2 amplitude and continuum. This tests shape transfer, not completely parameter-free prediction.

## Results
- 2.5pixel aperture, constant baselines: fittedcenter4.26312995um, sigma0.02146367um; training chi-square essentiallyzero (zeroDOF). Heldout chi-square44.7991 for2DOF. January2 residuals in order of increasing wavelength are -1.10,+0.17,+5.30,-3.93 formal standard deviations. Shape transfer fails.
- 3.5pixel aperture: center4.26499871um,sigma0.02093057um; heldoutchi-square91.6074/2DOF. Failure is not removed by increasing aperture.
- An assumed10% fractional error floor in quadrature changes heldoutchi-square to2.9305/2DOF. This is sensitivity analysis ONLY; it neither measures nor justifies a10% errorfloor.
- Linear-continuum training is underidentified and cannot select a physical Gaussian shape. Its optimizer-selected example has heldoutchi-square44.0029/1DOF.
- January2 continuum-only fits: constantchi-square1707.64/3DOF; linearchi-square249.06/2DOF. A simple smooth continuum is inadequate under the extraction/noise assumptions, but the molecular identity and instrumental alternatives are unresolved.
- Five-pixel sensitivity was not fitted because at least one required aperture value is null. Missing extractions were not replaced with zero or lower-quality estimates.

The nominal six-point training fit has an apparently small residual but cannot establish an accurate width. Even a coarse region with training delta-chi-square<=6 spans center4.262–4.266um and sigma0.00939–0.02618um. This grid region is NOT a calibrated confidence interval.

## Remaining alternatives
Finite molecular-band structure (rather than Gaussian); residual sky structure; aperture/PSF differences across detector positions; wavelength-dependent extended-coma sampling; calibration error; and real temporal changes all remain viable. Source selection and exploratory model choice preclude treating formal chi-square differences as discovery significance. Do not claim unusual CO2/H2O composition, production rate, temporal gas change, or first CO2 detection from this test.

## Portable reproduction
Python dependencies: numpy,scipy,astropy. Calibration guide and exact tar route are in RESPONSE_VALIDATION.md; calibration_hashes.json records SHA256s. Data input snapshots retain source URLs and FITS SHA256s.

```sh
curl -L 'https://irsa.ipac.caltech.edu/data/SPHEREx/docs/spectral_files_20260812/SPHERExSpectralFiles_202607.tar' -o calibration.tar
mkdir -p calibration
tar -xf calibration.tar -C calibration --wildcards 'ChannelMaps/*17*' 'FidSpecSummary/SubChannel17/Band5/*' 'MiddlePixel/SubChannel17/Band5/*' 'ParameterFiles/*Band5*'
python examples/swan_response.py --input results/spherex-dr1-moving-spectra/swan-full-d5.json --output swan_response_full.json --calibration-root calibration --cache-root FITS_CACHE
python examples/fit_empirical_feature.py --input swan_response_full.json --output empirical_feature_fit.json --calibration-root calibration
```
The frozen fit verifies all ten expected exposure IDs. The response replay checks the archived checksum of each acquired FITS and reports acquisition or calibration exclusions by exposure ID. Each FITS_CACHE filename is SHA256(exact URL UTF-8)+'.fits'. Download the URL verbatim and validate bytes against its recorded sha256. Input snapshot may include extra unneeded epochs; the empirical script explicitly selects December22,25 and January2 withbandcenter>4.1um. Archive snapshots before updating partial survey inputs.

Do not commit the436MB tar unnecessarily. The source route, exact hash and selected calibration products preserve provenance.
