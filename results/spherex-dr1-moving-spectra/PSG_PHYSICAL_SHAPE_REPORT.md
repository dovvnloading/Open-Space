# C/2025 R2 SWAN: physical CO2 template sensitivity test

## Status
A documented NASA PSG cometary CO2 fluorescence spectrum explains the repeated Band5 excess substantially better than a featureless continuum. This is a **CO2-consistent candidate**, not a confirmed new molecular discovery, production-rate measurement, or independent observing confirmation.

## Physical source and explicit assumptions
Official public API: https://psg.gsfc.nasa.gov/api.php
API documentation: https://psg.gsfc.nasa.gov/helpapi.php
Original official configuration: https://psg.gsfc.nasa.gov/data/templates/Comet_67P_from_Rosetta.txt
Model reference: Villanueva et al. 2018, 2022, as identified in returned PSG headers.

The original 67P configuration explicitly supplies GSFC[CO2] fluorescence and a 500000-second lifetime. Saved configurations replace its atmosphere with pure CO2 and use a generic 2-au heliocentric distance, 1-au observer distance, zero Doppler velocities, 500 m/s expansion, arbitrary 1e25 molecules/s normalization, and 15.5-arcsec beam. These are **shape-sensitivity assumptions**, not SWAN ephemerides or inferred physical properties. Temperatures 20, 60, 100, 200 K are a discrete sensitivity grid, not a temperature retrieval. Original unrelated object fields remain in the saved template, so these configurations must not be interpreted as an actual SWAN physical model.

All four successful responses identify Cometary Emission Model, CO2, and GSFC[CO2]. Only their final **CO2 column** is used; nucleus and dust columns are excluded. Output is F_lambda in W/m²/µm at R=3000 from 4.10 to 4.45 µm. Each spectrum is area-normalized before fitting. The model is integrated through the official Band5 MiddlePixel response shifted to each target's exact calibration-map wavelength. Photon weighting is integral(R F_lambda lambda)/integral(R/lambda), appropriate to an ideal flat-Fnu calibration; independent radiometric color-correction verification remains necessary. Only centroid response is used, not a full spatially weighted aperture convolution.

## Split and results
Dec22 and Dec25 are exploratory training dates, three usable aperture measurements each. For each fixed physical spectrum, each date fits a signed constant continuum plus nonnegative amplitude. The four January 2 measurements were held out for the original December-trained Gaussian test. Its failure on January 2 was known when the CO2 and other model-family comparisons were designed. Here January measures later-epoch shape transfer with amplitude and continuum refitted.

| Fixed CO2 temperature | Dec training chi-square | Jan2 chi-square (2 residual DOF) |
|---|---:|---:|
|20 K|161.684|18.588|
|60 K|52.934|13.003|
|100 K|6.473|4.879|
|200 K|59.276|10.369|

The 100-K grid member is selected using Dec training, not Jan2. Training has two residual DOF **before accounting for temperature-grid selection**; do not turn this into a temperature confidence interval. Jan2 observed aperture sums are [0.9503,2.0768,11.5577,11.9043], versus [1.2290,1.7686,11.3701,12.1005] predicted after its two nuisance parameters are fit. Residuals are [-1.36,+1.45,+0.66,-0.70] formal sigma.

Jan2 constant-only continuum has chi-square1707.64/3 DOF; linear-only continuum249.06/2. The previously tested shared empirical Gaussian gave44.8/2. A physical template is therefore materially more explanatory than those alternatives, but this comparison does not establish uniqueness against every molecule, contamination, or calibration artifact.

At aperture3.5 pixels the same100-K shape gives training28.55 and Jan2 8.35/2: aperture dependence remains. An assumed10% fractional floor gives Jan2 2.47/2 at aperture2.5, but **this floor is not a calibrated uncertainty**. A linear continuum plus feature has three parameters per three-point training epoch and fits exactly; it cannot select temperature from these data. Its small residuals must not be presented as independent evidence.

## Limitations that block discovery claims
- Core OUTLIER flags were deliberately retained diagnostically; a production-quality moving-object extraction is still required.
- Aperture-response wavelength spread is roughly1.2–1.4 nm at radius2.5 pixels, not included in this centroid-response convolution.
- Correlated pixels, sky structure, wavelength-dependent coma morphology, and radiometric uncertainties are not fully modeled.
- No measured SWAN ephemerides, solar spectrum/velocity sensitivity, optical-depth sensitivity, molecular-template competition, or production-rate retrieval are supplied here.
- This uses the same SPHEREx exposures as other analyses, not independent validation.
- Literature search did not locate earlier SWAN CO2 spectroscopy but does not prove priority.

## Reproduction
Retain all four `psg/co2_shape_*K.cfg` and successful `*_response.txt` files and original `psg/comet_template.txt`. Submit requests **sequentially**; simultaneous public API calls produced one empty response and one wrong-atmosphere error, both replaced by successful sequential retries. Always validate response headers.

```sh
for t in 20 60 100 200; do
 curl -A 'Mozilla/5.0' --max-time 150 -sS -d type=rad \
  --data-urlencode file@psg/co2_shape_${t}K.cfg \
  https://psg.gsfc.nasa.gov/api.php -o psg/co2_shape_${t}K_response.txt
done
python fit_psg_feature.py --input psg_fit_input_snapshot.json \
 --output psg_feature_fit.json --calibration-root calibration --model-root psg
```

See `psg_hashes.json` for exact artifact SHA256 and `psg_feature_fit.json` for all fitted predictions, residuals, assumptions and input hash. Calibration download/extraction is described by `RESPONSE_VALIDATION.md` and its accompanying scripts.
