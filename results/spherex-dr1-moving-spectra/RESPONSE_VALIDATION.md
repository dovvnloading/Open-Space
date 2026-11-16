# SWAN spectral-response validation: provisional forward-model diagnostic

Official calibration retrieved in full from https://irsa.ipac.caltech.edu/data/SPHEREx/docs/spectral_files_20260812/SPHERExSpectralFiles_202607.tar (436705280 bytes). Product guide https://irsa.ipac.caltech.edu/data/SPHEREx/docs/SPHERExSpectralFiles_202607.pdf, draft July22 2026. Hashes in calibration_hashes.json. Extracted Band5 MiddlePixel curves, nominal17-channel geometric map, fiducial summaries and per-pixel BC/BW parameter maps. No top-hat bandpass is used.

## Coordinate and response treatment

L2 sky-coordinate image alternate WCS A locates target in active pixels. Band5 parameter maps explicitly use DETCOORD=DEVICE, unlike channel map DETCOORD=sky. MWIR image is flipped in x to match SWIR, so evaluate DEVICE map at x_device=2039-x_sky, y_device=y_sky. Bilinear interpolation gives BC/BW. Nominal channel is from rounded sky-coordinate pixel. Shift middle-pixel curve by exactBC minus corresponding fiducial BC, per guide section3. No bandwidth broadening with MeanAllPixel curves. Approximate embedded WCS is retained solely as a comparison.

A delta f_lambda source at lambda0 gives relative photon-weighted response lambda0*R(lambda0)/integral[R(lambda)/lambda dlambda] under an ideal flat-Fnu calibration. This is a **shape diagnostic**, not a fluorescence model, production rate, absolute line flux, final color correction or proof of CO2. Actual source has finite spatial extent; this script evaluates the nucleus position and not a coma-weighted aperture bandpass. Near-zero negative transmission samples in delivered curves are retained as measured and must not be interpreted as negative physical response.

## Numerical diagnostic

December25 exact BCs 4.215217,4.251459,4.288034 micron; raw inclusive aperture sums5.70,17.94,13.44.
- Response to delta4.260:1.692,426.268,71.458 (arbitrary common units).
- Response to delta4.267:1.292,307.978,183.529.
Thus the red-side response is much larger in the measured sequence than a4.260 delta alone would predict. This does not refute CO2: molecular band structure, continuum/background and extended-source spectral averaging require proper modeling. It DOES preclude identifying the entire feature by a naive single-wavelength spike.

January2 exact BCs4.239820,4.276048 have delta4.260 responses203.821,288.848, compared to measured sums11.56,11.90. Relative response pattern differs from December25, so retain separate per-epoch models and require a physical band template rather than forcing one delta line with pooled data.

Missing or rejected aperture values remain null, not zero; notably December22 blue-band image is available but its measured aperture is invalid under the engineering extractor and must not be silently used.

## Reproduction

Use the portable CLI commands in EMPIRICAL_FIT_REPORT.md, with the complete frozen swan-full-d5.json input. Preserve input URLs/hashes and calibration_hashes.json. Do not commit the436MB tar.

Unresolved: exact radiometric convention/color correction; detector-frame flip independently verified against calibration source code; per-aperture bandpass averaging; full gas fluorescence template; absolute solid angles; PSF/extended coma corrections. Exact map centers agree approximately with visualization WCS (typically within few1e-4micron), supporting but not independently proving frame convention.
