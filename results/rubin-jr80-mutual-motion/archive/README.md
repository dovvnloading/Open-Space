# Historical relative astrometry of (556068) 2014 JR80

A search by Solar System Object Image Search (SSOIS) at CADC recovered MegaCam exposures at 13 dates between 2020 and 2022. Fitting the two components on a single detector independently at each visit yielded usable relative astrometry on nine dates. The discovery circular had already reported six ground-based epochs between 2019 and 2023, with separations of 0.5–1.0 arcsec. It also reported four 350-s HST exposures on 2026 June 13, separation 0.846 arcsec, position angle 16°, secondary 0.13 ± 0.01 mag fainter. The circular does not identify the six ground-based dates. The contribution here consists of measured vectors and their evolution. Credit for identifying the system and the archival binary detections belongs to Proudfoot et al. [CBET 5739, 2026 September 14](http://www.cbat.eps.harvard.edu/iau/cbet/005700/CBET005739.txt).

## Historical measurements

| MJD range | Result |
|---|---|
| 59027–59028 | Three separate frames per night place the axis at 173–175° and the centers 0.958–0.987 arcsec apart. |
| 59052 | Three frames give 0.582–0.600 arcsec at 40.3–43.2°. |
| 59373–59379 | Five nights resolve approximately 0.49–0.54 arcsec. The axis advances from 98.7° to 124.2°. Each night has one selected exposure. |
| 59403 | Three different filters give 0.37–0.39 arcsec, with a position-angle spread of 66.5–79.6°. |
| 59053–59054 | Interference from a bright unrelated source, an extreme residual, or loss of detector coverage on multiple frames. Withheld. |
| 59371 | The two-source fit and excess second moment agree on an axis near 72°. This single exposure gives separation 0.77 arcsec, discrepant with the nearby cluster of measurements. Withheld from the primary orbit fit. Including the exposure at the nominal 0.04″ floor raises the 103-day fit χ² from 2.26 to 54.46. [Two-family refit](../orbit/excluded_cfht.csv) |
| 59725 | Three frames have peak signal ≈2.2–2.4 times the local pixel noise, and none constrains an axis. |

Angles are measured east of north modulo 180°. The two components have nearly equal brightness, so the two possible assignments differ by 180°. The bright 2020 sequence obtains JPL topocentric sky displacements of about 0.16–0.18 arcsec over a 205-s exposure, substantially smaller than the resolved June separations. Seeing, relative motion, differential refraction, local blending and a non-Gaussian PSF affect the astrometry. Comparison stars experience essentially zero motion. The orbit fit must allow systematic uncertainty larger than a formal centroiding error.

## Reduction and checks

Each filter image is fitted separately. The search uses the SSOIS position and extracts a 17×17-pixel box centered on its nearest pixel. Up to 28 bright nearby local peaks define the stellar PSF; unresolved, approximately round stars supply individual elliptical Moffat fits. The robust median stellar shape is used for both single- and double-source models. The latter permits the centroid, background, amplitude, displacement and flux ratio (0.65–1.0) to vary. Sky-noise and positive Poisson variance weights are shared by both models. No angle or separation forecast enters a source fit. The reported fit improvement is the difference in weighted residual sums, without conversion into a significance. Comparison stars within the PSF-building population give a check on shape mismatch. Their gains are not statistically independent null trials.

`moments.py` independently estimates position angle directly from each pixel stamp by iterated Gaussian-weighted second moments, then subtracts the median field-star covariance in sky coordinates. It uses neither the Moffat fits nor their angles. Residual eigenvectors agree within roughly 5° on the accepted exposures, including the 2021 rotation sequence. A negative or nonfinite moment result is retained; the moment code does not assign component photometry or separation. This check uses the same pixels.

`inputs.csv` gives the ephemeris position, epoch, exposure length, source URL and checksum for each product. `measurements.csv` and `comparison-stars.csv` preserve fit quality and exclusions. The query URLs are in `coverage.json`. CADC SSOIS positions refer to predicted footprints; a listed image can fall in a detector gap. The field-search and acceptance rules do not provide a completeness estimate.

To replay an exposure, install NumPy, SciPy and Astropy, then run:

```sh
python replay.py 2497047p 2497069p --cache /tmp/jr80-data
python moments.py /tmp/jr80-data/2497047p.fits /tmp/jr80-data/2497069p.fits
```

A CFHT request downloads the small published cutout. A DECam request downloads the full multichip FITS file (approximately 300 MB), checks its archive MD5, transforms the predicted position through every chip WCS, and extracts the surrounding field when the target is at least 12 pixels from the edge. This geometric screen precedes photometric inspection. Checksums for CADC cutouts describe those downloaded at the time of analysis; a subsequent pipeline release or change to the server's FITS packing can change the bytes.

## Blanco holdouts

Three early exposures with long integration times were rejected by the detector WCS screen. On 2022 August 29 (MJD 59820.997) the ephemeris position was 66 pixels beyond the nearest chip edge. On August 30 it was only 8.4 pixels inside. The first 2023 May 19 pointing (MJD 60083.3534) missed by 1,196 pixels. At MJD 60423.3675 the 2024 April 23 exposure missed by 59 pixels. These epochs supply no fit.

Both 60-s z-band images on 2024 April 25 had the target on-chip. The peak/sky-noise ratios were 3.28 and 5.46, and the binary vs single improvements were 9.35 and 3.45. The 35-s r-band 2025 March 8 twilight exposure gave 3.98 and 3.39, respectively. None of these three frames supports an independent component separation. Their formal best vectors remain in the table for reproducibility and have no orbit-fitting weight.

After the first three weak frames, additional high-coverage pointings with longer integrations were selected using primary FITS headers, before examining their target pixels. For the six added on-chip frames the assessor recorded this quality screen before seeing the first five remaining fits: peak/sky pixel noise ≥10, improvement ≥25, reduced weighted residual sum <2, and improvement larger than all six comparison stars. These empirical thresholds were chosen after the first successful 2022 fit and before viewing orbital forecasts. They are not a calibrated false-alarm probability. Forecasts were preserved before the image measurements were shared; the assessor received no predicted PA or separation.

| UT date | Detector/filter | Separation | Axis PA | Moment PA | Target Δχ² | Largest comparison Δχ² |
|---|---|---:|---:|---:|---:|---:|
| 2019 June 2 | DECam g, 90 s | 0.654″ | 30.92° | 30.71° | 32.85 | 6.58 |
| 2022 August 6 | DECam i, 90 s | 0.660″ | 147.85° | 148.20° | 40.40 | 3.21 |
| 2023 May 19, second pointing | DECam r, 175 s | 0.660″ | 42.09° | 39.04° | 112.44 | 8.91 |

The 2022 flux ratio reached its fitting bound, so the useful measurement is the axis. A separate June 8, 2019 DECam exposure gave Δχ² 7.72, below a local star's gain of 69.54. On July 27, 2021 and July 28, 2023, peak/sky pixel noise was 4.01 and 7.54. These three fitted images join the April 2024 and March 2025 data in the insufficient-evidence group. Six of the nine on-chip holdouts therefore supply no independent binary axis, and all four geometry failures remain in the audit. Any coincidence of the published six epochs with these three detections would credit the archival binary identification to the original team. The individual quantitative relative-astrometry products and the withheld-vector comparison can be examined here.
