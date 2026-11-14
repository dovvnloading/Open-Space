# A pre-TESS timing test with KELT-South

KELT observed TIC 96100077 on 665 nights from May 2013 to October 2018. Before collecting the four light curves we recorded the comparison against the previously published Gaia orbit and the signed amplitude from our TESS timing analysis. The archive provides 6,028 detrended exposures (TFA) and 6,045 raw exposures. The two reductions share observations. We use the TFA set for the primary test. [Test design](protocol.md) · [KELT archive](https://exoplanetarchive.ipac.caltech.edu/docs/KELT.html) · [Source metadata](catalog.csv) · [Exact products and digests](sources.json)

The Gaia radial displacement predicts 443.7 light-seconds peak to peak for its RV tracer. A signed timing scale of −0.80164 predicts a 355.7-second peak-to-peak shift of the eclipses, in the opposite direction. At the individual exposure level this is a small phase-dependent change to an otherwise large and repeatable modulation.

| Quadratic inner ephemeris | Detrended | Raw |
| --- | ---: | ---: |
| Held-out RMS, zero outer delay (parts per thousand of reference flux) | 17.806 | 22.008 |
| Held-out RMS, fixed signed delay | 17.700 | 21.927 |
| Held-out RMS, reversed sign | 18.087 | 22.254 |
| Reduction in held-out mean squared error (ppt²) | 3.75 | 3.56 |
| Paired 21-day resampling interval for reduction (95%, ppt²) | 0.81 to 6.98 | −0.02 to 7.42 |
| Free scale, full sample, 21-day cluster interval (95%) | −0.851 [−1.220, −0.483] | −0.769 [−1.149, −0.388] |
| Other phase shifts as good as nominal, of 99 | 6 | 1 |
| Omitted years with improved predictions, of 6 | 5 | 4 |

Each night is withheld once in five interleaved folds. We fit a four-harmonic light curve per orientation and a common linear or quadratic phase correction to the other nights, then predict the withheld fluxes. The amplitude, eccentric outer orbit, phase and sign in the fixed test come from the earlier investigation. Relative flux uses a fixed archive magnitude per orientation. Each exposure gets equal weight. The cluster interval for free scale linearizes the local fit. Intervals on paired errors resample nights or 21-day blocks conditional on the fitted predictions. The largest interleaved-night gains occur in 2016 and 2017. Resampling six whole years gives a wider conditional interval that includes zero. Whole-season exclusions provide a separate test of extrapolation. Alternate outer phases overlap strongly and their rank has no direct false-alarm interpretation. Daily and sidereal frequency aliases fit the inner eclipses worse. [Both ephemerides](fits.json) · [Phase controls](phase-controls.csv) · [Omitted seasons](review-omitted-seasons.csv) · [Inner aliases](review-aliases.csv)

## Source and timestamp verification

NASA's source-search manifest for the KELT field S36 record at 264.52645°, −49.34803° pairs eastern file number 040595 with western number **038138**. The western file numbered 040595 belongs to another source. The analysis checks the SHA-256, exposure count, minimum, maximum and mean magnitude and time span against the target metadata for each product. The TFA ASCII timestamps include an extra constant of 2,450,000 days. Subtracting it matches every retained TFA timestamp to a raw exposure within 0.000002 d. The primary analysis uses those archive-reconciled timestamps. [Archive field definitions](https://exoplanetarchive.ipac.caltech.edu/docs/API_kelt_columns.html) · [Original download manifest](download-manifest.txt)

An annual timing component with semiamplitude 473 s in TFA and 482 s in raw is almost exactly opposite the Earth's projected barycentric displacement (correlations −0.991 and −1.000). Adding one target-direction geocentric light-time correction leaves annual semiamplitudes of 46 s and 40 s. We retain this as a retrospective sensitivity test because the archive labels its times BJD and its exact pipeline convention for S36 remains unresolved. A constant clock offset is absorbed by the inner waveform. Under the alternative time convention the detrended free scale is −0.959, the raw scale is −0.903, and the fixed signal improves five of six omitted years in each reduction. The differing absolute scales quantify one real timing systematic. [Clock audit](review-clock.csv) · [Prediction errors with both clocks](review-calibration.csv)

From this directory, Python with Astropy, NumPy, pandas and SciPy can download the four original ASCII files and reproduce all results.

```sh
OPENBLAS_NUM_THREADS=1 python kelt.py --fetch
OPENBLAS_NUM_THREADS=1 python review.py --raw-dir raw
```

The data are from KELT-South and the NASA Exoplanet Archive at Caltech/IPAC. [KELT Survey Light Curves](https://doi.org/10.26133/NEA8) The measurements and interpretation use Open-Space revision 062d4c61aa5b41608a61cf6aed9f82331b548f30 for the TESS and Gaia prediction.
