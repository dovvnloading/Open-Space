# Independent exposure checks

The high-priority native tests covered 23 source–wavelength pairs in 21 distinct sources. Of 447 attempted aperture measurements, 239 passed the native quality checks. Thirteen pairs had the core and both sidebands in the independent QR3 release. One is the 3.3 µm aromatic-band candidate in UGC 8856 NED01 described below. The other twelve have absolute native contrast scores below 0.65. The remaining ten pairs lack a clean three-bin QR3 comparison. The supplemental raw-trigger table is reported separately.

The native-image measurements test whether the source aperture has a local spectral contrast in separate exposures. DR1 contains constituent observations of the discovery cubes. QR3 images were taken after DR1 and give the independent test. The two splits of each release are formed at the median MJD of the selected images. A missing contrast means that one or both adjacent samples, or the central sample, lacked usable data.

## Sample and method

The two `targets.csv` files identify the selected sources, wavelengths and reasons. The sample includes the first five full-field template ranks, every restricted-field template score ≥5, both leading spectral-shape outliers and every source whose maximum absolute leave-one-out channel-adjusted score reached 5. A source has two tests if two rankings select different wavelengths. For the shape ranking we test its highest-scoring local residual. This tests one aspect of an unusual continuum.

IRSA TAP returned the detector footprints intersecting each 0.5-degree cone. The queries include the deep and wide observations in the release directories. The index receipts give the precise query, row count and raw response SHA-256. We used the four ordered corners to estimate detector coordinates and an exposure's WCS-WAVE table to predict wavelengths before downloading any target fluxes. The initial high-priority check samples at most six distinct UTC MJD days per bin and release. Choose the closest predicted wavelength on each day, then space the selected days by rank over the available sequence. For the primary wide-field list we subsequently tried other available images, in order of predicted wavelength proximity, when a whole release/bin still had zero usable measurements. This stops at the first usable image and does not select by flux. Such attempts have selection_reason=missing_bin_alternate. We used one detector per feature. DR1 and QR3 detector IDs are the same at a fixed wavelength.

The central sample is within ±0.55 nominal cube-channel widths of the target wavelength. Two sidebands extend 1.5–3 channel widths on either side. Each cutout recomputes sky-to-detector position using the distorted celestial WCS. We bilinearly interpolate the official native-resolution CWAVE and CBAND arrays at that position and recheck bin membership. The input ETags and hashes for these maps are in `science_wcs_D*.json`; for each detector IRSA distributes the same map in both releases. Each exposure CSV includes both the approximate WCS-WAVE values and the science-map values, observation date, native detector position, cutout URL and SHA-256.

A fixed 2.25-native-pixel radius aperture is integrated with 5 × 5 subpixel sampling. We subtract the exposure ZODI model and the residual local median in the 4–7 pixel annulus. Every contributing aperture pixel must have finite data, positive variance, and no flag except SOURCE. The sky annulus needs at least 28 valid pixels. The aperture flag OR and good area fraction remain recorded for rejected data. We repeat the photometry at eight positions nine pixels from the source. At least three off-source measurements must pass. The single-exposure noise is the larger of the propagated variance including the local sky and the 1.4826 × MAD of these blank apertures. Reported flux is the sum of surface brightness × fractional native-pixel area (MJy sr⁻¹ pixel units); absolute point-source calibration and encircled-energy correction were not applied.

The reported contrast is the median central flux minus the linearly interpolated blue and red medians, evaluated at the median measured wavelengths of all three bins. For each bin the error estimate uses the larger of the inter-date MAD scale and median exposure noise, times 1.253/√N, with a 3% relative floor added in quadrature. We propagate these three errors to the contrast. Signed scores are local consistency diagnostics. They omit correlated systematics and small-sample uncertainty in the estimated scatter. Bins can be observed on different dates; intrinsic variability and wavelength-dependent blending can also affect the contrast. A test spanning a few channel widths cannot establish or exclude other spectral shapes. The same detector and local optics are used in both releases. Flag masks and offset apertures give limited protection against unmodelled ghosts, residual persistence, PSF variation and distant bright stars.

## Reproduce

From the study directory, run the following with Python 3, NumPy, SciPy and Astropy. The cache path is disposable. Existing receipt CSVs support a replay against exactly the selected exposures, independent of subsequent QR3 archive growth.

```sh
python visits/choose_targets.py
python visits/science_wcs.py /tmp/spherex-native-wcs visits 1 2 3 4 5 6
python visits/index_visits.py /tmp/spherex-nep-index visits/nep 270 66.56 --detectors 2 5
python visits/index_visits.py /tmp/spherex-wide-index visits/wide 209.53125 30 --detectors 1 2 3 4 5 6
python visits/measure.py visits/nep/targets.csv sources.csv /tmp/spherex-nep-index visits/nep visits/nep /tmp/spherex-native-wcs --replay-dir visits/nep
python visits/measure.py visits/wide/targets.csv wide/sources.csv /tmp/spherex-wide-index visits/wide visits/wide /tmp/spherex-native-wcs --replay-dir visits/wide
python visits/summarize.py
```

Dropping `--replay-dir` performs a new selection from the available footprints. For a new selection run `python visits/complete_bins.py visits/wide/targets.csv wide/sources.csv /tmp/spherex-wide-index visits/wide visits/wide /tmp/spherex-native-wcs` to reproduce the wide-field missing-bin rule before summarizing. A live TAP query can grow as IRSA releases additional exposures. No source photometry is used to pick an exposure within a band or date.

The active-detector coordinate and local pixel-area audits use the selected cutouts. The wing fits compare two unsaturated annuli for the flagged stars and a control. After the supplemental target replay, run the following from the study directory.

```sh
python visits/audit_area.py visits/nep /tmp/spherex-nep-index/cutouts
python visits/audit_area.py visits/wide /tmp/spherex-wide-index/cutouts
python visits/audit_area.py visits/wide_remaining /tmp/spherex-wide-index/cutouts
python visits/wing_check.py visits/nep /tmp/spherex-nep-index/cutouts
python visits/wing_check.py visits/wide /tmp/spherex-wide-index/cutouts
```

Archive definitions and exposure access are documented in the [SPHEREx DR1 release overview](https://irsa.ipac.caltech.edu/data/SPHEREx/docs/overview_dr1.html), [Explanatory Supplement §§2.5, 3.2.5 and 4.3.1](https://irsa.ipac.caltech.edu/data/SPHEREx/docs/SPHEREx_Expsupp.pdf), and [IRSA's native-image cutout tutorial](https://caltech-ipac.github.io/irsa-tutorials/spherex-cutouts/).

## Resolved example in UGC 8856 NED01

All 55 available native D4 footprints containing J135607.93+300452.8 were extracted with the same masks. The published target list selected the source for its negative residual relative to a galaxy template. The data show an **emission** feature at this wavelength in both releases. Fitting each native scan separately in 3.07–3.73 µm gives observed centers of 3.390 and 3.394 µm (DR1, QR3). At the catalogued redshift 0.0310622 they correspond to 3.288 and 3.292 µm. The core centroids of the selected clean native images are within 0.25 arcsec of the catalog position. In each fit the continuum is linear and the intrinsic band is a Gaussian of fixed rest FWHM 0.04 µm convolved with the per-exposure Gaussian CBAND approximation. The numerical fit precision exceeds the spectral resolving power. Detection in both releases supports a real celestial feature compatible with the known aromatic band. Neighboring-source light and non-Gaussian instrumental response remain possible contributions to the profile. Comparable 3.3 µm aromatic bands are observed across star-forming galaxies. [AKARI galaxy measurements](https://academic.oup.com/pasj/article/78/2/454/8441758)

```sh
python visits/profile.py /tmp/spherex-wide-index /tmp/spherex-native-wcs visits/profile_J135607.csv
python visits/profile_fit.py
python visits/profile_position.py /tmp/spherex-wide-index
```
