# Rubin sample, image fits and controls

## How the target entered the sample

On 2026-11-16 we searched the Rubin moving-object stream for persistent excess flux outside a point-spread function (PSF), initially as a test for dust around main-belt asteroids. The search used ALeRCE’s Rubin-specific detection table and survey ID 2. The frozen screen required at least three detections with S/N > 25, pipeline extendedness > 0.8 and aperture flux > 1.25 times the PSF flux, with the dipole, PSF, aperture, saturated-center, bad-pixel and shape flags clear. The grouped query searched the live database and returned 12 objects, below its 2,000-row ceiling. This exploratory cut does not measure completeness for dust or companions.

The broker’s MPC cross-match places ten of these objects at perihelion distances below 1.3 AU. The main-belt asteroid 2016 FL8 met the morphology criteria seven times on a single observing night. The remaining object, (556068) 2014 JR80, met them four times over three observing nights at a heliocentric distance near 40.36 AU. Its binary status was already known from Hubble. The follow-up question was whether Rubin measures the pair’s relative motion. We retrieved all five associated ALeRCE detections, including the detection that missed the morphology cut. The selection and known-object match preclude inferring a population-level false-alarm rate or occurrence estimate.

The four ADQL files preserve the exact queries. Each has its frozen CSV response. ALeRCE’s lookup tables map sid=2 to LSST SS Object, band=3 to i and band=4 to z. Calendar dates in the astrometry table use the integer day of the detection MJD. These are the UTC calendar dates of the observations; the Rubin visit prefixes instead identify the local observing nights. The study does not require sub-minute absolute timing. The 2026 Nov 16 snapshot may change as alerts or orbit associations are reprocessed. The broker’s lsst_ss_object summary table returned zero rows on a five-row probe. We used populated detection, orbit and geometry tables throughout.

## Pixel measurements

Using DIA source IDs from ALeRCE, we downloaded FITS images from the Fink LSST cutout API. The original Rubin packet supplies each image, its two-dimensional variance and the normalized local PSF. Different brokers host copies of the same Rubin exposures. We include 15 target cutouts and 40 science cutouts for control SSOs, losslessly compressed. The stamp manifest gives the direct URLs and SHA-256 digests of the original FITS bytes. The primary and variance planes have units nJy and nJy² respectively.

The image model includes a constant sky and a centroid, using either one or two copies of the measured image-specific PSF. The pair adds a two-dimensional separation and flux fraction. Cubic splines interpolate the PSF. The fit divides pixel residuals by the square root of their variance and uses finite positive variances only. The nominal radius is 8.5 pixels (~1.7 arcsec). The centroid is limited to ±4 pixels around the packet's WCS reference pixel, offsets to ±6 pixels per coordinate, and the component fractions to 0.1–0.9. All target fits fall far from the fraction bounds. Five starting configurations include the coincident-pair solution. Sensitivity fits use a six-pixel aperture or broaden the local PSF by convolution with a 0.4-pixel Gaussian.

The FITS linear tangent-plane WCS converts the offset to east and north. We report position angle modulo 180° because the comparable fluxes make component labels ambiguous. Tabulated vectors point east by convention. Parameter covariance comes from the weighted Jacobian around the solution, scaled upward when residual variance exceeds the diagonal noise model. These conditional errors omit correlated pixels and optical or night-shared systematics. The improvement in χ² over one PSF is a diagnostic. At zero separation the extra parameters are nonregular, and blind selection introduces a trials factor.

## Empirical checks

- The five science exposures yield Δχ² of 181–904 in favor of a pair. The faint-component fractions are 0.447–0.497. The object’s centroid moved ~13 arcmin across three observing nights. Intra-exposure motion estimated from the displacement between same-night pairs is ~0.02 arcsec over 30 s, much smaller than the recovered separation.
- At the first and last i-band epochs the PA changed by 49.46° in 11.963 days. The six-pixel and broadened-PSF endpoints differ from the nominal change by <0.06°. Difference images give 51.41° and individual PAs up to 5.75° from the science estimates. Difference and science share exposure pixels. Each same-night i/z pair used one detector. Across the three nights the object landed on three different detectors.
- Same-night i/z separations differ by ~0.04–0.06 arcsec even though the angles agree within ~1.4°. The exposures are less than 42 minutes apart. PSF systematics, filter-dependent sensitivity and noise should be considered alongside formal errors. An assumed 0.06-arcsec independent error in each coordinate still yields χ² = 62.05 for 8 degrees of freedom under a fixed relative vector, after minimizing over all relative label assignments. With a 0.10-arcsec floor, χ² = 23.81. These are sensitivity calculations and do not establish a survey-wide significance.
- Each control pool comprises other solar-system detections in the target visit with S/N > 7, reliability > 0.9 and the listed flags clear. The five pools contain 600–642 peers each. We selected up to three same-detector peers by closeness in log S/N and log PSF size, filling to eight per visit from all detectors. All 40 images are included; 11 share the target CCD. Their nominal science Δχ² had median 4.81 and reached 931. Two visits each contained one peer with a larger improvement than the target. These controls may include motion trails, blends and bright detections. A per-visit improvement alone does not establish a binary.
- Template photometry within a four-pixel radius gives ~37–417 nJy versus 19,000–22,000 nJy in the science two-component fit. The last template has a positive residual (417 ± 65 nJy diagonal-noise estimate), about 2.1% of the science total. This can bias relative astrometry and is one reason to use the science frames as the primary measurement.
- Separate exploratory cutout checks on previously unconfirmed targets did not supply a new binary identification. For example, 2014 OH394 showed excess morphology in two successive difference images, but their science-image fits gave Δχ² ~1.7 and 5.8, and the next night's detection was consistent with one PSF. 2005 JA175 gave small improvements or inconsistent morphology across dates. Neither is advanced as a companion discovery. The frozen quantitative release concerns JR80 and its controls.

## Reproduce the frozen measurements

The tested image pipeline used Python 3.12, NumPy 2.3.5, SciPy 1.17.0 and Matplotlib. Run these commands in this directory.

    python verify_data.py
    python cutouts.py
    OPENBLAS_NUM_THREADS=1 python fit_images.py --workers 2
    python diagnostics.py

The pipeline runs offline and writes the astrometry, full fits, constant-vector stress test, model-sensitivity table and figure to the output directory. On other BLAS or SciPy builds, nearly degenerate solutions may vary in the least significant digits. The cutouts script accepts --refresh to download and verify the original packet pixels. Running verify_data.py with --refresh-dir /tmp/jr80-live saves the current ADQL responses outside the frozen sample. The database can change after retrieval.

[ALeRCE Rubin services](https://science.alerce.online/alerce-lsst/) · [ALeRCE TAP schema](https://tap.alerce.online/browse/alerce_tap/multisurvey_lowercase) · [Fink LSST cutout service](https://doc.lsst.fink-broker.org/services/api/solar_system/)

## Field rotation and empirical stellar PSFs

The camera orientation changed from 14.46° to 42.72° to 70.46° across the three Rubin nights. This aligns closely enough with the sky PA change to make a persistent detector anisotropy a serious alternative to relative orbital motion when only catalogue shapes are used. After averaging bands within a night and applying a systematic floor only once per night, a fixed sky vector has χ² = 39.04 for 4 degrees of freedom with a 0.06-arcsec floor. A fixed detector-frame vector has χ² = 5.15 for 4 degrees of freedom under the same test. Thus the time sequence alone cannot exclude an instrument-fixed pattern with that floor.

To check the actual local PSF, we independently rebuilt it from images of nonmoving sources on the **same visit and same CCD** as JR80. Fink reference-epoch flux identifies persistent sources. Five 500-arcsec cones in the Fink LSST alert database supplied the candidates. Fink's date selection also returns objects active across the time interval, so we imposed an exact visit and detector match on the returned records. We required science S/N ≥ 35, science flux ≤ 300,000 nJy, and template flux between 0.5 and 1.5 times the science flux. Within each visit we retained the twelve highest S/N candidates at most. The 52 retrieved images have S/N 38.8–421.7 and separations of 79–495 arcsec from the moving target. They are persistent point-source candidates, without independent catalogue classification as individual stars.

Pixel cuts check valid coverage, nearby secondary peaks, center offsets and aperture concentration. Two of 52 images failed the one-pixel centering cut. For the 50 retained references, we estimated central positions within 5.5 pixels using the packet PSF, sky-subtracted, resampled with cubic splines, normalized each star by its own aperture flux, and took the pixel-wise median on a common centered grid. The empirical profile's **shape** comes directly from the field pixels; the packet PSF enters only when estimating the position of each comparison star. Every target was then fit with the empirical PSF, and again with each star left out in turn. Each held-out star was fit as a putative double using the remaining reference stars. The input cone responses, full URLs, candidates, FITS, per-image digests and rejection flags accompany the results.

| UTC date in 2026 | Band | Reference stars | Empirical separation (arcsec) | Empirical PA (deg) | Δχ², one PSF minus two |
|---|---|---:|---:|---:|---:|
| June 27 | i | 12 | 0.632 | 55.69 | 813.5 |
| June 27 | z | 12 | 0.672 | 55.81 | 325.0 |
| July 7 | i | 8 | 0.545 | 94.51 | 306.7 |
| July 7 | z | 11 | 0.594 | 94.41 | 163.2 |
| July 9 | i | 7 | 0.580 | 104.84 | 821.7 |

The empirical i-band endpoint displacement is 49.15° in PA; the supplied packet PSFs gave 49.46°. All leave-one-star-out target PAs stayed within 0.65° of the corresponding full empirical result. This is a stability check against individual references, not a complete uncertainty interval. Local empirical PSFs preserve the displacement despite the camera geometry. One of 50 held-out reference sources mimicked a close pair under descriptive cuts of separation >0.45 arcsec, faint fraction >0.35 and Δχ² >100. It was on July 7 in i (DIA source 170635504327328187, fitted PA 98.4°, separation 0.456 arcsec, Δχ² 163); the target that night had Δχ² 307. None of the 19 reference sources from the i-band endpoint images met all three criteria. A few other sources show high fit improvement from crowding or anisotropy. The test does not prove that any comparison source is single, or rule out position-dependent effects between stars and the target. The independent HST measurement is especially important given the changing camera orientation.

Run the additional offline diagnostics:

    python stellar_sample.py
    OPENBLAS_NUM_THREADS=1 python empirical_psf.py
    python epoch_checks.py

