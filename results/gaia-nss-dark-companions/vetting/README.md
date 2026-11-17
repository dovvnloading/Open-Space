# Companion light and spectral tests

**HD 264291 (Gaia DR3 3378588057203660160) has a well-measured outer orbit and an unresolved companion classification.** The photocentre and RV amplitude disfavor a single, ordinary 1.77-solar-mass dwarf. A close pair of late dwarfs with the same total mass is compatible with the broadband SED. The spectra show weak cool-template responses whose strength depends on the assumed kinematics.

## Independent dynamical light test

For each of 32,768 Gaussian draws from the complete DR3 NSS orbit covariance, we fit a *free signed amplitude* and an offset to the LAMOST nightly RVs. Marginal likelihoods include the amplitude-normalization determinant. Each night has its formal uncertainty plus an independent 1 or 2 km/s floor in quadrature. We sample the conditional amplitude error and the asymmetric IsocLum primary mass error, and infer the dynamical mass from the independent RV amplitude and inclination. The separate node ambiguities leave two possible photocentre branches. Let B = M2/(M1+M2) and β be the G-band flux fraction of the companion system. Then β = B(1 ± a_phot/a1). Signed intervals are shown before imposing the physical interval [0,1]. They are diagnostic intervals, with no photometric likelihood or truncation prior.

| Nightly floor | M2, median [16th, 84th] in solar masses | a_phot/a1, median [16th, 84th] | β, primary-tracing branch, 2.5th to 97.5th | β, reversed branch, 2.5th to 97.5th |
|---|---|---|---|---|
| 1 km/s | 1.770 [1.707, 1.836] | 1.019 [0.939, 1.100] | −0.087 to 0.063 | 0.890 to 1.041 |
| 2 km/s | 1.774 [1.666, 1.892] | 1.045 [0.954, 1.141] | −0.113 to 0.063 | 0.894 to 1.062 |

The mass and light tests assume the measured RV tracks the visible star and retain the formal catalogue astrometry. Blends, chemical peculiarity and unmodeled systematics can affect the RV centroid and the primary mass. The reversed branch would have the measured star contribute little of the system's G light. [Calculation and effective sample sizes](photocentre.json)

## Optical through infrared

SIMBAD identifies the visible star as an A1 IV Sr Si chemically peculiar star. We use catalogued B, V and Gaia G, clean 2MASS JHKs and AllWISE W1. WISE W2 and W3 are consistent within 0.03 mag of W1; W4 is an upper limit. MAST returned zero GALEX catalogue sources within 3 arcsec for all three targets. An empty cone supplies no UV sensitivity limit. IR coordinates are propagated from Gaia to the approximate survey epoch. The full input tables retain quality flags and uncertainties. [Photometry and query receipts](photometry/inputs.json)

We compare an empirical dwarf-color primary with adjustable temperature, luminosity normalization and extinction to that primary plus one or two fixed-mass dwarfs. The empirical sequence is the author's 2022-04-16 tabulation based on Pecaut & Mamajek (2013). The flux scale uses the NSS parallax for HD 264291, 1.568 mas. The discrepant five-parameter single-star parallax is 2.161 mas with RUWE 4.70. The photometric uncertainties include 0.035 mag model/systematic floors, and 0.06 mag for B,V. The adopted approximate R_V=3.1 law and the primary's chemical peculiarity limit the interpretation of the color residuals.

| HD 264291 companion model | χ² / degrees of freedom | Best temperature, A_V | Companion light in V | In G | In Ks |
|---|---|---|---:|---:|---:|
| Zero light | 2.5 / 4 | 8100 K, 0.17 mag | 0 | 0 | 0 |
| One 1.77-solar-mass dwarf | 2.6 / 4 | 9300 K, 0.34 mag | 34.0% | 34.8% | 43.5% |
| Two 0.885-solar-mass dwarfs | 3.3 / 4 | 9237 K, 0.40 mag | 2.6% | 3.2% | 12.5% |

The single ordinary companion conflicts with both photocentre branches under these model assumptions. The lower-mass pair fits the primary-tracing branch. An SED alone is degenerate with primary temperature and extinction. For HD 192288 (1830677037678583680), adding a 4.05-solar-mass dwarf gives χ² ≈ 4.7 for three degrees of freedom, the same quality as a primary-only fit. That dwarf would contribute approximately 13% in V, G and Ks. This analysis cannot establish a dark companion in HD 192288. Its W1–W3 bands have contamination flags; only G, B, V, J, H and Ks enter its fit. [Model fits and magnitude residuals](sed_fits.json) · [Empirical dwarf sequence](https://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt)

## LAMOST spectra and injections

The actual vacuum, heliocentric blue spectra from two nights near opposite RV extrema have one dominant cross-correlation peak over 4950–5340 Å, at +33.11 km/s, with normalized correlation 0.987. Six subwindows give +30.0 to +35.5 km/s. Secondary local peaks in two windows do not reproduce in the full interval. At R ≈ 7500, one resolution element is roughly 40 km/s. This test measures differential line motion. It cannot establish a single star spectrum. [Window results and input product hashes](spectral_motion.json) · [Pixel cross-correlation curves](blue_ccf.csv)

The injection experiment uses eight DR11 spectra comprising four independent pairs of coadds. Each pair contains a night near each RV extreme. The synthetic component is a solar-metallicity PHOENIX 5400 K dwarf with log g = 4.5, broadened to the instrument resolution and either 6 or 40 km/s projected rotation. The two putative inner stars have equal flux and spectra. Their mean velocity follows the outer orbit and their illustrative component velocities relative to that mean fall in four scenarios, expressed as (low-night, high-night) half-separations in km/s: (0,0), (40,40), (0,80), and (80,120). These scenarios sample possible line shapes; they do not fit an inner orbit.

We align the observed primary lines in each pair, difference the normalized spectra, and project out ten local continua, slope, primary line strength, and first/second line derivatives to reduce wavelength and line-spread residuals. The red-arm repeat uses 6400–6800 Å with 6510–6615 Å excluded around Hα. At each trial velocity we mix the synthetic spectra with both real spectra at the specified *local continuum flux fraction*, and apply the same differencing and nuisance projection. The empirical detection threshold is the leave-one-out 95th percentile at 54 velocity offsets displaced at least 200 km/s from the predicted systemic velocity. This is a known-template, known-velocity sensitivity test. Correlated residuals at the true stellar line wavelengths can be much larger than the off-velocity residuals. The nearby velocity trials illustrate that limitation. No multiple-hypothesis detection significance is assigned.

Across four blue pairs, both rotations and all four velocity scenarios, 2.6%- and 3.6%-flux injections clear this threshold at every off-source trial. In the red arm 2.6% injections clear it in 50–100% of trials depending on scenario, and 5% injections clear it in all trials. The mass-pair SED predicts 2.6% in broadband V, which only approximates local blue-arm light. The actual unsplit 6 km/s template yields positive blue fractions of **1.8%, 0.5%, 1.1%, 2.0%** across the four pairs, and **1.8%, 1.3%, 2.4%, 2.8%** in red. At ±40 km/s splitting in both nights these change to **4.2%, 2.5%, 1.3%, 3.0%** in blue. For the widest splitting the target response is negative in all four blue and red pairs. Shared metal lines in the peculiar primary, epoch-dependent line shapes and wavelength calibration can give false responses near the systemic velocity. A specified synthetic signal can be recovered in shifted controls while the nature of the actual target remains ambiguous. [Initial blue experiment](injections.json) · [Other blue pairs](replicated_injections.json) · [Red repeats](red_injections.json) · [PHOENIX models](https://doi.org/10.1051/0004-6361/201219058)

We also masked every pixel within about one resolution element of a line exceeding 5% depth in either primary-aligned spectrum, leaving 1040–1159 blue pixels and 1636–1669 red pixels per pair. The unsplit, slow-rotation blue fits remain positive, 0.5–1.7%, with comparable on-target fractions in the red. This removes the strongest common features but leaves weak common metal lines. These conditional template signals prevent a clean spectral exclusion of the hierarchical alternative. [Weak-primary-line experiment](weak_line_injections.json)

## Catalogue cross-checks and reproduction

The 983 primary-screen stars yield four exact DR3-ID matches to the 12,426 LAMOST MRS SB2 candidates of Kovalev et al. (2024). The match table lists their identifiers and LAMOST designations. None exceeded the 3-solar-mass point-estimate screen. A two-arcsecond position match to the DR17 APOGEE double-lined VAC yields no matches among its 7,204 distinct names. Unobserved stars and undetected companions are part of both catalogue non-match sets. [Spectral identifications](sb2_hits.csv) · [Sample sizes and input hashes](sb2_depth.json) · [Kovalev et al.](https://doi.org/10.1093/mnras/stad3222) · [SDSS DR17 double-lined VAC](https://www.sdss4.org/dr17/data_access/value-added-catalogs/?vac_id=double-lined-spectroscopic-binaries-in-apogee)

Run after the Gaia screen, LAMOST extraction and the orbit fits, from this directory. The pixel scripts use the sibling FITS reader. They require NumPy and SciPy. Each raw spectral/model product URL and SHA-256 is recorded in the experiment outputs. Cache files are reusable, and the products are fetched from the originating survey or model grid.

```sh
python fetch_photometry.py
python sed.py --cache /path/to/cache
python sb2.py --cache /path/to/cache
python photocentre.py
python spectra.py --cache /path/to/cache
OPENBLAS_NUM_THREADS=1 python inject.py --cache /path/to/cache
```
