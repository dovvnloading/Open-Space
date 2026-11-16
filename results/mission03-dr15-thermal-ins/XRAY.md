# Pointed X-ray observations

## 5XMM J085618.5+380711

The October 2005 XMM-Newton observation 0302581801 targeted the cluster RXJ0856.1+3756. The serendipitous source lies 9.36 arcmin off-axis. The independently generated PPS source list places source 14 at (134.077036°, +38.119683°), 0.82 arcsec from the DR15 position. The data here are observation-level, pipeline-calibrated imaging event lists and extracted source and background pulse-height spectra, with their associated responses. No raw ODF reprocessing was performed.

Event selection used per-CCD standard good-time intervals, FLAG=0, single/double pn events or MOS patterns 0–12, a 15 arcsec source circle and a 25–45 arcsec background annulus. The band-1 exposure map weights both apertures; it is appropriate as an approximation for the 0.2–1.0 keV test because the great majority of aperture events have energies below 0.5 keV. The table reports total source-aperture counts. The original fit and event localizations derive from the same observation.

| Camera | 0.2–1.0 keV | Exposure-scaled local background | Counts below 0.5 keV | CCD good time | Four equal-good-time bins |
|---|---:|---:|---:|---:|---|
| pn | 133 | 28.36 | 114 | 24.92 ks | 33, 28, 34, 38 |
| MOS2 | 35 | 4.66 | 30 | 28.86 ks | 9, 15, 7, 4 |
| MOS1 | No exposure | | | | |

The pn counts occupy 28 detector pixels, including 22 double-pattern photons in the aperture. The MOS2 counts occupy 35 distinct pixels. Their event centroids fall 0.17 and 0.86 arcsec from DR15, respectively. The 8 arcsec aperture contains 85 pn and 21 MOS2 soft events. The concentration and agreement across two detectors support a compact source; a calibrated extent confidence bound would require off-axis PSF fitting. A four-bin Pearson test of the total aperture rate gives p=0.676 for pn and 0.060 for MOS2. Background variability and low counts limit these tests. Only one XMM epoch is available.

### Spectra and Galactic absorption

We fit the un-subtracted pn PHA with its own ARF, the exact redistribution matrix named in the header and an independent background PHA. Only QUALITY=0 channels enter the fit, giving 0.3–5.0 keV coverage in 28 predefined bins. There are 100 counts in the source PHA and 1,935 in the background PHA. The exposure and BACKSCAL ratio is 0.0175895. A separate Poisson background intensity is profiled in every bin. The reported deviances include both Poisson samples. These are exploratory model comparisons with the same data, binning and free parameter count. Deviance differences between these non-nested models have no automatic chi-square significance.

| Model | Spectral shape | Fitted column, 10²⁰ cm⁻² | Deviance |
|---|---:|---:|---:|
| Spherical blackbody | kT 54.3 eV (37.3–73.6) | 5.81 | 25.56 |
| Solar abundance MEKAL, zero redshift | kT 26.4 eV | 27.62 | 24.22 |
| Optically thin free-free continuum, constant Gaunt factor | kT 74.0 eV (45.8–123.6) | 8.61 | 25.74 |
| Power law | Γ 8.0 at upper fitting bound; lower limit 5.58 | 10.93 | 26.59 |

Parentheses give one-parameter approximate 90% profile intervals (Δ deviance 2.706). The MEKAL fit uses the original XSPEC additive table, each elemental abundance set to one solar, with all continuum and line emission included. Other fits use the Morrison–McCammon (1983) effective photoelectric absorption cross sections as does MEKAL here. The blackbody's absorbed 0.2–2.0 keV flux is 2.46 × 10⁻¹⁴ erg cm⁻² s⁻¹. The 0.2–0.3 keV contribution is extrapolated from the quality-good PHA, although this band contributes to the event and catalogue detections.

Bilinear interpolation of the 16.2 arcmin FWHM HI4PI map yields a Galactic atomic column of 3.08 × 10²⁰ cm⁻². This beam-averaged atomic measurement omits molecular and intrinsic absorption. Fixing the model column to this value gives kT = 57.9 eV and deviance 25.66 for the blackbody and kT = 67.0 eV and deviance 28.30 for solar MEKAL. The data admit both a thermal surface and a cool line-emitting plasma. A typical hard AGN continuum is disfavored by the steep spectral slope; soft-excess AGN, thermally emitting accretors, and a very soft plasma remain relevant. Spectra alone do not distinguish the broad population classes.

For the blackbody, the intrinsic photon spectrum is K E²/[exp(E/kT)−1] with E in keV. Equating K to 2π(R∞/d)² (1 keV in erg)³/(h³c²) gives R∞/d = 1.76 × 10⁻¹⁶, or 5.43 km/kpc, at the free-column best fit. A 10 km apparent emitting radius corresponds to 1.84 kpc; a 1 km cap to 184 pc. A 7,000 km full white-dwarf surface would require about 1.29 Mpc at this temperature and normalization. Smaller accretion-heated regions are possible. At fixed HI4PI column, the inferred radius per kpc falls to 3.58 km; the free-column best fit is model dependent. These are normalization examples, not distance measurements. NS atmospheric color corrections, emitting fraction and gravitational redshift affect the physical radius.

### eROSITA

The public eRASS1 0.2–2.3 keV upper-limit FITS tile 135051 gives one aperture count, 0.315 background counts, 72.47 s exposure and positional flag 0 at nested, NSIDE=65536 pixel 5360802870. The one-sided 99.87% UL_B is 1.485 × 10⁻¹³ erg cm⁻² s⁻¹ assuming the release's power-law conversion. The survey upper limit on the PSF-corrected count rate is 0.160 s⁻¹. The catalogue and original pixel do not furnish independent confirmation. Comparing reference-model survey fluxes to an ultrasoft EPIC flux requires a spectral conversion. The cumulative DR2 results and release-specific radii appear in the main investigation.

### Comparison fields

- **J210351.3−181038.** The spectral audit uses pn PPS source 2 from 0655090432, with a 0.39 arcsec positional offset. The associated optical source is documented separately. Among the three continuum models power law has the smallest deviance (22.98 vs 31.20 blackbody). The three DR15 observation IDs cover 0.84 days within the comet P/2005 U1 mosaic.
- **J180019.4+080504.** In 0109465301, the pipeline source lies 0.50 arcsec away, and the pn blackbody fit gives 77.4 eV (59.1–101.1 eV) or a free-free fit with similar deviance. There are 160 fitted aperture counts and 40.80 scaled background counts. The shorter XMM exposures have large flux uncertainties and the MOS coverage is partly truncated.
- **J190116.4−220057.** Observation 0127720501 places this location in the immediate field of the extremely bright supersoft nova V4743 Sgr. The nearest PPS source is 93.57 arcsec away. There are soft-band excesses in both imaging cameras above a flat-annulus background estimate. The local background is highly structured and the source lies at a pn chip boundary, with the three most populated detector pixels at RAWY 198–199. The uniform-background p-values in the JSON are inapplicable to the structured nova field. Independent point-source and local PSF-wing subtraction would be needed to validate this detection.

## Reproduction

Use Python 3.12 with NumPy, SciPy, Astropy and astropy-healpix. From the mission directory the following commands retrieve the byte-identical public products and recompute the lead's principal results. Each archive product URL, size and SHA-256 is in the [archive manifest](analysis/archive_manifest.csv). Numerical results are in `results/`. The product manifests, spectra and event files carry the original observation IDs.

```sh
python analysis/fetch_manifest.py products
python analysis/fetch_manifest.py products --group mekal
python analysis/fetch_manifest.py products --group hi4pi
python analysis/fetch_manifest.py products --group erass1
python analysis/event_apertures.py 134.077254 38.119831 products
python analysis/exposure_timing.py 134.077254 38.119831 products
python analysis/fit_pha.py products/P0302581801PNS003SRSPEC000E.FTZ
python analysis/fit_plasma.py products/P0302581801PNS003SRSPEC000E.FTZ products/mekal.mod
python analysis/fit_pha.py --fixed-nh .30810075 products/P0302581801PNS003SRSPEC000E.FTZ
python analysis/fit_plasma.py --fixed-nh .30810075 products/P0302581801PNS003SRSPEC000E.FTZ products/mekal.mod
python analysis/hi_column.py products/hi4pi.fits 134.077254 38.119831
python analysis/erass_pixel.py products/em01_135051_024_UpperLimitTab_c010.fits.gz 134.077254 38.119831
```

For a comparison field, specify its uppercase source tag and use a separate working directory. The spectral acquisition tool `pps_spectra.py` also searches public PPS lists dynamically by position and retains its match distance.

XMM pipeline and response documentation: [ESA Data Files Handbook](https://xmm-tools.cosmos.esa.int/external/xmm_user_support/documentation/dfhb/pps.html), [HEASARC pipeline response calibration](https://heasarc.gsfc.nasa.gov/docs/xmm/xmmhp_caldb_mirror.html). Absorption and plasma models: [Morrison & McCammon 1983](https://doi.org/10.1086/161102), [MEKAL table model](https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/help/mekal.html), [HI4PI Collaboration 2016](https://doi.org/10.1051/0004-6361/201629178). Survey calibration and limits: [Merloni et al. 2024](https://doi.org/10.1051/0004-6361/202347165), [Tubín-Arenas et al. 2024](https://doi.org/10.1051/0004-6361/202346773).
