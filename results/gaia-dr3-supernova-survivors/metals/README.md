# Composition-first search in DESI DR1

Three of the 1,452 quality-selected, rapidly moving DESI white dwarf candidates with metal, carbon or hot-helium spectral designations pass the 350 km/s transverse-velocity screening statistic. All three were already published. They are J1109+0001, J1240+6710 and the carbon-rich white dwarf LP 93−21. The type-Iax interpretation of LP 93−21 had already been proposed and subsequently challenged using pressure- and gravity-corrected radial velocities. This branch adds no new survivor.

## Completed selection

The parent contains 63,968 visually inspected Swan et al. DR1 catalogue entries, including 63,152 with exact Gaia DR3 IDs. A WD classification containing Q, Z or O anywhere after the initial D is retained. This includes DZ, DQ, DO and every mixed type such as DAZ, DQA and DBO. The O in DO conventionally denotes ionized helium. Class membership alone supplies no oxygen-abundance inference. There is no catalogue type-confidence or parallax-S/N floor. This exploratory branch was chosen to retain chemically interesting slower objects and uses known survivors for calibration. [Swan et al. 2026 and catalogue](https://arxiv.org/abs/2609.04314)

| Successive gate | Entries |
|---|---:|
| Qualifying spectral types across all designations | 2,917 |
| Exact DR3 ID and parallax/proper-motion solution | 2,913 |
| Proper motion >35 mas/yr | 1,612 |
| Also RUWE <1.4, ≥9 visibility periods, no duplication, ≤2% multi-peak fraction, IPD harmonic amplitude <0.1 | 1,452 |
| Positive raw parallax + 2σ and 4.74047 µ/(parallax + 2σ) >350 km/s | 3 |
| Repeating the last cut with valid individual zero-point calibration | 3 of 1,435 |

For the four additional objects whose catalogue designation is DR2, the official DR2–DR3 neighbourhood table returns one closest match per object at <1.7 mas, with proper motion 6.9–19.3 mas/yr. One DR2 source also has a more distant DR3 neighbour without a proper motion. None of these five possible links passes the proper-motion gate. The primary counts retain exact DR3 designations, with the DR2-only audit in dr2_check.csv.

The statistic uses a two-standard-error upper value for measured parallax and ignores proper-motion uncertainty. It serves as a screening quantity. It gives neither a 3D escape probability nor a calibrated posterior velocity quantile. The gate retains both known comparison stars actually present in the exact-ID catalogue, 2/8 LP 40-365-like reference objects overall. Neither the 11 D6 controls nor US 708 has an exact entry in the catalogue. The 17 matching-quality stars outside the individual zero-point calibration domain stay in the raw test and have a maximum raw statistic of 155.3 km/s. There is no all-Galaxy completeness or rate limit. Completeness also depends on type assignment, the DESI target selection, extinction, proper motion and measurement errors.

| Object | Swan label | Raw statistic, km/s | Corrected statistic, km/s | Prior identification |
|---|---|---:|---:|---|
| LP 93−21 | DQ | 482.4 | 481.8 | Carbon-rich massive WD. A published Iax hypothesis was challenged by later spectroscopy and halo-orbit analysis. |
| J1109+0001 | DZ | 412.8 | 408.8 | Published LP 40-365 star. DESI spectra are fitted in the sibling spectra directory. |
| J1240+6710 | DZ | 380.0 | 383.5 | Published oxygen-dominated candidate remnant. Its original analysis argued for an origin distinct from the LP 40-365 class. El-Badry et al. 2026 include it in their comparison group. |

Exact SIMBAD identifiers, including J1109's historical 2QZ designation, resolve for all three. Older SIMBAD spectrum types differ from DESI in two cases. Literature classifications, atmospheric studies, and the El-Badry reference table determine the novelty verdict. [Ruffini & Casey on LP 93−21](https://doi.org/10.1093/mnras/stz2176) · [Kawka et al.'s reanalysis](https://arxiv.org/abs/1910.13053) · [Gänsicke et al. on J1240+6710](https://arxiv.org/abs/2006.07381) · [El-Badry et al. 2023](https://arxiv.org/abs/2306.03914) · [El-Badry et al. 2026](https://arxiv.org/abs/2606.11293)

## DESI spectral checks

The Swan DR1 archive supplies one exposure and one corresponding three-arm coadd each for LP 93−21 and J1240+6710. Each spectrum package and its six component files have SHA-256 hashes in spectral_receipts.json. Coadds alone enter the fit to avoid counting the same exposure twice. The spectra use the native DESI vacuum barycentric wavelength solution. The fits retain this calibration without applying a second velocity correction. We fit the observed shift and distinguish it from the star's systemic velocity. No DESI catalogue RV is supplied for either object. [Source spectra and catalogue description](https://cygnus.astro.warwick.ac.uk/phsdaj/WhiteDwarfs/) · [DESI DR1](https://data.desi.lbl.gov/doc/releases/dr1/)

**J1240+6710.** The Mg II 4481 Å feature gives an apparent shift of −191 ± 37 km/s. Jointly fitting the weak Si II 6347/6371 lines gives −178 ± 35 km/s. The O I 6158 and 7773 blends give −128 ± 35 and −200 ± 33 km/s under single effective-wavelength approximations. Independent narrow regions show Mg II equivalent width 3.92 ± 0.38 Å and O I 7774 equivalent width 2.51 ± 0.61 Å. The Ne I 6402 region is 0.76 ± 0.43 Å in this single DESI exposure. Gänsicke et al.'s higher-quality UV, optical and atmosphere analysis remains the basis for Ne and other abundances. They find an apparent RV of −158 ± 12 km/s and a systemic RV of −177 ± 10 km/s after correcting the gravitational redshift. These values have different physical meanings. [Gänsicke et al. 2020](https://arxiv.org/abs/2006.07381)

**LP 93−21.** Separate C I fits in fixed ±25 Å windows around vacuum 4933.40, 5053.56 and 5381.83 Å give apparent shifts of +276 ± 26, +242 ± 25 and +181 ± 32 km/s. Varying the half-window to 13 and 40 Å changes the latter shift to +167 and +180 km/s. The 5053 region is especially sensitive to the broader continuum and neighboring structure. The difference between the 4933 and 5382 centroids is consistent with the known wavelength-dependent pressure shifts. Kawka et al. infer a +15 ± 21 km/s systemic velocity after correcting pressure and gravitational shifts. A single carbon-centroid RV or default survey-template RV would be unsuitable for an orbit. The object is a known massive DQ on a retrograde halo orbit. [Kawka et al. 2020](https://arxiv.org/abs/1910.13053)

The custom C I fitter estimates a local linear continuum, Gaussian depth, centroid and width, minimizing inverse-variance-weighted residuals from multiple starts. At each wavelength, only unmasked finite pixels with positive inverse variance contribute. A formal per-line uncertainty comes from the Jacobian and residual variance. J1240 uses the same masked line and equivalent-width fitter as J1109, converting air laboratory wavelengths to vacuum. Pressure shifts, continuum placement, unresolved multiplets and WD gravitational redshift are systematic uncertainties outside those errors. Detailed non-LTE abundance modelling and a new origin inference for LP 93−21 are outside this measured catalogue result.

## Reproduction

First populate the Swan and Gaia cache with the sibling desi/fetch.py. Run the following with the same cache directory, Python, NumPy, pandas, SciPy and gaiadr3-zeropoint 0.1.0.

```sh
python select.py /path/to/cache
python identify.py
python fit.py /path/to/cache
```

The selector verifies the catalogue hash and all 288 Gaia batch hashes from desi/receipt.json before merging exact integer-string IDs. Query text and receipts for both SIMBAD and the DR2–DR3 check are included.
