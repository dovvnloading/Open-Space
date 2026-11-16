# DR15 source selection

The input is the native 5XMM-DR15 version 1.0 **source** table. It has one row per fitted source, including single-observation fields. Source IDs are reassigned between releases. Identify the same object across earlier releases using sky coordinates and the observation record. All fluxes below are observed catalogue estimates in erg cm⁻² s⁻¹. They depend on the release's energy conversion factors. Refit a count spectrum to measure a model-specific flux.

## Reproduction

Run from the mission directory. NumPy and Astropy are required. The uncompressed file occupies about 2.04 GB. The native FITS remains the authority for every published catalogue number in these outputs.

```sh
curl -fL -o 5XMM_DR15cat_source_v1.0.fits.gz https://heasarc.gsfc.nasa.gov/FTP/xmm/data/catalogues/5XMM_DR15cat_source_v1.0.fits.gz
echo '3b58ee2c7c8dbe86009b65e23c35a534  5XMM_DR15cat_source_v1.0.fits.gz' | md5sum --check
gzip -dc 5XMM_DR15cat_source_v1.0.fits.gz > 5XMM_DR15cat_source_v1.0.fits
python analysis/select_dr15.py 5XMM_DR15cat_source_v1.0.fits results
python analysis/control_recovery.py 5XMM_DR15cat_source_v1.0.fits analysis/controls_coords.csv results
python analysis/simbad_screen.py results/measured.csv results/boundary.csv results/rescue.csv results/bright.csv results/simbad.json
python analysis/erass1_screen.py results/erass1.json results/measured.csv results/boundary.csv results/rescue.csv results/bright.csv
```

## Selection function

The source-row cut flow counts successive intersections. We accept SUM_FLAG 0 or 1, EXTENT=0, DIST_NN ≥20 arcsec, STACK_DET_ML ≥30, EP_CTS ≥100, and 0 < RADEC_ERR ≤3 arcsec. Each of the first three EPIC band fluxes must be finite and nonnegative and their sum, F_SOFT, must be ≥10⁻¹⁴. These cuts select bright, point-like, comparatively isolated objects with constrained formal positions. They lose confused, young embedded, and low-count sources. A nearby bright point source and diffuse structure can still produce false entries at greater separations.

The measured-colour branch requires –1 < HR1 < 1 and –1 < HR2 < 1 with positive reported errors, HR1 + σ₁ < 0.4 and HR2 + σ₂ < –0.65. For adjacent bands with rates A,B, HR = (B–A)/(B+A). These particular boundaries emphasize declining emission above 1 keV and favor low absorption. The nominal one-error offsets are selection margins and do not have calibrated confidence coverage near a boundary. The secondary boundary branch uses the same HR1 selection but requires HR2 exactly –1. Its value denotes a censored measurement when the band-3 rate reaches zero. Lack of an observed 1–2 keV rate does not imply zero physical flux.

Two supplementary branches were defined after inspecting early colour-selected leads. The **rescue** branch admits objects with HR1+σ₁ < –0.5 (both strictly interior and positive σ₁) or HR1 = –1 irrespective of measured HR2. Band 1 must have flux ≥3 times its positive uncertainty in at least two of pn, MOS1, and MOS2. All three band 3–5 errors must be finite and positive, and the sum of (flux + 1.645 × reported error) for bands 3–5 must be less than half F_SOFT. This is a sensitivity screen. It assumes catalogue flux errors can approximate a margin, including where rates reach zero. It is not a joint Poisson confidence bound. The **bright**, absorption-tolerant branch requires F_SOFT ≥1.5×10⁻¹³, and either an interior HR2 < –0.35 with HR2 + σ₂ < –0.25 and positive σ₂, or HR2 = –1. If HR3 is interior and its error is finite and positive, require HR3 + σ₃ < 0. At least two cameras must have ≥3σ in either band 1 or band 2. The rescue queue excludes objects eligible for the first two colour cuts. The bright queue excludes objects eligible for any prior colour branch. The same hard-band and variability tests apply in each branch.

For both bands 4 and 5 we require finite nonnegative flux, positive finite uncertainty and flux <3σ. The test suppresses a detected hard tail. A large σ can still conceal one. The VAR_FLAG true value or a finite VAR_PROB below 10⁻⁵ vetoes a source. An absent temporal statistic is retained as untested. N_CONTRIB counts observations covering the position, including upper limits. Any APPROX_SOURCE_VAR value published for a source is a separate cross-mission warning that requires individual examination because passbands, spectral shapes, and confused beams can produce apparent variations.

Within each branch, records with a published Gaia, OM, or WISE match probability ≥0.8 sort after other records. We sort each group by decreasing F_SOFT and break ties by integer SRCID. Missing probabilities do not demonstrate blank sky. F_SOFT_ERR is the quadrature sum of the reported band errors and omits covariance and spectral calibration systematics. Bands 1–5 cover 0.2–0.5, 0.5–1, 1–2, 2–4.5, and 4.5–12 keV. The release STACK_NH and STACK_GAMMA describe a power-law detection model. SPEC_*_PL contains an absorbed power-law pipeline fit when available. None of these parameters identifies thermal emission.

## Rejection and cross-identification limits

The output records the number of entries after every cut. Test the pipeline against known, thermally emitting neutron stars using coordinates in controls_coords.csv. For each of these eight control names, control_recovery.py looks for the nearest 5XMM source within 60 arcsec. The seven historical isolated stars and the previously published XMM source 4XMM J022141.5–735632 are illustrative positive controls. Their targeting, spectrum, variability and count-rate distribution differ from the unknown population. This empirical recovery fraction cannot estimate an all-sky completeness or occurrence rate. The 5XMM footprint is targeted and heterogeneous. A population completeness model must account for exposure, Galactic absorbing column, detector response and background at each sky location.

SIMBAD searches retrieve every indexed object within 15 arcsec and retain the actual angular separations. For eROSITA, search the public eRASS1 main catalogue within 30 arcsec in the western Galactic hemisphere. Record all sources in the search cone, including provenance and uncertainties. Verify positional consistency, source density, angular resolution, exposure, spectral response, and XMM source multiplicity before attributing the eROSITA flux. Outside the public footprint the script records absence of coverage. A position without a catalogue source in the footprint requires an exposure-map-derived count upper limit before claiming a flux constraint. eRASS1 flux uses 0.2–2.3 keV and an absorbed power-law conversion.

DR15's pre-stacking astrometric correction was not propagated into its final positions. Its release guide estimates a systematic σ per coordinate of 0.88 ± 0.01 arcsec from a high-latitude sample. When an astrometric reference has negligible uncertainty, use σ² = RADEC_ERR²/2 + 0.88² as a first-pass circular approximation. The 95% radius is 2.448σ. Per-field shifts, PM offsets, bad pixels and crowded fields can invalidate this approximation. A blank Gaia/WISE record alone cannot supply an optical or near-IR upper limit. Independent local image limits, a search for bright saturated and moving counterparts, solar-system contamination when relevant, and source/background spectral extraction are required before an isolated neutron-star interpretation.

- [SSC release, product checksums and watchouts](https://xmmssc.irap.omp.eu/Catalogue/5XMM-DR15/5XMM_DR15.html)
- [SSC user guide, sections 3.4 to 3.10](https://xmmssc.irap.omp.eu/Catalogue/5XMM-DR15/5XMM-DR15_Catalogue_User_Guide.html)
- [HEASARC XMMSSC schema and release details](https://heasarc.gsfc.nasa.gov/W3Browse/xmm-newton/xmmssc.html)
- [SIMBAD TAP](https://simbad.harvard.edu/simbad/sim-tap)
- [eROSITA-DE DR1 footprint](https://erosita.mpe.mpg.de/dr1/)
- [HEASARC eRASS1MAIN column definitions](https://heasarc.gsfc.nasa.gov/W3Browse/erosita/erass1main.html)

The native run recovered two of eight controls. RX J0420.0–5022 and 4XMM J022141.5–735632 enter the measured branch. All six omitted controls satisfy its hardness cuts and have small but highly significant long-term flux changes, with FRATIO from 1.057 to 1.467 and VAR_PROB below 10⁻⁵. Three also exceed 3σ in a hard band. The variability and hard-band cuts systematically penalize extremely high-count objects. These eight controls also illustrate real offsets from static catalogue coordinates.

All 321 final source rows have a SIMBAD and eRASS1 search record. SIMBAD returned 382 neighboring entries for 257 rows. Of the 195 rows inside the published eRASS1 hemisphere, 121 have at least one eRASS1 entry within 30 arcsec, for 122 pairings. These counts are positional neighbors.
