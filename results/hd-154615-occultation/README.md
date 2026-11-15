# HD 154615 shows a resolved companion candidate and a changing eclipse clock

On 7 July 2022 the asteroid (32772) 1986 JL briefly cut off light from two sources toward HD 154615. The four geometries allowed by Dave Herald’s observation put the two sources 21.0–21.9 milliarcseconds apart. At the Gaia parallax that is 11.33–11.81 au in projection. The star is also a 1.851-day eclipsing binary. Four TESS visits show both eclipses arriving late in the middle of the six-year baseline. A third star orbiting the eclipsing pair could account for both the spatial and timing measurements. The outer orbit and whether the two resolved sources are bound remain open.

The occultation and eclipsing-binary classification were published separately. This contribution connects them, tests whether the resolved stars can be the eclipsing pair, and measures correlated timing changes in independent space photometry. [Herald et al., NASA PDS occultations](https://doi.org/10.26033/ehqs-jp27) · [Mowlavi et al., Gaia eclipsing binaries](https://doi.org/10.1051/0004-6361/202245330) · [MIT Quick-Look Pipeline](https://tess.mit.edu/qlp/)

## The geometric test

The selection used every entry in the PDS v4 asteroid event table with a first double-star position-angle solution. The 9,645 events contain 136 qualifying occultations of 135 distinct Gaia sources. We joined these by exact source ID to the Gaia DR3 eclipsing-binary variability table. All three matches are below. The 132 other IDs include sources absent from Gaia’s variability catalog. This selected sample cannot estimate the frequency of triple stars.

| Star | Occultation | Gaia eclipse period | Allowed projected separation at nominal parallax | Interpretation |
|---|---|---:|---:|---|
| HD 154615, TIC 126293345 | (32772) 1986 JL, 2022-07-07 | 1.85117 d | 11.33 or 11.81 au | Candidate triple with independent timing curvature |
| ATO J265.2097−11.9936, TIC 417476835 | (1309) Hyperborea, 2021-05-28 | 1.43860 d | 8.95 or 30.27 au | Candidate triple if the cataloged orbital interpretation holds |
| HD 167971, MY Ser | (713) Luscinia, 2019-09-20 | 3.32166 d | 40.15 or 40.70 au | Known hierarchical triple, positive control |

For a bound Keplerian double with eccentricity below one, the instantaneous separation is less than twice the semimajor axis. Thus even the most permissive mass required to put both occulted points in the short-period orbit is greater than $(s/2)^3/P^2$ solar masses, for projected separation $s$ in au and period $P$ in years. The smallest HD 154615 solution gives **7.07 million solar masses**. The corresponding floor for the ATO source is **5.77 million**. These bounds rule out one ordinary stellar orbit accounting for both measurements. The gap survives large errors in the measured separation or distance. The PDS formal separation errors are missing for several geometries, including all four for HD 154615. The alternatives are correlated geometric solutions of a single positive chord per target. [All measured alternatives](separations.csv) · [Independent interferometry of the positive control](https://doi.org/10.1051/0004-6361/201629260)

Dave Herald recorded two separated dips toward HD 154615 at 0.04-second cadence. The fainter component’s dip was short, with only a few deep samples. The flux almost entirely recovered between them. Wayne Hanna recorded a miss. For ATO J265.2097−11.9936, Dave Gault recorded nested drops at 0.32-second cadence. The inner interval reached background, a different profile from the two distinct HD dips. The original VizieR time series and PDS contacts are included. The measured count levels contain retained background, so the counts by themselves do not give a calibrated stellar flux ratio. [HD light curve](inputs/occ_32772.txt) · [ATO light curve](inputs/occ_1309.txt) · [Observer contact times](inputs/contacts.csv)

## Independent eclipse timings in HD 154615

TESS observed the target in sectors 12, 39, 66 and 93, approximately two years apart from 2019 to 2025. We measured each well-covered minimum in the MIT QLP light curves with a Gaussian dip and a local sloping baseline. The main aperture yielded 41 primary and 41 secondary times. The median primary depths by sector were 4.8–5.1% of aperture flux and the secondary depths were 2.7–2.9%. The measured mean period across the baseline is approximately 1.85108 d. For **each** eclipse type we set a straight clock through its first and last sector medians, then compared its two middle visits.

| TESS sector | Primary delay | Secondary delay | Conditional 95% intervals, primary / secondary |
|---|---:|---:|---|
| 39 (2021) | +9.63 min | +8.83 min | +6.58 to +10.46 / +3.43 to +11.80 min |
| 66 (2023) | +6.40 min | +6.99 min | +3.87 to +8.39 / +1.62 to +10.75 min |

The intervals resample whole measured eclipses 20,000 times within each sector and eclipse type, holding the sector mean cycle numbers fixed. This captures event scatter while conditioning on the adopted light curves, fits and selected visits. Across the raw simple-aperture flux and three detrended aperture choices, all primary delays lay between 9.52 and 9.91 minutes in sector 39 and between 6.40 and 6.53 minutes in sector 66. The secondary delays ranged from 8.55 to 8.89 and 6.52 to 6.99 minutes respectively. An independent shape check learned an empirical eclipse template from sector 12 alone and shifted it against each later visit, fitting the local flux slope and depth. It gave +8.38 and +5.54 minutes for the primary and +7.01 and +7.55 for the secondary. These aperture products and shape tests reuse the same photons. The primary and secondary fits move together, as expected for a light-travel delay or a change in the inner period. Apsidal motion by itself generally moves the two timings in opposite directions. [Individual minima](timings.csv) · [Every aperture comparison](timing_summary.csv)

We also formed primary and secondary out-of-eclipse minus in-eclipse images from the original 9 × 9 TESS cutouts. Fits lie 0.07 pixel from the target in sector 39 and 0.21–0.23 pixel in sector 93. Within 108 arcseconds the brightest other Gaia DR3 source with a measured G magnitude is 4.31 mag fainter than HD 154615 and 51.2 arcseconds away. The closest is 4.83 arcseconds away and 9.91 mag fainter. These checks disfavor the resolved field neighbors as the source of both eclipses. The Gaia catalog itself does not resolve a 21 mas double. Its parallax of 1.8539 ± 0.0296 mas comes from a single-source solution with RUWE 1.71. [Pixel positions and source URLs](localization.csv) · [Difference images](difference_images.csv)

ASAS-3 and Gaia recover eclipses in earlier observations extending back to 2001. The inferred mean phase shifts depend on the period carried across the intervening years. [Historical photometry, timing estimates and alias limits](early-photometry/README.md)

A chance optical alignment, an error in the single-chord double-star interpretation, and an intrinsic change in the eclipsing binary’s period remain possible. The four timing epochs do not determine an outer period or companion mass. A new resolved position near the reported 21 mas separation, combined with repeated systemic radial velocities, would test the physical association and distinguish outer motion from intrinsic period change. The four published 2022 position angles are 232.5°–264.9°, east of north.

## Second sightline

Gaia and ATLAS already identified ATO J265.2097−11.9936 as variable. Gaia’s broad, similar minima give the 1.438603-day classification. In TESS sector 92 the signal at twice that orbital frequency localizes 0.13 pixel from the Gaia position. Other periodic structure is present in this crowded aperture. A Gaia source 2.85 arcseconds away is 4.04 mag fainter, and a star 67 arcseconds away is brighter than the target. The broad modulation provides weaker eclipse evidence than the distinct unequal minima of HD 154615. A single TESS visit supplies no multi-year timing test here. The 6.8 and 23.0 mas occultation solutions therefore remain separate alternatives on a candidate sightline. [Gaia measurements](inputs/gaia_epoch_4162228724694601600.csv) · [Field sources](inputs/neighbors_4162228724694601600.csv)

## Coverage, provenance and replay

The PDS event table covers asteroid occultations through 2023, with partial late-year coverage. The PDS bundle total also includes planetary events, which this join excludes. The VizieR B/occ light curves come from Herald, Gault and colleagues. The photometry is from ESA/Gaia/DPAC and NASA/MIT TESS QLP via MAST. The two sets of TESS reductions share spacecraft exposures. Time stamps in the QLP extracts are BJD(TDB) minus 2457000. Gaia epoch photometry uses BJD(TCB) minus 2455197.5. Sky coordinates are ICRS and occultation dates are UTC. The small catalog extracts and observation series preserve their archive identifiers in the inputs or manifests. The Gaia membership queries in the manifest sort the IDs for readability. A server can return matching rows in a different order on replay. [Archive URLs and downloaded byte hashes](sources.json) · [Input checksums](inputs.sha256)

Run the following from this directory using Python, NumPy and SciPy. The first command verifies every frozen input file, including the membership of the search sample. The next commands regenerate the tables and the conditional timing checks.

```sh
sha256sum -c inputs.sha256
python analyze.py
python checks.py
```

For the independent pixel check, install Astropy and run `python pixels.py --cache /path/to/cutout-cache`. This downloads three TESScut requests and rebuilds both pixel tables. Generated FITS headers and ZIP timestamps can change on retrieval. The per-cutout hashes refer to the retained quality-zero FLUX array in big-endian float32 order. Pixel positions are zero-based. The manifest lists the original QLP FITS files. The compressed CSVs project TIME, QUALITY, SAP_FLUX, and the default, small and large detrended flux columns with all rows retained. QLP renamed KSPSAP_FLUX to DET_FLUX between these releases. The PDS screen projects the identity and four geometric alternatives from every nonempty DOUB_PA_SOLN_1 in Asteroid_2024Mar.psv. The two VizieR records retain the five complete source lines starting at the dated observation header.

Prior-work checks on 14 November 2026 (US Eastern) included SIMBAD bibliographies and coordinate searches within 30 arcseconds in the VizieR editions of the Multiple Star Catalog (J/ApJS/235/6/catalog) and Washington Double Star Catalog (B/wds/wds). The published multiple-system catalog recovered HD 167971 and returned no entry at the other two positions. Searches used the HD, TIC, Gaia and UCAC identifiers with “triple”, “tertiary”, “occultation” and “eclipse timing”. The two occultations, source variability and Gaia periods are prior measurements. We found no earlier report joining the two scales for these candidates or measuring HD 154615’s common timing change. This search cannot establish publication priority. [Multiple Star Catalog](https://doi.org/10.3847/1538-4365/aaa1a5) · [Original B/occ light curves](https://cdsarc.cds.unistra.fr/viz-bin/cat/B/occ)
