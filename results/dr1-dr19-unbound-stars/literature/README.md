# Published candidates and independent measurements

**Ten of the 2,708 screened Gaia sources have a possible match in the literature index.** Six have an exact identifier explicitly labelled Gaia DR3 in a source. The other four have an earlier-release identifier or a positional association for review. All eleven sources in `radial_followup_gaia.csv` lack a match within this index. Absence from this finite index leaves their discovery priority open.

The union contains the 1,020 extreme DESI entries, the additional DESI moderate-speed screen, the 501-source Gaia-matched DR19 tail, and individual follow-ups. See `matches.csv` for the hits. Some indexed papers discuss targets later found to be bound, pulsating or contaminated. A hit establishes prior attention at the level specified by the cited paper.

## Coverage and reproducibility

`manifest.tsv` gives the source, release, URL and hashes for 21 published tables, compilations or full-text ID extracts. `index.csv` has 4,723 records, 4,346 distinct nonempty identifier strings and 3,741 records with coordinates. The index includes the following data.

- The archived Open Fast Stars Catalogue at commit `503548d5`, with 492 literature objects including Brown's HVS objects, position pairs and available historical aliases. Some objects have multiple catalogue identifiers.
- The 591-star Li et al. LAMOST/DR2 set, the Li et al. 2023 target and comparison tables, Liao et al.'s 519-source catalogue, Montanari's 1,649 candidates, Marchetti's DR2 tables, and Verberne et al.'s 600 GC-ejection prospects. Katz et al.'s high-RV validation sample includes refuted velocities.
- All 53 Gaia identifiers in Bhat et al.'s DESI candidate/contaminant Table D.1, DESI-312, DESI-HVS1, and ten J2016 positions for LAMOST-HVS5 through HVS14. Li et al.'s comparison table also includes S5-HVS1.
- Printed IDs from the Marchetti EDR3 and DR3, Igoshev, El-Badry, Scholz, Prudil, Hattori, Carballo-Bello and Han papers. These full-text extracts include comparison stars. They cover the machine-readable identifiers printed in the available PDFs, and can miss objects referenced only by name, coordinates, spaced digits or an online appendix.

For reproducibility, the following commands assume the working directory is the containing results directory. The index rebuild requires Python and Poppler's `pdftotext`. The matcher requires NumPy and SciPy.

```sh
python literature/build_index.py --cache /tmp/hvs-literature
python literature/crossmatch.py desi_screen.csv moderate/desi_screen.csv sdss/gaia_tail.csv literature/radial_followup_gaia.csv --out literature/matches.csv
```

Gaia identifiers are read as decimal strings. The matcher compares exact IDs first. For coordinate-bearing records with a known epoch, it propagates catalogue proper motions to J2016 where available and searches within 3 arcsec. For archival records with uncertain epoch it admits a review radius of 5 arcsec plus the candidate's proper motion accumulated over 25 years. Coordinate and legacy-ID leads need individual verification. Ambiguous DR2 split/merge associations and entries lacking both usable coordinates and a corresponding DR3 identifier can escape this search. Some papers have extended online selections beyond their printed targets. The index is a documented novelty screen, with no claim to exhaust all historical catalogues.

VizieR responses contain a query timestamp. The manifest preserves each observed response hash plus a hash of its nonempty, noncomment tabular content, which the rebuild checks. Other input hashes apply to complete downloaded bytes.

## Follow-up measurements

`followup/receipts.tsv` records the archived extracts and spectral file with hashes and requests. A 5.04-arcsec SIMBAD cone and an 8-arcsec LAMOST DR11 v2 search covered all eleven prioritized objects. LAMOST returned three observations of one target. The specialist RR Lyrae catalogue reports a systemic RV of +67.59 ± 2.37 km/s for Gaia DR3 1170190610327078144, alias LINEAR 13754273 (RRc). Its flag 0 indicates a measurement obtained by Baeza-Villagra et al. The independent LAMOST spectra have catalogued RVs of +67.24 ± 2.50, +89.31 ± 9.39 and +82.89 ± 7.11 km/s. The SDSS archive identifies Gaia DR3 2752545725694798336 as SDSS J001357.40+090521.6. Plate 4536, MJD 55857, fibre 740 has ZWARNING=0, an F3/F5V template and z=−0.0002808511 ± 0.0000161373 (about −84.2 ± 4.8 km/s). The public SSPP row has empty atmospheric parameters for this object. The archived DR17 spectrum accompanies the DR18 metadata.

The full AIP StarHorse EDR3 table returned four photo-astrometric estimates, in `followup/starhorse.csv`. The underlying query and XML response are preserved. The Gaia DR3 2752545725694798336 distance quantiles (16/50/84%) are 4.39/6.47/9.81 kpc. The pipeline used Gaia parallax with photometry and population priors. For the RRc star the StarHorse result is a generic stellar fit without an explicit pulsating-star luminosity calibration.

### Parallax-free photometric sensitivity

`followup/photometric_inputs.csv` contains clean stellar SDSS PSF photometry and/or the nearest Pan-STARRS DR2 mean PSF source within 3 arcsec of each of six DESI stars, with at least five PS1 detections. `followup/dwarf_distances.py` computes conditional distances from the Ivezić et al. (2008) single-main-sequence calibration (Eqs A2/A7) and DESI metallicity, capped at the metal-poor calibrator limit of −2.5 dex. The SDSS measurements use the catalogue's band extinctions. For PS1, the quadratic PS1-to-SDSS transformations use the observed PS1 g−r from Tonry et al. (2012) Table 6. The converted g/r/i extinction coefficients are 3.303/2.285/1.698 times the SpecDis E(B−V), from Schlafly & Finkbeiner (2011). These distance estimates depend on stellar class, reddening, and the assumed metallicity. They do not use the target parallax or a neural-network absolute magnitude. The equal-luminosity binary column reports a distance increased by a factor of √2 at fixed observed flux. Near the blue turnoff, the calibration allows substantial scatter in absolute magnitude. A subgiant can be brighter still.

| Gaia DR3 source | SDSS single-dwarf d (kpc) | PS1 single-dwarf d (kpc) |
|---|---:|---:|
| 1717410430709663616 | — | 6.67 |
| 2752545725694798336 | 3.63 | 3.46 |
| 4435527217525195264 | 0.86 | 0.87 |
| 4444095982221534720 | — | 9.00 |
| 4550108286326235520 | — | 0.98 |
| 6340918677332110720 | 2.93 | 3.06 |

Sources [Ivezić et al. 2008](https://doi.org/10.1086/589678), [Tonry et al. 2012](https://doi.org/10.1088/0004-637X/750/2/99), [Schlafly & Finkbeiner 2011](https://doi.org/10.1088/0004-637X/737/2/103), [StarHorse](https://gaia.aip.de/metadata/gaiaedr3_contrib/starhorse/), [Baeza-Villagra et al. 2026](https://doi.org/10.1051/0004-6361/202659088).

### White-dwarf reclassification

Swan et al.'s DESI DR1 white-dwarf catalogue matches 52 of the 1,020 extreme-screen objects, seven within the 13-object nearby subset. Four nearby stars have magnetic DAH spectra without a specialist RV. Three DA stars have published line-fit RVs of −66 ± 274, −36 ± 547 and +151 ± 42 km/s. White-dwarf apparent Doppler shifts also include gravitational redshift. These specialist fits use the same DESI photons as the MWS pipeline. The catalogue classification and fit uncertainties provide direct checks on the choice of stellar template. The downloaded source hash and reproduction are in `followup/white_dwarfs.py`; the returned rows are in `followup/white_dwarf_matches.csv`.

```sh
python literature/followup/dwarf_distances.py
python literature/followup/white_dwarfs.py desi_screen.csv --catalogue /tmp/DESI_DR1_catalogue_1.0.csv.gz
```

[DESI white-dwarf catalogue](https://zenodo.org/records/23453342) · [Survey methods and RV caveats](https://arxiv.org/abs/2609.04314)
