# Archive, identification and prior-report checks

Checked 2026-11-15. Positions use ICRS. OGLE-IV identifies OGLE-LMC-RRLYR-09022 at 05:15:25.99, −70:29:51.4, with identifiers LMC502.23.1177 and LMC113.2.62070. Its paired Gaia DR3 source is 4652104345088779264. EROS source lm0561m6445 measures a blended early light curve.

| Check | Finding | Primary source |
|---|---|---|
| OGLE-III LMC identification, RRc and remarks files | RRc, 0.3604030 d. No target entry in the remarks file. | [OGLE-III source directory](https://www.astrouw.edu.pl/ogle/ogle3/OIII-CVS/lmc/rrlyr/) |
| OGLE-IV LMC identification, RRd, anomalous RRd and remarks files | RRd, listed periods 0.3603430 and 0.4852720 d. Absent from the anomalous RRd list and remarks. The release supplies the second period. | [OGLE-IV source directory](https://www.astrouw.edu.pl/ogle/ogle4/OCVS/lmc/rrlyr/) |
| Gaia DR3 variability analysis and epoch data | RRab, 0.485765553 ± 0.000006959 d. AIP TAP supplies 41 raw G records, 39 after variability filtering. | [Gaia variability table](https://gaia.aip.de/metadata/gaiadr3/vari_rrlyrae/) · [Epoch photometry](https://gaia.aip.de/metadata/gaiadr3/epoch_photometry/) |
| EROS-II/390 | EROS Red/Blue fluxes provide separate pre-OGLE clocks from one instrument. | [EROS release](https://doi.org/10.1051/0004-6361/202557321) |
| SIMBAD at the coordinates and by alias | The catalog matches both OGLE aliases and the Gaia ID. Its seven bibliographic entries identified in this check are survey releases or population studies. | [SIMBAD object and bibliography](https://simbad.cfa.harvard.edu/simbad/sim-id?Ident=OGLE+LMC-RRLYR-9022) |
| Primary literature and exact-identifier searches | Queries used the quoted five-digit and four-digit OGLE designations, its field identifiers and Gaia source ID, combined with “mode switching”, “fundamental”, “period change”, and Magellanic RR Lyrae. The published LMC mode switch in the comparison papers is OGLE-LMC-RRLYR-13308. No source-specific report of the measured transition in 09022 was located. Publication priority remains open. | [OGLE comparison in section 7.4](https://ftp.astrouw.edu.pl/ogle/ogle4/OCVS/blg/rrlyr/pap1.pdf) · [Other measured period increases](https://acta.astrouw.edu.pl/Vol67/n4/pdf/pap_67_4_2.pdf) |
| Repository search | Main result paths and open PR titles were checked before the search. A later all-PR search for “09022” returned no match. | [Repository](https://github.com/dovvnloading/Open-Space) |

## Independent early and late data

OGLE, EROS and Gaia are independently acquired measurements. The two OGLE releases share the survey, and EROS Red/Blue and OGLE I/V share the exposures or instrument depending on the survey. The fit code does not pool timing or amplitudes across these systems. See the raw query list and Gaia neighbor cone in [measurements.json](measurements.json).

The NOIRLab Source Catalog DR2 reports object `190930_85409`. Its later 2018–2019 data are affected by crowding and extraction flags. SMASH DR2 resolves target `47663.1729184` and nearby `47663.1724344`, approximately 1.17 arcsec apart. The reported magnitude difference is about 0.85 mag in g and 0.17 mag in r. The typical NSC FWHM exceeds their separation. The post-2017 data do not supply a reliable later mode classification. DELVE DR3 MC includes the target at 0.039 arcsec in overlapping tiles, with duplicate object identifiers `0792m705.229800` and `0784m705.213513`. Each lists ten measurements from 2014–2016, with three g, two r, two i and three z. The two tiles are duplicates. [NOIRLab Source Catalog](https://datalab.noirlab.edu/data/nsc) · [SMASH](https://datalab.noirlab.edu/data/smash) · [NOIRLab TAP](https://datalab.noirlab.edu/tap)

At NOIRLab TAP, use the following for an independent check of the sparse DELVE detections.

```sql
SELECT mjd, mag, err, filter, exposure, chi, sharp, ra, dec, objid
FROM delve_dr3.mc_meas
WHERE objid = '0792m705.229800'
```

## Attribution

OGLE-III LMC RR Lyrae were published by Soszyński and collaborators (2009), and the Magellanic OGLE-IV collection by Soszyński and collaborators (2016), with subsequent catalog updates. EROS-2 light curves are from the EROS Collaboration via CDS/VizieR. Gaia observations are from ESA's Gaia mission and processed by the Gaia Data Processing and Analysis Consortium. The Gaia data were accessed through the AIP Gaia archive. This packet contains extracted measurements for reanalysis and retains the original source locations.
