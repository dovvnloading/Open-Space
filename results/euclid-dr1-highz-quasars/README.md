# Euclid DR1-Foundation: high-redshift quasar dropout search

**Claim (provisional; before data acquisition): no newly discovered or spectroscopically confirmed z > 7 quasar is established by this investigation.** This folder records a fresh DR1-Foundation search and its honest coverage, candidate and contamination limits. No `[REMARKABLE]` claim is warranted without discovery-level evidence.

## Verified release and prior-art correction

Euclid **DR1-Foundation**, released **2026-11-12**, is public via ESA and IRSA. Approximate public Wide Survey footprint is **1,900 square degrees** (about **500 north / 1,400 south**); DR1 includes Deep Survey and calibration areas. The complete DR1 science release is planned for 2027. Published release: [ESA explanatory supplement](https://euclid.esac.esa.int/dr/dr1/expsup/), [IRSA DR1](https://irsa.ipac.caltech.edu/data/Euclid/docs/overview_dr1.html), [IRSA DR1 MER Wide Survey](https://doi.org/10.26131/IRSA674). This is a release-wide figure, **not** a claim that the scripted search covers 1,900 deg².

The historical z=7.64 record was superseded in **July 2026**: [Yang et al. 2026](https://doi.org/10.1051/0004-6361/202658883) confirmed **31** Euclid-selected quasars at 6.6 < z < 7.8, including **EUCL J172902.75+641018.1, z≈7.77** and **EUCL J125308.55+705432.3, z≈7.69**. Any match to these must be classified **known / positive control**, never new. Yang et al.'s pre-DR1 search drew from ~3,000 deg² over 18 months; the public DR1-Foundation Wide Survey is smaller and overlaps, rather than constituting a blind previously unexamined sky.

## Status and provenance

- **Completed:** verified public release, DOI, published record and scientific-prior-art correction; read [remarkable bar](../../docs/REMARKABLE_FINDINGS.md).
- **Pending:** run frozen Euclid DR1 MER TAP selection on networked GitHub Actions, retain raw-query receipt and SHA-256, measure actual queried coverage, calibrate selection response versus brown dwarfs / galaxies, test matched known quasars and crossmatch independent archival data.
- **Not yet established:** a blind full-DR1 screen, secure Y/J/H colour-discriminating completeness, any novel quasar candidate, or any independent spectroscopy of a new object.
- **Reproduction:** scripts and exact query/config will be committed before the first science run; generated outputs are committed separately by the single scoped GitHub Actions job.

**Note:** A VIS non-detection is not itself a valid veto for NIR-only sources: DR1 `POINT_LIKE_PROB` is null for NIR-only objects and requires an independent NIR compactness check. Asymmetric optical limits, unmatched archive catalogues, and missing mid-IR coverage must never be treated as confirmed nondetections.

This work has made use of the Data Release (DR1) data from the Euclid mission of the European Space Agency (ESA).
