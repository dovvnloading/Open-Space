# 1LHAASO dark/UHE counterpart census (2026 public-catalogue audit)

**Claim under test (NOT an identification):** do post-2024 public catalogues produce any counterpart excess or energetically credible new accelerator class for the 1LHAASO sources listed without a TeV counterpart in Cao et al. (2024)? Until the job writes verified output, **no new source identification or significance is claimed.** All assignments are *candidates for follow-up*, not established PeVatrons.

**Observing unit:** the 90 unique 1LHAASO sources, with independent KM2A/WCDA *component* rows. The high-energy subset requires a recorded >100 TeV detection (TS100), not just a large fitted flux. **Dark selection:** no source-level TeV association in the frozen Cao table; the stricter “no TeV, GeV, pulsar or SNR” label must be tested against named catalogues, not presumed.

**Controls / falsification:** for every candidate counterpart *class* and dark source, replicate the identical angular/extension matching rule on Galactic-latitude-matched longitude-shift controls. Report finite-sample p=(1+k)/(1+N), coverage masks, multiplicity and multiple-test issues; veto a physical claim if associations or energetics cannot be verified. A catalog coincidence is not an accelerator identification.

## Receipts and inputs (freeze on successful run)
- Primary 1LHAASO: Cao et al. 2024 ApJS **271, 25**, DOI [10.3847/1538-4365/acfd29](https://doi.org/10.3847/1538-4365/acfd29); CDS VizieR [J/ApJS/271/25/catalog](https://cdsarc.cds.unistra.fr/viz-bin/cat/J/ApJS/271/25); **180 component rows != 180 physical sources**. The analysis should persist the raw table and its SHA-256, then use the 90 unique names.
- Fermi [4FGL-DR4](https://fermi.gsfc.nasa.gov/ssc/data/access/lat/14yr_catalog/): the 2024-07-24 v35 correction must be distinguished from 2023 v32. Latitudes, angular sizes and flags matter.
- eROSITA-DE [eRASS1 DR1](https://erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/Catalogues_dr1/): **only the released German half-sky**, not global non-detections; DR1 version 1.2 released 2026-01-12.
- ATNF, 2024-25 Galactic SNR (Green/SNRcat), MeerKAT SMGPS, ASKAP/EMU, Gaia DR3 clusters, HAWC/H.E.S.S.: individual provenance and usable coverage must be demonstrated before quoting matches. Absence of access or a public authoritative table = **pending**, never zero matches.

## Completed checks
- [x] Confirmed primary paper/catalogue identifiers and that component-level rows must not be double-counted.
- [x] Reviewed repository discovery-label rules: `docs/REMARKABLE_FINDINGS.md`. **No [REMARKABLE] designation without validated discovery-level evidence.**

## Pending checks
- [ ] Freeze exact 1LHAASO bytes, hash, source/member counts and dark definition; confirm the 43 UHE sources.
- [ ] Load each counterpart catalogue with coverage, epoch, metadata, hash, class and position checks.
- [ ] Compute latitude-stratified shift controls with finite-N p, background density and multiple-test correction.
- [ ] Check pulsar spin-down/age/distance, physical extent and energy budget; compare literature proposals since 2024.
- [ ] Run independent referee pass (including failed/missing catalogues); update this README only with checked results.

## Reproduce
GitHub Actions workflow `.github/workflows/lhaaso-dark-counterparts.yml` will run on branch push changes under `results/lhaaso-dark-counterparts/`, excluding generated `output/`; **one job** accesses data and computes in the runner, and commits run outputs in bounded batches. Scripts and receipts will live in this folder. No private token is required to fetch public catalogues.
