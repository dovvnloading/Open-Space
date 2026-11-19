# Euclid DR1-Foundation wide-survey VIS-dropout follow-up

**Claim — no new high-redshift quasar is established by this investigation.** This is a systematic, partition-complete attempt to expand [the 2.01 deg² pilot](../euclid-dr1-highz-quasars/README.md) to the public Euclid DR1-Foundation MER **Wide Survey** HATS catalogue. The scan and downstream validations must complete before reporting any survey counts, recovery fraction or new candidates. A photometric candidate list alone is **not** a discovery.

## Dataset, fixed selection and coverage

- Public **2026-11-12 Euclid DR1-Foundation** release, approximately **1,900 deg²** nominal imaging area; IRSA Wide Survey MER catalogue, **DOI [10.26131/IRSA673](https://doi.org/10.26131/IRSA673)**; [IRSA HATS data](https://irsa.ipac.caltech.edu/euclid-archive-documentation/data-products/). Deep Survey **10.26131/IRSA675** is *not* added to this wide-only denominator. The HATS wide catalogue comprises 888,152,165 published rows, not 888 million independent clean sky sources.
- The pilot's *unchanged* strict catalogue selection: J ≥ 6σ, H ≥ 4σ, |VIS| < 2σ, Y−J ≥ 0.40 (with specified 2σ Y lower-limit handling), J−H ≤ 0.45, VIS_DET=0, NIR J and H clean flags, ellipticity < 0.30, with the pilot's broad prefilter and no added lower J−H cutoff. Report the very-blue J−H tail separately, rather than retroactively changing selection.
- **Morphology:** the MER `POINT_LIKE_PROB` model applies only to VIS detections; using it for VIS dropouts would be invalid. Record the independently published NIR-source shape diagnostic `-3.25 ≤ MUMAX_MINUS_MAG ≤ -2.65` plus ellipticity < 0.2 from [Euclid's ultracool-dwarf study](https://doi.org/10.1051/0004-6361/202556980), explicitly as a *diagnostic*, not validated PSF morphology for these images.
- **120 dated known-quasar controls:** 31 spectroscopically published discoveries from Yang et al. (2026), [A&A 711 A104](https://doi.org/10.1051/0004-6361/202658883), frozen in [controls.tsv](controls.tsv), **plus 89 previously known z > 6.5 quasars** with exact published RA/Dec from [Bañados' public 18 March 2026 quasar census](https://github.com/banados/reionization-era-quasars), frozen in [census.tsv](census.tsv). The dated census plus July discoveries cannot establish that *every* November 2026 known quasar in the DR1 footprint was included; test all catalog-matched controls and report this limitation. Matching an HATS partition alone does not establish imaged-pixel completeness.

## Reproduction and evidence

A single path-filtered branch-only GitHub Actions job uses the runner's network, scans HATS parquet **row groups** in sequential resumed batches, and commits checksums, metadata, strict-selected row receipts and control matches under [output/](output/). Generated outputs never trigger a push scan. The job may dispatch its own next bounded run, never a repository-wide pull-request event. The worker halts rather than claiming survey completion on a failed partition.

After the scan, [vet.py](vet.py) ranks cautious colour/morphology survivors and queries public unWISE, CatWISE, Gaia and spectrum metadata when reachable; per-object query errors are recorded as errors, not nondetections. Image-level Euclid VIS ghost/persistence/edge checks require actual image pixels. No external query, model posterior, crossmatch, projected luminosity-function yield, or visual artefact verdict is asserted without its corresponding output receipt.

**Expected abundance benchmark, not a new QLF fit:** Yang et al. (2026) cite an approximate **one z≳7 quasar per 100 deg² to J≈23**, implying order **19 in a genuinely usable 1,900 deg²** before strict-selection completeness and masks. This is not a calibrated expectation for this selection (whose flux, colour and signal-to-noise selection function differs from the published reference sample), and cannot be compared to raw dropout counts as quasars.

Per [the discovery-level standard](../../docs/REMARKABLE_FINDINGS.md), this work is **not [REMARKABLE]**; no log entry is justified. Known discoveries retain their publication priority.

This work has made use of the Data Release (DR1) data from the Euclid mission of the European Space Agency (ESA).
