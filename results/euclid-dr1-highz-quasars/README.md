# Euclid DR1-Foundation: high-redshift quasar dropout pilot

**Result — no new high-redshift quasar is validated.** A frozen small-area selection in the **public 2026-11-12 DR1-Foundation** MER catalogues retrieved **37,443 broad VIS-faint/NISP colour rows**, of which **290** satisfy additional flux-significance, roundness-proxy and flag cuts. These are *photometric selection artefact/dwarf/galaxy/quasar possibilities, not 290 quasar candidates*. An independent replay reproduced both counts. The top five rows sorted by J significance have extremely blue **J−H ≈ −2.3 to −2.7**, a strong warning against interpreting them as high-z quasar discoveries.

**Positive control:** the already spectroscopically known **EUCL J125308.55+705432.3 (z≈7.69)** is found within **0.0512 arcsec** of its published position and passes the same strict cut. The published z≈7.77 and z≈7.61 objects were not matched in these particular public DR1-F cones; that is *not* a failed recovery or evidence against their published identity.

## Release, identifiers and prior art

- **Euclid DR1-Foundation**, public **November 12, 2026** via ESA/IRSA; about **1,900 deg²** of Wide Survey imaging, plus deep/calibration data. This is the public *release footprint*, not the area inspected here. [IRSA overview](https://irsa.ipac.caltech.edu/data/Euclid/docs/overview_dr1.html), [Euclid release](https://euclid.caltech.edu/page/euclid-dr1-foundation).
- Queried MER product tables: `euclid_dr1_mer_catalogue_wide_survey` — **10.26131/IRSA673**; `euclid_dr1_mer_catalogue_deep_survey` — **10.26131/IRSA675**. **Erratum:** `freeze.json` initially recorded IRSA674, which actually identifies Wide *Mode*, not Wide *Survey*. The prospectively frozen file is retained unchanged for SHA-256 audit; the query table name and both correct product DOIs are recorded in `s3_output/summary.json`.
- The user's historical z=7.64 record has already been surpassed: **Yang et al. (2026), [DOI:10.1051/0004-6361/202658883](https://doi.org/10.1051/0004-6361/202658883)** spectroscopically reported 31 new Euclid quasars, including **EUCL J172902.75+641018.1 at z≈7.77** and the z≈7.69 positive control. These are *known*, not new finds. Their survey examined about 3,000 deg² beyond the new DR1-F catalog footprint.

## What was actually computed

A single GitHub Actions job, triggered **only on pushes to `research/euclid-dr1-highz-quasars` that modify the frozen inputs or its workflow**, performed all networked data access and calculations; changes under generated `output/` and `s3_output/` are deliberately excluded from its trigger. The original IRSA TAP attempt timed out (7/7 science queries, and 6/6 control queries failed or returned 502), preserved in [TAP receipts](output/requests.json). This is **not** a null observation.

The successful fallback used the public IRSA **DR1 HATS Parquet dataset on S3**, retrieved and checked its two partition indexes, and scanned **eight actual Parquet source partitions**. It applied the fixed `freeze.json` seven-cone search, including two deep cones; **six cones intersected HATS partitions; the intended Fornax deep cone had none**. Sum of requested cone areas was **2.014698 deg²**, but the effective imaged area is **not independently measured** and no full-DR1 quasar abundance bound is justified.

The broad cut required NISP J flux >0.65 µJy, redder Y/J, blue J/H, VIS_DET=0 and low catalogue spurious probability. The fixed strict screen required J≥6σ, H≥4σ, |VIS|<2σ, Y−J≥0.40 (or lower limit), J−H≤0.45, clean NIR flags, ellipticity<0.3. No lower J−H bound was frozen, so the leading rows are pathologically blue. An **exploratory, post hoc** J−H diagnostic found **41/290** below −0.5 mag and **249/290** within −0.5 to +0.45; this is *not* a calibrated quasar filter.

## Reproduce and audit

- [Freeze](freeze.json), [S3 HATS pipeline](run_s3.py), [original attempted TAP pipeline](run.py), [selection tests](tests.py), [independent compact-data referee replay](referee.py).
- [Full S3 summary](s3_output/summary.json), [reviewer replay](s3_output/referee.json), [all 37,443 broad records](s3_output/broad_rows.csv), [five *unvetted* SNR-ranked proxy rows](s3_output/shortlist.json), [partition/object receipts](s3_output/partitions.json), [source metadata SHA-256](s3_output/source_receipts.json), [output hash manifest](s3_output/hashes.json).
- Prospectively frozen `freeze.json` SHA-256: `4e427c201cc7ca4f643414e361a1b71bdde2202625fa0024551ad1f95c0ee882`. The publisher's HATS partition-index SHA-256 values are in source receipts, and individual S3 Parquet object identifiers/ETags are in the partition log. **Multipart S3 ETags are not equivalent to full-file SHA-256;** the full Parquet bytes were not independently rehashed. Deposited outputs and strict-count replay are SHA-256 checked.
- [Independent referee report](REFEREE.md). Runner attempts are linked in [PR #206](https://github.com/dovvnloading/Open-Space/pull/206).

## Evidence gates: completed vs pending

**Completed:** public release and product verification; publication/record correction; frozen colour screen in new DR1-F catalogues; actual partition reads and receipts; Y/J/H and signed VIS selection with independent replay; one z≈7.69 known-object positive control; basic shape/flag screens.

**Not completed / no discovery claim allowed:** actual NISP point-source PSF validation; simulated filter-response quasar tracks calibrated to **real** L/T/Y brown-dwarf spectra with measurement scatter and base rates; red galaxies, ghosts, persistence, cosmic rays and multi-exposure image inspection; an effective footprint/injection-recovery model; W1/W2 unWISE forced photometry, Gaia/CatWISE proper motion, verified Legacy/DES/Pan-STARRS image limits; all-source literature vetting and **DESI DR1, SDSS, ESO, JWST/MAST** spectra, plus **eRASS1, VLASS, LoTSS, RACS** crossmatches. The toy likelihood field `toy_P_Q_given_Q_or_dwarf_UNCALIBRATED` in the convenience shortlist uses invented illustrative dwarf centres and a 1:1000 prior: **it is not a defensible quasar-versus-dwarf probability**. The calibrated probability for each row remains **unknown**.

Per [docs/REMARKABLE_FINDINGS.md](../../docs/REMARKABLE_FINDINGS.md), **no `[REMARKABLE]` tag or log entry**: the discovery bar is not met.

This work has made use of the Data Release (DR1) data from the Euclid mission of the European Space Agency (ESA).
