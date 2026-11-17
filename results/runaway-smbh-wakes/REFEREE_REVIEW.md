# Referee record — runaway SMBH wake search (PR #155)

## Decision and scope

**Accept as a scoped null, not a discovery or occurrence-rate constraint.**
No new line met the runaway-wake discovery bar. A published-control morphology was recovered, but the strict endpoint-to-host test **missed the known host in both filters**. The catalogue shape proxy can split real wakes into multiple entries or miss extended, low-surface-brightness emission altogether. **No population upper limit follows.** No [REMARKABLE] label and no docs/REMARKABLE_FINDINGS.md addition.

## Recomputed evidence from deposited inputs

From checkout root, execute the standard-library-only compact-data referee command:

    python results/runaway-smbh-wakes/referee_replay.py

The script performs independent group-bys, joins, counts and arithmetic and stops if a recorded conclusion does not match.

| Reproduction check | Recomputed number | Source data |
| --- | --- | --- |
| Freeze/footprint | 352 MER tile IDs, 8 Orion; **344** extragalactic tile-centre stamps, **272** with valid VIS/RMS | euclid/q1_tiles.csv; detector/tables/coverage.csv |
| Pixel denominator | **73,152,512** valid pixels = **0.056444839506 deg²**; conservative full-stripe center **24,610,032** pixels = **0.018989222222 deg²** | detector/tables/coverage.csv; detector/tables/summary.json |
| Positive/sign-inverted screen | **8,264** positive raw, **651** positive profile passes; **2,529** inverted raw, **450** inverted profile passes | detector/tables/coverage.csv |
| Injection recovery | **840** trials in **27** stamps, ten 84-trial cells, FWHM 0.20 arcsec; 6-arcsec outcomes **62,54,40,11,0** and 10-arcsec **64,58,47,23,2** at AB 23–27 | detector/tables/injections.csv; detector/tables/recovery.csv; detector/tables/recovery_by_sample.csv |
| Full parent-catalogue ledger | **114** HATS partitions, **29,953,430** rows, **22,349,951** VIS detections, **5,134** unique entries passing the shape cuts; **5,120** deep + **14** LDN1641 | catalog/hats_partitions.csv; catalog/hats_summary.json; catalog/hats_selected.csv; catalog/gaia_alignment.csv |
| Star-axis confounding | Deep selections **4,490** aligned to G<16 Gaia stars versus **92** perpendicular matches; **198** deep entries with zero quality flag | catalog/gaia_alignment.csv; catalog/hats_selected.csv |
| Individual VIS review | **68** reviewed segments; **18** diffraction spikes; **32** weak/confused and unresolved | detector/review_rows.csv |
| Blind archival area | HST **5.24288 arcmin²** plus **0.26205 arcmin²** of JWST *inside* the HST footprint, not additional blind area | archive/tables/field_totals.csv |
| HST/JWST outcomes | HST **637** raw F606W / **629** raw F814W, **104** cross-band pairs, **one** two-fit continuity pass that fails fixed-axis follow-up; JADES: **38** raw, **3** cross-band, **zero** bilateral continuity passes | archive/tables/field_totals.csv; archive/tables/archive_pairs.csv; archive/tables/western_followup.json |
| Published RBH-1 control | Recovered in both ACS bands; strict host distance **2.53 arcsec** F606W and **1.14 arcsec** F814W, both beyond 1 arcsec endpoint rule; F814W also fails peak-dominance veto | archive/tables/rbh1_measured.json |

### Checksum and receipt integrity

* Frozen tile-list SHA-256: bfd3724847646662e651363d4624728b69cb9a50577548320e5d7fddf89494b6.
* Frozen **700-row** cutout-receipt SHA-256: a987e71517e3513dabfa6b6b3dae46ccaf1b69cbb3ab6356bb0dd98c557c1e7b.
* Receipts cover exactly 344 VIS, 344 RMS, 4 Y, 4 J and 4 H cutouts. Every (tile, product) key is unique and includes byte count, geometry, source URL and hex SHA-256.
* All **688** VIS/RMS hash references in the 344 coverage-table rows match their individual cutout receipts; dimensions and tile identities also join without gaps. An earlier Windows-checkout review checked sampled frozen file bytes.

**Integrity boundary:** the git deposit does *not* contain the original full FITS cutout bytes or the full parent HATS/Parquet data. Matching recorded hashes establishes internal consistency of deposited records, **not** that remote upstream files were independently re-fetched and rehashed in this referee pass. This script replays the compact results, not native pixel/image or partition reductions.

### Limitations and decisive next tests

* The nominal **63.1 deg²** Euclid deep-field footprint is **not** the area searched blindly for linear features. The tile-centre stamps measure far less. The catalogue screen does not substitute for full-pixel search of every tile.
* A matched line, apparent nearby galaxy endpoint, or occurrence in sign-inverted images is not a wake classification or a false-alarm probability. Real-source sidelobes and chance alignments contaminate the controls.
* The straight synthetic filaments have fixed width/brightness/length templates. Recovery omits galaxy-host selection, individual-exposure persistence, rejection of disks/streams/jets, spectroscopic confirmation, fragmented/diffuse wakes, and unobserved sky. **No population or occurrence-rate limit** is possible.
* All **32** unresolved individual segments appear with their exact review image/panel, ambiguity, and case-specific deciding test in [detector/unresolved_followup.md](detector/unresolved_followup.md). Four priority cases failed to appear with >=3 local-contrast strength across 16 usable single VIS exposures; that is not exclusion of faint coadd structures.
* The archival HST western feature fails fixed-axis F775W, F850LP and targeted JWST continuity; RBH-1 demonstrates the host-selection false-negative. The astrophysical interpretation of RBH-1 itself remains contested.
* Remaining valuable tests are wider/overlapping stamps, independent exposures, deeper multiband imaging, resolved spectra and multi-model classifications. These are outside the completed measurements.

## Replay commands and CI

From checkout root:

    python results/runaway-smbh-wakes/referee_replay.py
    python -m open_space.verify --worktree

From results/runaway-smbh-wakes/detector, with NumPy/SciPy installed:

    python test_wake.py

The [referee CI](../../.github/workflows/runaway-smbh-wakes-referee.yml) additionally compiles deposited scripts and runs guard tests under Windows with a fresh D: temporary directory. Full native reductions require external caches/retrieval and fresh output directories and are not asserted by this verification.

**Verdict:** merge as a bounded, negative search with unresolved diagnostic cases. No [REMARKABLE] marker or remarkable-findings entry.
