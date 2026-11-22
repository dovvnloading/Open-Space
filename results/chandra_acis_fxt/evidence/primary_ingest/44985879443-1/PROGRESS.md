# Chandra primary search progress

**All event products have been scanned in 3,680 of 6,608 dated ObsIDs. The accumulated eligible start time for 60-second flares is 59.174125 Ms.**

The deposited files contain 3,680 authenticated observation ledgers. Catalog checks and required photometry finished in 3,468 ObsIDs. The tally spans the immutable checkpoint and 108 of 192 primary worker artifacts. The latest included source material is dated 2026-11-22T10:15:31Z.

The completed scans produced 23,579 real trigger rows. 221 passed the instrumental tests before counterpart review and cross-product deduplication. The observation table keeps trigger rows from partial scans in separate columns.

| Trial duration (s) | Eligible flare-start live time (s) |
|---:|---:|
| 60 | 59,174,125.5 |
| 180 | 57,868,354.1 |
| 600 | 53,311,859.7 |
| 1800 | 41,124,804.9 |
| 3600 | 25,326,451.7 |

Support is the union of available start intervals within each ObsID with its live-time fraction, summed across ObsIDs. Counts and support above require scans of every acquired product, including the 24 scrambled streams and the 12 scheduled injection draws per product.

2,928 dated ObsIDs remain outside this exposure. Which fields appear in a checkpoint depends on worker scheduling, run time and archive or catalog availability. These totals describe the deposited fields. The 14 supplemental ObsIDs with missing release dates are reported separately.

[Observation status and trigger counts](progress_obsids.csv.gz) · [Input URLs and checksums](progress_inputs.csv.gz) · [Machine summary](progress.json) · [Artifact coverage and gaps](coverage.json)

[Prespecified 128-ObsID cohort](../../provisional/moderate-count-v1/44985879443-1/selection.json) · [Separate 14-ObsID supplement](../../supplemental_search/44985879443-1/README.md) · [Original Actions run](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443)
