# Chandra primary search progress

**All event products have been scanned in 1,042 of 6,608 dated ObsIDs. The accumulated eligible start time for 60-second flares is 16.351395 Ms.**

The deposited files contain 1,043 authenticated observation ledgers. Catalog checks and required photometry finished in 991 ObsIDs. The tally spans the immutable checkpoint and 31 of 192 primary worker artifacts. The latest included source material is dated 2026-11-22T07:18:02Z.

The completed scans produced 12,108 real trigger rows. 74 passed the instrumental tests before counterpart review and cross-product deduplication. The observation table keeps trigger rows from partial scans in separate columns.

| Trial duration (s) | Eligible flare-start live time (s) |
|---:|---:|
| 60 | 16,351,395.3 |
| 180 | 15,983,226.1 |
| 600 | 14,700,894.0 |
| 1800 | 11,250,486.2 |
| 3600 | 6,805,610.3 |

Support is the union of available start intervals within each ObsID with its live-time fraction, summed across ObsIDs. Counts and support above require scans of every acquired product, including the 24 scrambled streams and the 12 scheduled injection draws per product.

5,566 dated ObsIDs remain outside this exposure. Which fields appear in a checkpoint depends on worker scheduling, run time and archive or catalog availability. These totals describe the deposited fields. The 14 supplemental ObsIDs with missing release dates are reported separately.

[Observation status and trigger counts](progress_obsids.csv.gz) · [Input URLs and checksums](progress_inputs.csv.gz) · [Machine summary](progress.json) · [Artifact coverage and gaps](coverage.json)

[Prespecified 128-ObsID cohort](../../provisional/moderate-count-v1/44985879443-1/selection.json) · [Separate 14-ObsID supplement](../../supplemental_search/44985879443-1/README.md) · [Original Actions run](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443)
