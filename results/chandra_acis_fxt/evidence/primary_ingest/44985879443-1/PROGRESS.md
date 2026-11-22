# Chandra primary search progress

**All event products have been scanned in 1,382 of 6,608 dated ObsIDs. The accumulated eligible start time for 60-second flares is 21.965433 Ms.**

The deposited files contain 1,383 authenticated observation ledgers. Catalog checks and required photometry finished in 1,314 ObsIDs. The tally spans the immutable checkpoint and 41 of 192 primary worker artifacts. The latest included source material is dated 2026-11-22T08:02:14Z.

The completed scans produced 14,621 real trigger rows. 98 passed the instrumental tests before counterpart review and cross-product deduplication. The observation table keeps trigger rows from partial scans in separate columns.

| Trial duration (s) | Eligible flare-start live time (s) |
|---:|---:|
| 60 | 21,965,433.3 |
| 180 | 21,476,368.0 |
| 600 | 19,771,514.2 |
| 1800 | 15,169,758.6 |
| 3600 | 9,248,926.3 |

Support is the union of available start intervals within each ObsID with its live-time fraction, summed across ObsIDs. Counts and support above require scans of every acquired product, including the 24 scrambled streams and the 12 scheduled injection draws per product.

5,226 dated ObsIDs remain outside this exposure. Which fields appear in a checkpoint depends on worker scheduling, run time and archive or catalog availability. These totals describe the deposited fields. The 14 supplemental ObsIDs with missing release dates are reported separately.

[Observation status and trigger counts](progress_obsids.csv.gz) · [Input URLs and checksums](progress_inputs.csv.gz) · [Machine summary](progress.json) · [Artifact coverage and gaps](coverage.json)

[Prespecified 128-ObsID cohort](../../provisional/moderate-count-v1/44985879443-1/selection.json) · [Separate 14-ObsID supplement](../../supplemental_search/44985879443-1/README.md) · [Original Actions run](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443)
