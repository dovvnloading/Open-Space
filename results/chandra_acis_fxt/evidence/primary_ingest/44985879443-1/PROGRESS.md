# Chandra primary search progress

**All event products have been scanned in 2,002 of 6,608 dated ObsIDs. The accumulated eligible start time for 60-second flares is 32.038457 Ms.**

The deposited files contain 2,003 authenticated observation ledgers. Catalog checks and required photometry finished in 1,894 ObsIDs. The tally spans the immutable checkpoint and 59 of 192 primary worker artifacts. The latest included source material is dated 2026-11-22T08:36:27Z.

The completed scans produced 16,867 real trigger rows. 134 passed the instrumental tests before counterpart review and cross-product deduplication. The observation table keeps trigger rows from partial scans in separate columns.

| Trial duration (s) | Eligible flare-start live time (s) |
|---:|---:|
| 60 | 32,038,457.2 |
| 180 | 31,329,205.6 |
| 600 | 28,855,682.2 |
| 1800 | 22,218,662.2 |
| 3600 | 13,680,570.2 |

Support is the union of available start intervals within each ObsID with its live-time fraction, summed across ObsIDs. Counts and support above require scans of every acquired product, including the 24 scrambled streams and the 12 scheduled injection draws per product.

4,606 dated ObsIDs remain outside this exposure. Which fields appear in a checkpoint depends on worker scheduling, run time and archive or catalog availability. These totals describe the deposited fields. The 14 supplemental ObsIDs with missing release dates are reported separately.

[Observation status and trigger counts](progress_obsids.csv.gz) · [Input URLs and checksums](progress_inputs.csv.gz) · [Machine summary](progress.json) · [Artifact coverage and gaps](coverage.json)

[Prespecified 128-ObsID cohort](../../provisional/moderate-count-v1/44985879443-1/selection.json) · [Separate 14-ObsID supplement](../../supplemental_search/44985879443-1/README.md) · [Original Actions run](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443)
