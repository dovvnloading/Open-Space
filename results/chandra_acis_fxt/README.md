# Chandra ACIS fast X-ray transient search

**The scan covered 1.591199 Ms of eligible 60-second flare-start live time in 128 ObsIDs and retained 2 unresolved leads.**

The frozen release inventory contains 6,608 ObsIDs. All data products were scanned in 128. Catalog calls and needed photometry completed in 124. The release-date or mode ledger contains 14 exceptions.

Across 24 full-survey scrambles the mean number passing the same instrumental and counterpart screen was 0.167. The real fields produced 18 raw triggers. Of 1,297 executed injections, 939 passed instrumental selection. 239 additional draws lacked usable observing time or geometry.

| Trial duration (s) | Eligible field-start live seconds |
|---:|---:|
| 60 | 1,591,198.7 |
| 180 | 1,546,140.0 |
| 600 | 1,392,053.9 |
| 1800 | 1,022,433.4 |
| 3600 | 572,183.3 |

[Input checksums](evidence/provisional/moderate-count-v1/44985879443-1/event_receipts.csv) · [Coverage](evidence/provisional/moderate-count-v1/44985879443-1/obsid_coverage.csv) · [Scrambles](evidence/provisional/moderate-count-v1/44985879443-1/control_trials.csv) · [Injection grid](evidence/provisional/moderate-count-v1/44985879443-1/completeness.csv) · [All unresolved leads with positions, durations, fluxes, hosts and input hashes](evidence/provisional/moderate-count-v1/44985879443-1/leads.json)

Each compressed ObsID ledger records all raw triggers, vetoes, controls, catalog receipts and injection outcomes. Leads require examination of the target and proposal, independent visits, the detector image, and contemporary notices before a claim of novelty.

[Methods and reproduction](METHODS.md) · [Actions run](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443)

The early cohort selected 128 distinct observations in four catalog event-count strata before opening their event lists. It excludes observations above 200,000 catalog counts and observations without a count estimate. The denominator includes only completed product scans, and all pending observations remain assigned in the full survey.

[Selected ObsIDs and selection rule](evidence/provisional/moderate-count-v1/44985879443-1/selection.json) · [Per-observation ledgers](evidence/observations) · [Field notes](FIELD_NOTES.md)
