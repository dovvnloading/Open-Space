# Chandra ACIS fast X-ray transient search

**No FXT survived in 3 scanned ObsIDs covering 0.034828 Ms of eligible 60-second flare-start live time.**

The frozen release inventory contains 6,608 ObsIDs. All data products were scanned in 3. Catalog calls and needed photometry completed in 3. The release-date or mode ledger contains 14 exceptions.

Across 24 full-survey scrambles the mean number passing the same instrumental and counterpart screen was 0.000. The real fields produced 0 raw triggers. Of 30 executed injections, 24 passed instrumental selection. 6 additional draws lacked usable observing time or geometry.

| Trial duration (s) | Eligible field-start live seconds |
|---:|---:|
| 60 | 34,827.7 |
| 180 | 33,772.7 |
| 600 | 30,080.0 |
| 1800 | 19,529.6 |
| 3600 | 4,727.3 |

The one-sided 95% Poisson upper bound is 86 detected events per Ms for this field mixture, using 60-second start support. Intrinsic sky rates require effective survey area and a source population model.

[Input checksums](event_receipts.csv) · [Coverage](obsid_coverage.csv) · [Scrambles](control_trials.csv) · [Injection grid](completeness.csv) · [All unresolved leads with positions, durations, fluxes, hosts and input hashes](leads.json)

Each compressed ObsID ledger records all raw triggers, vetoes, controls, catalog receipts and injection outcomes. Leads require examination of the target and proposal, independent visits, the detector image, and contemporary notices before a claim of novelty.

[Methods and reproduction](../../../METHODS.md) · [Actions run](https://github.com/dovvnloading/Open-Space/actions/runs/44984928643)

This pre-fan-out checkpoint measures only the three pilot observations selected by catalog event count before the scan. Every other ObsID, including pending observations and recovered shards, is outside this measurement.
