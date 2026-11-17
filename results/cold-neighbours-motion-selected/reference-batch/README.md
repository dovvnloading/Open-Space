# The held-out single-exposure test loses sensitivity at the faint end

**Eleven of 29 wide-sample positions that coincided with published sources passed the independent NEOWISE motion criterion. None of the 11 at W2 = 15.0–15.5 passed.** The original 29 are positional literature alerts, so these fractions describe this reference cohort and do not by themselves measure purity or completeness for a population of brown dwarfs. Their measured W2 range matters when interpreting a missing single-exposure track.

| CatWISE W2 (Vega mag) | Literature-alert positions | Pass the held-out motion test |
|---|---:|---:|
| Brighter than 14.0 | 3 | 3 |
| 14.0 to <14.5 | 9 | 6 |
| 14.5 to <15.0 | 6 | 2 |
| 15.0 through 15.5 | 11 | 0 |

All 29 positions already passed the wide CatWISE color, motion and brightness selection. The tests use the same per-exposure quality requirements, ≥8 matches, ≥4 held-out visits and comparison with four alternative velocity vectors as the [wide unmatched sample](../widebatch/README.md). The counts reflect proper motion, crowding and survey geometry as well as brightness. The separately identifiable T8 positive control is one row in this reference cohort and one explicitly labelled comparator in the wide packet. None of the 29 contributes to the 988 literature-unmatched denominator.

The [per-source results](summary.csv) give the CatWISE coordinates' identifiers, literature matches/aliases, W2 magnitudes, retained exposure counts and alternative counts. [Inputs](inputs.json) records archive query URLs, raw and stable science-row hashes and the per-track match statistics. From this directory run `python controls.py --selected ../wide/selected.csv --matches ../catalogues/wide_matches.csv` with the sibling [NEOWISE program](../widebatch/screen.py). The [brightness bins](recovery.csv) use the catalogue Vega magnitudes of the selected rows.
