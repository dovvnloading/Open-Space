# Later-date test of the Congo stage response

We fitted the Theil–Sen line to the 17 quality-screened pass-111 observations with simultaneous Chenal stage in the original RiverSP 2.0 analysis. We froze the direction and fitted relationship before downloading the later RiverSP D time series. We used the median-based intercept implemented by SciPy and the mapped lower, upper and Chenal reach IDs from the official SWORD v17b-to-v16 translation. All three mappings have zero boundary-change flag. The lower reach has a slightly different centroid in D. The December 2023–January 2024 flood observations are part of the fit in this particular test. The separate flood-window test in the parent directory omits them.

Seven new visits from 14 August 2025 through 24 August 2026 passed all three reach filters after the v3 time window ended on 2 May 2025. Spearman ρ was −0.964. Their median absolute difference from the line fitted to earlier visits was 0.638 m. Chenal stages spanned 274.17–276.18 m, within the original 273.53–279.93 m range. Only these seven of the 19 later paired outlet visits also supplied the upstream stage under the filters. The holdout therefore tests recurrence at lower and intermediate stages. It contains no second flood of the January 2024 magnitude. On seven same-period D visits the correlation was −1 and median absolute prediction error was 0.44 m. RiverSP D uses the same satellite and places each observation in a single processing release. It provides a temporal holdout and processing-sensitivity check. All discharge values in the mission’s primary analysis come from SoS v3.

## Reprocessing the four original visits

The table gives upper minus lower reach elevation, before quality filtering. The dates and orbits are the same in each release. All four pairs pass the study criteria in RiverSP 2.0. In RiverSP D only the July pair passes.

| Orbit and dates | RiverSP 2.0, high / low | RiverSP D, high / low | RiverSP D quality |
| --- | ---: | ---: | --- |
| 111, Jan / Jul 2024 | 8.7872 / 19.5286 m | 8.8374 / 19.5628 m | Upper reach flagged degraded in January. |
| 570, Dec 2023 / Aug 2024 | 9.8550 / 19.5243 m | 22.3078 / 19.3858 m | Upper reach covered 78.26% of nodes on each date, below the 80% cutoff. |

The raw pass-111 contractions agree within 0.02 m. On pass 570, RiverSP 2.0 gives a 9.67 m contraction and RiverSP D a 2.92 m increase in reach-scale fall. Neither orbit has a quality-accepted high/low reach comparison in D. The seven later visits tested the lower-stage relationship on pass 111, a different question from reproducing the flood measurements. [Same-acquisition heights, flags, coverage and granules](event_by_release.csv) · [Fixed-node comparison](../profile/README.md)

From the mission directory run `python congo/holdout/test_d.py`. The default verifies the frozen D CSVs and regenerates the result tables. The optional `--retrieve` refreshes the three queries into the current checkout. Retrieval date for the archived files is 16 November 2026. [Scores](summary.csv) · [Visits](visits.csv) · [Source receipt](manifest.json) · [SWORD ID translation](https://doi.org/10.5281/zenodo.22259077)
