# Nearest-neighbor crater diameters in modern Martian impact fields

Across 259 independently cataloged Martian impact fields containing 13,034 mapped craters at least 1 m across, the crater nearest a given crater has a more similar diameter than expected after accounting for broad size sorting along the field. Mean absolute difference in diameter-percentile rank was **32.57 percentile points**, against **35.57** when sizes are exchanged within bands along the field's long axis. The site-weighted reduction was **3.00 points** (95% field-bootstrap interval 2.29–3.69), or **8.4%** of the expected gap. The reduction occurred in 191 of 259 fields. The fourth-nearest crater shows no detectable residual similarity.

This is a new analysis of existing measurements of individual craters. It measures an ensemble spatial association in the published catalog. Correlated fragment sizes from related breakup events, spatially varying crater formation conditions, and correlated mapping or visibility errors can all produce local size similarity. The measurements do not uniquely identify a fragmentation mechanism.

The [comparison with a published atmospheric-breakup simulation](MODEL.md) finds that it reproduces part of the nearest-neighbor effect. The observed and simulated fields differ most clearly at the two smallest measured diameter thresholds.

## Measurement and countertests

One impact field contributes one measurement. We retained fields with at least 12 craters above the selected size cutoff and ranked diameter within each field, assigning ties their mean rank, then rescaled ranks to 0–1. For every crater we measured the absolute rank difference to its nearest center. Mutual nearest neighbors contribute both directed edges. We computed uncertainty across impact fields so these edges do not count as independent impacts.

The primary comparison permutes diameters within up to four equal-count positional bands on the first principal axis of crater centers, subject to at least four craters per band. The null preserves all crater positions, the field's diameter distribution, the neighbor graph, and broad longitudinal size trends. It does not model smooth variation within each band. Its expected absolute difference is calculated exactly for each edge, including the without-replacement correction when both ends share a band. In 9,999 independent realizations across all fields, no site-weighted mean was as low as the observation (one-sided add-one permutation p = 0.0001). The central 95% range of simulated means was 34.92–36.26 percentile points. This probability is conditional on the selected analysis. The investigation is exploratory and was not preregistered.

| Selection and control | Fields | Observed mean gap (points) | Expected (points) | Reduction (points), 95% field-bootstrap interval |
|---|---:|---:|---:|---:|
| ≥1 m, unrestricted in-field shuffling | 259 | 32.57 | 36.66 | 4.09 [3.30, 4.87] |
| ≥1 m, along-axis comparison | 259 | 32.57 | 35.57 | 3.00 [2.29, 3.69] |
| Along-axis and cross-axis bands | 259 | 32.57 | 35.01 | 2.44 [1.85, 3.02] |
| Along-axis and local-density bands | 259 | 32.57 | 34.71 | 2.13 [1.51, 2.76] |
| Isolated endpoints and isolated exchange pool | 214 | 33.24 | 34.89 | 1.65 [0.73, 2.55] |
| ≥2 m, along-axis comparison | 99 | 33.79 | 36.08 | 2.29 [1.20, 3.40] |
| ≥3 m, along-axis comparison | 37 | 34.01 | 35.68 | 1.67 [−0.21, 3.62] |
| Second-nearest neighbor, ≥1 m | 259 | 34.64 | 35.71 | 1.08 [0.43, 1.71] |
| Fourth-nearest neighbor, ≥1 m | 259 | 36.43 | 36.13 | −0.30 [−0.90, 0.30] |

The isolated subset uses 4,280 directed edges for which **both endpoints** are farther than the field's largest retained diameter from every other retained center. Diameters are exchanged only among equally isolated positions in their own along-axis band. No reassigned retained rim can overlap another retained rim at an isolated position. We require six usable edges per field in the restricted tests. This weakens the signal, and spatial variations in terrain and image contrast remain possible.

The two-dimensional control bisects each along-track band by cross-track position. The density control instead bisects each band by distance to the fourth-nearest crater. The reductions are 3.08 [2.03, 4.10] points at MRO orbit numbers below 40,000 (117 fields) and 2.92 [1.90, 3.92] points at orbit 40,000 or later (142 fields). These are disjoint field subsets sharing a mapping pipeline. The published analysis reports 0.25 m/pixel imagery for most fields and 0.5 m/pixel in twelve, so the 1 m cut approaches the resolution limit. At 2 m the signal persists in 99 fields, and the 37-field result at 3 m is inconclusive. The 10,000-draw field bootstrap treats fields as exchangeable. Survey targeting and regional correlations limit extrapolation to every Martian impact site.

## Prior work and provenance

Daubar and colleagues explicitly used the placement of the largest craters to infer entry direction in 77 fresh impact fields. Neidhart and colleagues mapped 557 additional fields and analyzed dispersion, diameter distributions, breakup and elevation. Grant and colleagues measured a diameter gradient in the older Jau cluster. A 2026 study compared 89 Martian two-crater fields with lateral fragment-spreading models. These motivate conditioning on the broad size gradient. In a focused literature search on 2026-11-15, we found no prior report quantifying the within-field nearest-neighbor diameter difference across these multi-crater fields under this conditional comparison. Exhaustive priority and an independent repeat mapping of crater diameters remain open. [Daubar et al. 2019](https://doi.org/10.1029/2018JE005857) · [Neidhart et al. 2023](https://doi.org/10.1029/2022JE007611) · [Grant et al. 2025](https://doi.org/10.1029/2024JE008763) · [Lateral spreading of meteoroids fragments, 2026](https://doi.org/10.1016/j.icarus.2026.117288)

The input is Neidhart, Sansom, Miljković, Daubar, Eschenfelder & Collins, *New crater clusters on Mars* (2021, CC BY 4.0), including the earlier measurements cited there. The published ZIP contains 634 per-field Excel workbooks and a field summary. We reproduced the published counts in all 634 and read 15,091 mapped diameters. Thirteen were below 1 m. Geometry was checked independently by computing each field's median all-pairs separation from its member positions on a spherical Mars with radius 3,389.5 km. All but one agreed within 0.001%. In ESP_012047_1755 we recovered 39.58 m against the tabulated 45.77 m. Omitting this field increases the primary difference from 3.00 to 3.02 points. Fields and radii here come from the published mapping, without independent image remeasurement. [Source archive](https://doi.org/10.5281/zenodo.5760176) · [Source license](https://creativecommons.org/licenses/by/4.0/)

A [direct HiRISE image check](pixel-check/README.md) in field ESP_066902_1830 locates 71 of 76 dark centers. A diameter-blind pixel-width estimate correlates with catalogued size (Spearman ρ = 0.563), above all 144 translated controls. This checks one field on the original exposure.

## Reproduce

From this directory, with Python 3, NumPy, SciPy and openpyxl installed, run:

```sh
curl -fL 'https://zenodo.org/api/records/5760176/files/cluster_data.zip/content' -o cluster_data.zip
python analyze.py cluster_data.zip rerun
```

The program verifies SHA-256 `f32fd7ece38d8ebb9a8874651d14f09c450278f8e58a27858c12f88219651aed` and reads diameter and coordinates by header name. The workbooks use different local longitude origins, so only within-field differences enter the calculation. The outputs include every eligible field's observed and expected score for every sensitivity test, the full all-field geometry audit, and summary statistics. In the CSV, “n” counts all above-cutoff craters in the field and “edges” counts the directed edges used. In the JSON, gap statistics use a 0–1 rank scale. Multiply by 100 for the table above.
