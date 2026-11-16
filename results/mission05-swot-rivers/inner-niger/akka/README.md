# Akka and Diré flood-peak comparison

This follow-up tests whether the unusually high 2024 water level at Diré changed in proportion to the upstream high water at Akka. It was added after selecting the nine reaches in the principal screen and is a separate two-station comparison.

The SWORD v16 Akka Niger reach 14279000061 is a type-1 river reach at the Lac Débo outlet. Its prior centerline passes about 0.39 km from the approximate hydrometric location of 15.39° N, 4.23° W. Among the Niger River polylines examined, it was the closest to that point. The v17b translation is reach 14279000011. The Diré v16 and v17b IDs are 14277700101 and 14277700081. The per-query manifests and satellite CSV granule IDs identify both the data version and the orbit.

OPIDIN publishes annual maximum levels in centimetres above station-specific local datums. The published service distinguishes those historical levels from its display's separate ongoing-season forecasts. The API supplies no per-year observed/reconstructed indicator. The Sahel Rain Index investigators report using Akka only through 2018 and Mopti from 2019. A separate remote-sensing study reports an Akka daily gauge series ending in 2018. The provenance of the post-2018 Akka annual entries is unverified, which precludes their use as independent in-situ validation. [Sahel Rain Index methods](https://sahelrainindex.info/methods-data/) · [Akka gauge study](https://doi.org/10.1080/10106049.2024.2311203).

OPIDIN reports that DNH gauges began at Akka in 1955 and Diré in 1931. Comparisons use changes at each station relative to its own 2023 level. Flood crests can pass the two stations on different dates. Those annual peak differences cannot be used to calculate simultaneous longitudinal water-surface slope.

Run `python peak_matchup.py` to recompute the station and orbit results from the preserved data. A flood season runs from June through the following May. The 2023 vD series begins in July before the rising-season peak. A satellite annual maximum is reported only if at least one valid observation falls in October–January. The number of such observations quantifies peak-period sampling. At an 80% observed-node threshold, neither central reach yields an eligible sampled peak in any of the three seasons. The 75% threshold is a sensitivity calculation. Most of its high-stage returns at Akka have 45/57 nodes and at Diré 46/58 nodes. Retained data have a good or suspect reach flag, a complete reach flag and crossover calibration quality at most one. A further 50% sensitivity returns identical all-orbit peaks. Restricting to pass 363 tests dependence on orbit and swath geometry.

| Flood season | Akka OPIDIN, m local | Akka vD max, m | Diré OPIDIN, m local | Diré vD max, m |
|---|---:|---:|---:|---:|
| 2023 | 4.68 | 263.0078 | 5.02 | 261.6892 |
| 2024 | 5.60 | 264.2396 | 5.91 | 262.5131 |
| 2025 | 5.18 | 263.6487 | 5.28 | 262.1730 |

Both stations rose in 2024. The 2024-minus-2023 changes in sampled maxima were +1.232 m at Akka and +0.824 m at Diré. The OPIDIN series gives +0.92 m and +0.89 m. Restricting to orbit 363 alone gives +1.232 m and +0.860 m. The two series indicate a larger flood and differ by 0.38 m on the change in the difference between the two stations' maxima. Sampling, footprint, and the observations' QC limitations restrict interpretation. The main v3 product has no 2025 flood-season coverage. The Akka satellite height and satellite discharge series use some of the same underlying radar observations.

[Hydrocron discharge](https://podaac.github.io/hydrocron/user-guide/discharge/), [SWORD crosswalk](https://doi.org/10.5281/zenodo.22259077), [OPIDIN hydrometric methods](https://opidin.org/en/background/technical-background-opidin), [historical annual maxima](https://backend.opidin.org/opidin/rest/annualmaxlevels).
