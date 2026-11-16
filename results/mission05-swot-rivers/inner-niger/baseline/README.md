# Historical comparison for the July–January stage contrast

Comparable independent observations are too sparse to show that the 2024 contrast was unusual for these two river reaches. The complete results include each site's observed dates, levels, reported uncertainties and coverage failures. This is a retrospective comparison with the two calendar dates selected from the SWOT 2024 event.

The DGFI-TUM DAHITI public Mali index contains 155 targets. We used only public targets catalogued as “Niger, River” within 60 km of the two SWOT reach centroids, yielding four near Mopti and five near Diré. A station must span 2018–2026 to permit comparisons at the same crossing before SWOT launch and during the target event. Two near-Mopti and three near-Diré stations met that requirement. Primary stations were selected by geographic distance before inspecting seasonal differences. This is a regional comparison across distinct satellite crossings. Straight-line distance does not guarantee identical channel morphology or flood propagation.

The primary pair was Sentinel-3B target 17277 near Mopti (25.9 km) and Sentinel-3B target 15109 near Diré (22.1 km). **Zero of eight candidate seasons had all four admissible observations.** The nominal dates were July 22 of the flood year and January 26 of the following year, with a 15-day tolerance and reported uncertainty no greater than 0.30 m. A matched upstream stage required an absolute Mopti change no greater than 0.20 m.

The alternative downstream Sentinel-3A target 13763 (43.9 km from the Diré reach) produced the following comparisons with Mopti target 17277.

| Sample | 2024 downstream change | 2024 upstream change | Historical denominator and rank |
|---|---:|---:|---|
| Calendar dates alone | +5.010 m | +0.218 m | Two complete pre-SWOT seasons, 2019 and 2020. 2024 ranks first of three including the event. |
| Upstream match within 0.20 m | — | 0.018 m outside tolerance | One historical match in 2019. The target event is ineligible. |
| Upstream match within 0.50 m | +5.010 m | +0.218 m | One historical match in 2019 (+3.929 m downstream). 2024 ranks first of two; empirical upper-tail rank fraction 0.50. |
| Downstream station alone | +5.010 m | — | Four complete pre-SWOT seasons (2018–2021); range −0.293 to +3.993 m. 2024 ranks first of five; empirical upper-tail rank fraction 0.20. |

The two 2024 endpoints near Mopti were July 25 and January 30. Those near Diré were July 19 and January 24. The downstream-minus-upstream change was +4.791 m. The root-sum-square of the two provider-reported standard errors was 0.074 m upstream and 0.021 m downstream, assuming independent errors at the two endpoints. Those errors exclude spatial representativeness and the six-day offsets between satellite crossings. A tighter 0.15 m uncertainty limit and a seven-day maximum difference between crossings preserve the 2024 target pair and the lone 2019 comparison under the relaxed stage tolerance. For target 13763 alone, 2022–23 returned +4.615 m. January 2023 occurred after SWOT launch.

A third near-Diré target, 1610 (54.1 km from the reach, upstream of Diré), returned just +0.084 m during the same 2024 windows despite having four-date coverage with target 17277. The 0.20 m upstream-stage condition fails for both downstream alternatives. The spread across target locations limits any inference about a downstream regional anomaly. The earlier Envisat/SARAL tracks yielded seven complete station-pair seasons across five distinct years in 2002–2015, with zero qualifying 0.20 m upstream matches; they cannot resolve the 2024 rank either. Sparse samples, unusual proxy readings and retrospective choice of event prevent these empirical ranks from serving as formal discovery significance.

The product is the DAHITI water-level chart from the Technical University of Munich (DGFI-TUM). The two arrays plotted on each station page supply UTC timestamps, orthometric heights and reported errors. The manifest preserves the source URLs, response dates, response SHA-256 digests, numeric CSV digests, stations, missions and tracks. Use `python compare.py` and `python historical.py` to reproduce the tables offline. `python retrieve.py` refreshes the station data. Station IDs and the chosen cutoff rules are in `design.md`.

[Mali station list](https://dahiti.dgfi.tum.de/en/virtual-stations/africa/mali/) · [17277 near Mopti](https://dahiti.dgfi.tum.de/en/17277/water-level-altimetry/) · [15109 near Diré](https://dahiti.dgfi.tum.de/en/15109/water-level-altimetry/) · [13763 downstream comparison](https://dahiti.dgfi.tum.de/en/13763/water-level-altimetry/) · [1610 upstream of Diré](https://dahiti.dgfi.tum.de/en/1610/water-level-altimetry/).
