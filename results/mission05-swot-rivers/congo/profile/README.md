# Where the stage response changes along the Kinsuka rapids

We downloaded node time series for every nominal SWORD v16 node in the two focal reaches, 60 in the downstream reach 13217000061 and 46 in the upstream reach 13217000091. Four same-orbit acquisitions supplied 424 records. On pass 111, all 106 fixed nodes had usable water levels on both 8 January and 13 July 2024. On pass 570, 98 met the same criteria on 13 December 2023 and 20 August 2024. We required node quality 0 or 1, valid WSE, and positive WSE uncertainty ≤0.30 m on each date. Actual water-pixel centroids can move within a fixed node. The complete per-node displacement is archived. The median movement was 82 m on pass 111 and 51 m on pass 570.

We divided each reach into consecutive nonoverlapping bins of five SWORD node indices, giving roughly 1 km per bin at nominal 200 m spacing. Across the reach boundary, the contrast between two bin centers separated by ten nominal nodes (roughly 2 km) captures a large part of the reach-scale change. The bins were noticed after inspecting the profile. The complete spatial profile and second-orbit check show the effect of that selection.

| Five-node bin in downstream-to-upstream order | Low-water location near bin center | Jan minus Jul 2024, median rise | Dec 2023 minus Aug 2024, median rise |
| --- | --- | ---: | ---: |
| Downstream reach, nodes 51–55 | 15.216° E, 4.327° S | 13.52 m | 11.93 m |
| Downstream reach, nodes 56–60 | 15.224° E, 4.322° S | 9.24 m | 7.59 m |
| Upstream reach, nodes 1–5 | 15.232° E, 4.319° S | 5.46 m | 4.58 m |

The two outer bins give contractions of 8.06 m and 7.35 m on the two tracks. Each of the six table entries uses five nodes. With a 200 m limit on high-to-low centroid movement, the contractions become 8.15 m and 7.35 m (4/5 and 3/5 nodes in the outer downstream/upstream bins). A 100 m limit leaves just 1/3 and 2/2 nodes and gives 6.95 m and 7.48 m. A 50 m limit has no retained outer downstream nodes. The adjacent bins share water, calibration and geolocation conditions. The records identify a narrow zone of stage amplification and provide no independent sample of 106 floods. This is a water-surface measurement. Flow partitioning among channels and individual rapid crests require a field survey for attribution.

From the mission directory run `python congo/profile/profile.py`. To regenerate the frozen measurements from source run `python congo/profile/profile.py --retrieve /path/to/node-cache`. `manifest.csv` retains all 106 NASA Hydrocron URLs and hashes for both JSON and decoded CSV. The extracted rows preserve dates, orbit, cycle, CRID, node coordinates, uncertainties, quality masks and source granule IDs. [Four-visit extracts](measurements.csv) · [All fixed-node responses](profile.csv) · [Bins](bins.csv)
