# Where the stage response changes along the Kinsuka rapids

We downloaded node time series for every nominal SWORD v16 node in the two focal reaches, 60 in the downstream reach 13217000061 and 46 in the upstream reach 13217000091. Four RiverSP 2.0 acquisitions on two tracks supplied 424 records. On pass 111, all 106 fixed nodes had usable water levels on both 8 January and 13 July 2024. On pass 570, 98 met the same criteria on 13 December 2023 and 20 August 2024. We required node quality 0 or 1, valid WSE, and positive WSE uncertainty ≤0.30 m on each date. Actual water-pixel centroids can move within a fixed node. The complete per-node displacement is archived. The median movement was 82 m on pass 111 and 51 m on pass 570.

We divided each reach into consecutive nonoverlapping bins of five SWORD node indices, giving roughly 1 km per bin at nominal 200 m spacing. Across the reach boundary, the contrast between two bin centers separated by ten nominal nodes (roughly 2 km) captures a large part of the reach-scale change. The bins were noticed after inspecting the profile. The complete spatial profile and second-orbit check show the effect of that selection.

| Five-node bin in downstream-to-upstream order | Low-water location near bin center | Jan minus Jul 2024, median rise | Dec 2023 minus Aug 2024, median rise |
| --- | --- | ---: | ---: |
| Downstream reach, nodes 51–55 | 15.216° E, 4.327° S | 13.52 m | 11.93 m |
| Downstream reach, nodes 56–60 | 15.224° E, 4.322° S | 9.24 m | 7.59 m |
| Upstream reach, nodes 1–5 | 15.232° E, 4.319° S | 5.46 m | 4.58 m |

The two outer bins give contractions of 8.06 m and 7.35 m on the two tracks. Each of the six table entries uses five nodes. With a 200 m limit on high-to-low centroid movement, the contractions become 8.15 m and 7.35 m (4/5 and 3/5 nodes in the outer downstream/upstream bins). A 100 m limit leaves just 1/3 and 2/2 nodes and gives 6.95 m and 7.48 m. A 50 m limit has no retained outer downstream nodes. The adjacent bins share water, calibration and geolocation conditions. The records identify a narrow zone of stage amplification and provide no independent sample of 106 floods. This is a water-surface measurement. Flow partitioning among channels and individual rapid crests require a field survey for attribution.

## Same acquisitions in RiverSP D

We rechecked the three five-node bins on all four original acquisitions in release D. The official v17b–v16 node translation maps each of these 15 nodes directly, with zero shift and boundary flags. The comparison uses the same per-node quality and uncertainty cuts as the original analysis. The high- and low-water labels describe the independently documented hydrological periods.

| Track and release | Downstream median rise, nodes 51–55 | Middle median rise, nodes 56–60 | Upstream median rise, nodes 1–5 |
| --- | ---: | ---: | ---: |
| 111, D, January minus July | +13.57 m (5 pairs) | 0 pairs | 0 pairs |
| 570, D, December minus August | −10.15 m (3 pairs) | −14.68 m (4 pairs) | −17.77 m (4 pairs) |

On pass 111, release D reproduces the original 13.52 m rise in the downstream bin. All January middle and upstream nodes fail the D screen, leaving the 8.06 m local contraction unconfirmed in D. On pass 570, the eleven nodes accepted in both releases differed between releases by −22.79 to −21.45 m in December and −0.04 to +0.18 m in August (D minus 2.0). Consequently the absolute seasonal elevations reverse direction. Their spatial contrast survives much of this offset. Across the three downstream and four upstream nodes common to both releases, the contrast is 7.11 m in 2.0 and 7.62 m in D. At reach scale the two releases sample different portions of the reaches, and the unscreened D contrast changes sign. These acquisitions supply a processing-sensitive check of localization. They do not independently validate the seasonal elevations in all three bins. [Reach comparisons](../holdout/event_by_release.csv) · [Node counts and medians](d-node-check.csv) · [Common-node measurements](d-node-common.csv)

Two additional mapped upper-reach nodes had same-day ICESat-2 observations on 13 December. At upstream node 35 the laser median was 272.057 m (59 segments within 0.5 km of the 2.0 centroid). RiverSP 2.0 measured 272.106 m and D measured 272.146 m with an accepted node flag of 1. The centroid moved 28 m between releases. At node 12 the equivalent laser median was 270.590 m (86 segments), compared with 270.471 m in 2.0. D reported 270.506 m with node flag 3 and a 163 m centroid shift. These upstream footprints do not coincide with the three bins. Their close agreement cannot resolve the roughly 22 m December discrepancy near the reach boundary. [Laser and D matchups, both radii](d-laser-check.csv) · [Original laser matching method](../../laser/README.md)

From the mission directory run `python congo/profile/profile.py`. To regenerate the frozen measurements from source run `python congo/profile/profile.py --retrieve /path/to/node-cache`. `manifest.csv` retains all 106 NASA Hydrocron URLs and hashes for both JSON and decoded CSV. The extracted rows preserve dates, orbit, cycle, CRID, node coordinates, uncertainties, quality masks and source granule IDs. [Four-visit extracts](measurements.csv) · [All fixed-node responses](profile.csv) · [Bins](bins.csv)

Run `python congo/profile/compare_d.py` from the mission directory to verify and replay the D comparison. The `--retrieve` option refreshes the queries. The compressed payload archive retains the original per-query CSV text. [SWORD mapping, source member and hash](d-node-translation.json) · [NASA queries and payload hashes](d-node-receipts.json) · [Archived D measurements](d-node-payloads.json.gz)
