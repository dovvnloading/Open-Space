# Frozen IceCube track inventory

The historical search has 340 eligible spatial localizations. The 2024–2026 search has 23 native localizations among 69 active Gold/Bronze alerts. The observation cutoff is 2026-11-17 12:09:02 UTC. Inputs were retrieved between 12:08 and 12:26 UTC on that date. Detection times and circular publication times are cut at the observation cutoff. The last included event occurred November 13.

| Input | Frozen coverage | Geometry used |
|---|---:|---|
| IceCat-1, Harvard Dataverse V4.0 | 348 distinct run/event IDs, May 2011–October 2023, including 8 cosmic-ray veto flags | Every annual FITS archive was read. All 348 masks are in `icecat1_likelihood90.jsonl`. |
| IceCat-2, Harvard Dataverse V1.0 | 7 selected events, all present in V4 | Seven publisher-supplied 90% coordinate loops supersede the corresponding V4 regions in the combined sample. The release is preliminary. |
| GCN/AMON, January 2024–cutoff | 145 notice rows, 71 unique events, 2 retracted | 23 released multi-order maps. The other 46 active alerts have circular error bounds in `gcn_alerts.csv`. |
| VizieR paper Table 2 | 275 entries through 2020 | Initial-publication comparison only. Five repeated date labels in this secondary rendering make name-only deduplication unsafe. |

The 69 active GCN alerts comprise 18 from 2024 (5 Gold), 27 from 2025 (8 Gold) and 24 from 2026 (6 Gold). The two excluded track notices are 139977_2910365 (GCN 37794) and 141456_888784 (GCN 42172). The circular audit also records two cascade retractions, which have no matching track run/event. The 240105A error correction from circular 35498 supersedes 35485. The GCN joins use detector IDs, explicit retraction circulars and the latest refined position. Date-only fallback requires a single event on the date. Duplicate submissions containing no reconstruction cannot replace a refined position.

## Spatial footprints

`archival_maps.jsonl.gz` is the science input for the historical sample. Rows carry UTC, stream, source identification by run/event and geometry. The file drops the eight V4 cosmic-ray vetoes and joins the seven updated reconstructions by run/event. It contains 333 lossless nested HEALPix covers and seven publisher polygons. Each polygon is one closed ICRS loop with no missing vertices. Spherical triangle integration agrees with each publisher area at its 0.01 deg² reporting precision.

For IceCat-1 we used finite HEALPix cells satisfying the event FITS header's 90% delta-2LLH criterion. The criterion is 64.2 for 343 events. Five maps specify Wilks' theorem, which gives 4.60517 for the two fitted angular coordinates. The nearest-pixel likelihood minimum is zero in every map. Complete NESTED sibling quartets were recursively merged, with equality of the footprint size checked per map. `area_deg2` sums the original equal-area cells. The official V4 table disambiguates all event names.

For recent GCN maps we sorted cells by probability density per steradian and accumulated density times solid angle. The FITS headers specify highest posterior density and instruct users to ignore non-finite values. The JSON contains every selected NUNIQ cell, including the cutoff cell and equal-density ties. Native finite integrated density ranges from 0.978338 to 1.001875. The calculation normalizes each finite map to its own total. Differences between the resulting area and the separately rounded pipeline area reach 0.1451 deg² in absolute terms and 6.17% in relative terms (different events). Native continuous edges are unresolved at the pixel scale. No exact polygon for the remaining 46 active alerts was found in the published GCN/ROC track archive. The notice radii and circular extrema are retained for audit.

`contours.py` converts NUNIQ cells to merged index ranges for vectorized ICRS membership. It handles RA wrapping and polar geometry with astropy-healpix. The 23 modern best-fit directions and all included MOC cell centers tested inside their respective masks.

## Known-source checks

`native_control_membership.csv` applies the two native catalogue localizations to positions in the sibling historical-control files. TXS 0506+056 is inside both regions. AT2019aalc is inside V4 and outside the revised region. AT2019fdr and AT2019dsg are outside both native regions. The refined dataset and the original alert/counterpart selection procedures support distinct historical questions.

`figure_controls.csv` independently checks the four stars in the original IceCat-2 proceedings figures. It records the per-panel SHA-256, filled areas, star radius and signed pixel distances (positive inside). Dashes were joined at five radii (12, 16, 20, 24, 28 original image pixels). The old AT2019dsg star overlaps the dashed edge at the published image resolution. All four numerical results above rely on the native datasets. The publication-source archive SHA-256 is `49f9b1bf84a3a08e51e914e27fcf53da3b506e42d3d813d2c41723a9caf1e29e`.

## Reproduction and receipts

Run these commands from this directory in a Python environment with NumPy, Astropy, astropy-healpix, SciPy and Pillow. Use a new cache to retrieve the source data. Mutable notice/circular indexes can change on re-fetch; retain the committed derived inputs for exact snapshot replay.

```sh
python freeze.py --cache /tmp/ice-alert-cache
python dataverse.py --cache /tmp/ice-alert-cache
python archival.py
curl -fLsS https://arxiv.org/src/2507.06176v1 -o /tmp/icecat2-source.tar
python figure_controls.py /tmp/icecat2-source.tar
```

`source_receipts.csv` records SHA-256, bytes and source URLs for the GCN notice index, bulk circular archive, ROC map listing, all 23 GCN FITS files, publication Table 2 and two DataCite records. `gcn_circular_audit.csv` identifies relevant individual archive members and their hashes. `dataverse_receipts.csv` records hashes and direct URLs for the version metadata, both tables, the seven-event contour/FITS packages, and all 13 annual V4 TARs (323,431,424 bytes). All annual TARs and both IceCat-2 ZIPs passed the publisher's MD5. The API's ingested TSV representation is hashed as delivered. Its byte representation differs from the pre-ingest original. Each V4 compressed FITS and each IceCat-2 contour text member has an additional digest in its geometry row.

Primary sources: [IceCat-1 V4](https://doi.org/10.7910/DVN/SCRUCD), [IceCat-2 V1](https://doi.org/10.7910/DVN/RX28YT), [IceCat-2 proceedings and control figures](https://arxiv.org/abs/2507.06176), [IceCube GCN mission documentation](https://gcn.nasa.gov/missions/icecube), [GCN archive](https://gcn.gsfc.nasa.gov/amon_icecube_gold_bronze_events.html), [ROC multi-order map archive](https://roc-2.icecube.wisc.edu/public/alerts/).
