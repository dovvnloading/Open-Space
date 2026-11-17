# ATLAS shutter coverage of isolated TESS boxes

No public ATLAS exposure in the four queried CATCH holdings overlapped any of the 252 isolated TESS boxes on 37 stars. The closest shutter interval ended or began 33.5 minutes from a box edge. The result supplies no simultaneous ground-based flux constraint on these events.

| TESS sample | Boxes checked | Boxes with a footprint in their padded query | Closest shutter to a box edge |
|---|---:|---:|---:|
| 39-star FFI extension | 85 | 65 | 131.7 min |
| Earlier auxiliary FFI visits | 24 | 24 | 33.5 min |
| 20/120-s deep-dip screen | 142 | 114 | 156.1 min |
| TIC 55997020, partial Sector 110 | 1 | 1 | 448.3 min |

A box is one detector window. Simultaneous boxes from different widths, cadences or stars remain separate. The 59 queries returned 340 target-footprint rows from 336 distinct exposures. Haleakalā contributed 121 rows, Maunaloa 67, Río Hurtado 56 and Sutherland 96. There were 38 nonempty queries. A missing exposure can reflect the source's observational cadence or the archive's holdings. The independent tests of predicted periodic returns on other dates are in [ATLAS_EPOCHS.md](ATLAS_EPOCHS.md).

Each request uses the white dwarf's Gaia position propagated to the date group's midpoint. It covers at least one complete UTC day beyond either TESS box edge. The audit verifies that all four documented ATLAS source keys are in the reply and that CATCH returned the declared number of products. For every returned observation, the reported exposure time agrees with its shutter stop minus start. The code converts both UTC shutter endpoints to barycentric TDB using the documented site longitude, latitude and height, the Gaia direction propagated to each endpoint, and Astropy's solar-system ephemeris. It tests for strictly positive overlap with the TESS center ± half the reported duration. The direct Sector 110 measurement uses the archived FFI field timing. The nearest returned shutter for this box is 7.47 hours away. [CATCH API](https://catch-api.astro.umd.edu/openapi.json) · [ATLAS station coordinates, §3](https://pds-smallbodies.astro.umd.edu/data_sb/missions/atlas/ATLAS_Operations_and_Data_Processing.pdf)

The response bundle preserves all 59 original CATCH replies. The receipt gives each URL, retrieval time, byte count and SHA-256. The window ledger links each box to its query and records the nearest shutter gap in that query when one was returned. [Query receipts and counts](ground_windows/coverage.json) · [Frozen responses](ground_windows/responses.json.gz) · [Window ledger](ground_windows/windows.json.gz) · [Timing and extraction code](atlas_window_coverage.py)

```bash
python atlas_window_coverage.py short_cadence/roster_screen.json.gz atlas-cache atlas-audit --frozen ground_windows
```

Omit `--frozen ground_windows` and select a fresh cache directory for a new archive search. Overlapping images, when present in future releases, trigger same-site, same-filter, same-night cutout photometry.
