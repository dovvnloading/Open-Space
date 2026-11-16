# Deep-water Florida survey

We sampled 172 dates between March 5, 2024 and March 10, 2026, taking every fourth available date from the first 687 time indices in the USF IMaRS daily OCI grid. The lattice spacing is about 0.216° in each direction. The screen uses ten reflectance bands between 413 and 678 nm. The decline and the rebound in the May 2025 Atlantic event have additional spectra for all 149 bands published here, spanning 346–688 nm with gaps. The ultraviolet bands are archived for provenance. Only visible bands enter the metric. [USF product and metadata](https://erddap.marine.usf.edu/erddap/griddap/pace_oc_aop_1d_fl.html)

The observing unit is a fixed grid point with valid reflectance on both dates. Bounds are in degrees, inclusive. The depth column counts lattice points whose nearest sampled ETOPO1 elevation is at least 200 m below sea level. Shallow bottom interference is less likely there. Currents can still transport coastal water into these regions. [NOAA ETOPO1](https://upwell.pfeg.noaa.gov/erddap/griddap/etopo180.html)

| Region | Latitude °N | Longitude °E | Deep points | Passing pairs / 171 | Valid point-visits |
| --- | --- | --- | ---: | ---: | ---: |
| Eastern Gulf | 24.1–29.9 | −88.95 to −86.0 | 348 | 2 | 6,215 |
| West Florida slope | 24.8–29.9 | −85.99 to −82.5 | 112 | 3 | 2,149 |
| Florida Straits | 24.1–25.3 | −82.49 to −79.15 | 48 | 1 | 868 |
| Florida Atlantic | 25.31–29.9 | −80.50 to −79.15 | 84 | 3 | 1,487 |

USF identifies the source as NASA Ocean Biology Processing Group Level 2. Its regional composites report V3.1 and the ATMFAIL, LAND, CLDICE and HILT exclusions. We retain samples with zero stray-light mask, positive Rrs442, Rrs490 and Rrs560 below 0.05 sr⁻¹, Rrs560 > 0.00005 sr⁻¹, 0 ≤ Rrs665 < 0.0015 sr⁻¹ and Rrs665/Rrs560 < 0.5. The red thresholds are exploratory quality filters. Missing data, red brightness and atmospheric conditions affect coverage. These filters cannot establish cloud-free conditions in every accepted pixel.

For each pair at most eight days apart we require at least ten common points and 15% of the deep lattice in its box. Of 684 regional date pairs, nine pass. The curvature is ln(Rrs490) − [70 ln(Rrs442) + 48 ln(Rrs560)]/118. The largest ranked decline and rebound were May 13–17 and 17–21, 2025, respectively. Both signs were coherent at all 18 and 13 OCI points. The first footprint spans 27.16–29.74°N, with a median seabed depth of 575 m. Sparse point spacing leaves the intervening area unmeasured. The Sentinel-3 comparison uses the same locations on both dates. The May 2024 Gulf cruise matchups in the Sentinel-3 folder are a separate instrument check and do not validate the May 2025 Atlantic event.

The five V3.2 granules in `matched_l2_catalog.json` intersect the specified Atlantic search box. They identify possible overpasses for GOES timing. The Earthdata download and OPeNDAP returned HTTP 401 without an authenticated session. We did not reproduce the May event using V3.2 swath radiances or individual Level-2 flags. NASA’s V3.2 corrected the bidirectional reflectance correction in earlier processing. Aerosol, viewing geometry, regional processing, and time or footprint differences remain competing explanations for the OCI/OLCI disagreement. [NASA processing notes](https://nasa.github.io/oceandata-notebooks/notebooks/release_note_highlights.html)

## Reproduction

From this folder run the following commands. Python dependencies are NumPy, SciPy and pandas. `source_manifest.csv` and `event_receipts.csv` record the analyzed responses and their SHA-256 digests. ERDDAP can include the time of the request in its NetCDF metadata, so compare decoded arrays when a renewed download changes the response digest.

```sh
python fetch.py
python analyze.py
python footprints.py
python event_spectra.py
```

`regional_visits.csv` records date coverage and medians. `regional_changes.csv` records all tested intervals, including those failing coverage. `footprints.csv` records endpoint coordinates, depths, and reflectances for the nine accepted comparisons. `event_spectra.csv` gives separate endpoint medians for each May footprint. `usf_info.csv` preserves the provider's band inventory and processing metadata.
