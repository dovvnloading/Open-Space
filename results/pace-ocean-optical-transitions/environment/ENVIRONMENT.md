# Physical and atmospheric context

The southern California site is USC SeaPRISM (33.56371°N, 118.11782°W). The Florida comparison follows the individual OCI coordinates in the two May 2025 Atlantic pairs. The occurrence of different optical dates and a different number of usable samples in each pair changes the May 17 spatial medians.

## Florida Atlantic

| Paired sample | Dates in May 2025 | SST anomaly median, °C | Eastward geostrophic velocity median, m/s | Northward geostrophic velocity median, m/s |
|---|---|---|---|---|
| 18 OCI positions | 13 → 17 | +0.91 → +1.04 | +0.0199 → +0.1415 | +0.6513 → +0.6611 |
| 13 OCI positions | 17 → 21 | +0.91 → +1.62 | +0.0786 → +0.0373 | +0.6785 → +0.7194 |

The historical reference compares the difference in each pair's spatial medians with overlapping four-day changes in May 2017–2024 at the same positions. The eastward component increased by 0.1216 m/s on May 13–17. The largest positive change in 216 prior May intervals was 0.0896 m/s. Three historical changes were larger in absolute magnitude. These ranks describe an event and positions selected using the optical record. Overlapping time intervals and spatial correlation preclude an independent-sample significance calculation. The SST change ranked at the 62nd percentile on May 13–17 and the 94th percentile on May 17–21. Northward flow persisted during both intervals. Established observations describe weekly meanders and frontal eddies in this region. [Lee et al.](https://doi.org/10.1029/91JC02450)

The GOES-19 ABI column aerosol optical depth (AOD) at 550 nm rose and fell over the OCI positions. This table uses pixels with the highest quality flag (DQF 0) on both dates and reports medians over matched positions. The May 21 OCI daily product can contain more than one passage. Both GOES scan times are retained.

| OCI pair | GOES scan starts, UTC | Matched pixels | Median AOD, first → last | Median paired difference | Pixels changing in the indicated direction |
|---|---|---|---|---|---|
| May 13–17 | 17:21 → 18:06 | 10 | 0.0795 → 0.2397 | +0.1668 | 10 increased |
| May 17–21 | 18:06 → 17:11 | 7 | 0.2286 → 0.1299 | −0.1093 | 7 decreased |
| May 17–21 | 18:06 → 18:46 | 7 | 0.2286 → 0.1704 | −0.0459 | 6 decreased |

A fixed-clock sensitivity check at approximately 18:01 UTC gives an increase at 11 of 11 high-quality sites for the first pair, and a decrease at 7 of 8 for the second. The physical current change and the change in atmospheric loading coincide with the optical disturbance. Aerosol optical depth leaves aerosol composition and the size of any atmospheric-correction error unresolved. Daily NOAA Hazard Mapping System smoke polygons intersected zero paired sample positions on May 13, 17, or 21. These analyst-drawn polygons may miss diffuse or cloud-obscured smoke. [ABI product](https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.ncdc%3AC01511) · [NOAA quality flags](https://www.star.nesdis.noaa.gov/atmospheric-composition-training/glossary.php) · [HMS](https://www.ospo.noaa.gov/products/land/hms.html)

## Southern California

On March 11–16, 2024 the nearby 0.25° OISST cell averaged +0.98 °C relative to its seasonal climatology. On April 10–15 it averaged −0.18 °C. The mean eastward and northward geostrophic components were −0.007 and −0.071 m/s in March and +0.027 and +0.048 m/s in April. The mean magnitude of the daily 10 m wind vectors was 5.22 and 4.02 m/s, respectively. The OISST, current and wind grid cells span different spatial scales. Resolving circulation at the optical platform, tidal currents, wave resuspension or a nearshore discharge plume requires finer observations.

A separate March 8 oil sheen reported near Huntington Beach lay about 12 km from the optical station. The Coast Guard closed the response March 11 and distinguished the sheen from a reported produced-water discharge at Platform Elly. The archive evidence does not connect either incident to a station spectrum. [March 8 response](https://www.news.uscg.mil/Press-Releases/Article/3701744/unified-command-is-responding-to-oil-sheen-offshore-of-huntington-beach/) · [Final response](https://www.news.uscg.mil/Press-Releases/Article/3703679/update-final-unified-command-concludes-response-to-sheen-offshore-of-huntington/)

## Datasets and replay

NOAA/NCEI OISST v2.1 from NOAA PSL, Boulder, Colorado, provides daily temperature anomalies relative to the 1971–2000 reference. NOAA/AOML CoastWatch combines altimetry and climatology to estimate daily 0.2° surface geostrophic currents. NOAA/NCEP/NCAR Reanalysis 1, distributed by NOAA PSL, provides daily zonal and meridional 10 m wind components on its coarse Gaussian grid (approximately 1.9°). Components are positive toward the east and north. Each optical coordinate is assigned to its nearest native cell. Spatial summaries give each OCI position equal weight and record the number of distinct native cells. [OISST](https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.highres.html) · [Geostrophic currents](https://coastwatch.pfeg.noaa.gov/erddap/griddap/miamicurrents.html) · [Reanalysis](https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.html)

The scripts read the optical pair file supplied as an argument and write CSV evidence and SHA-256 source receipts in this directory. The ABI script requires Python GDAL bindings with netCDF and PROJ. The smoke script uses the GDAL command-line program ogr2ogr. The footprint script requires NumPy.

```sh
python physical_context.py
python footprint_environment.py /path/to/footprints.csv
python aerosol_context.py /path/to/footprints.csv
python smoke_context.py /path/to/footprints.csv
```

The footprint table includes the number of supporting native cells and the interquartile range on each day. The prior-May CSV supplies the daily footprint medians used to calculate historical ranks. Exact source URLs and byte sizes accompany each SHA-256 digest. The digest of an unmodified raw granule can be checked with sha256sum. Time-coordinate assertions prevent date-index shifts during replay. The GOES projection transforms longitude and latitude into the ABI geostationary fixed grid before sampling and applies the file's scale and offset.

## Adriatic and Black Sea stations

The daily Mediterranean and Black Sea analyses from Copernicus Marine give current, salinity and potential temperature in the shallowest modeled layer. Each nearest station cell is compared with its surrounding five-by-five grid of water cells. The model grids have nominal resolutions of 4.2 km and 2.5 km, respectively. The numbers below are daily means. [Mediterranean analysis](https://data.marine.copernicus.eu/product/MEDSEA_ANALYSISFORECAST_PHY_006_013/description) · [Black Sea analysis](https://data.marine.copernicus.eu/product/BLKSEA_ANALYSISFORECAST_PHY_007_001/description)

| Station and interval | Modeled salinity | Modeled temperature, °C | Modeled eastward velocity, m/s | Modeled northward velocity, m/s | OISST anomaly, °C |
|---|---|---|---|---|---|
| Acqua Alta, November 6–7, 2025 | 36.736 → 36.573 | 16.959 → 16.880 | −0.030 → −0.151 | −0.026 → −0.061 | −0.37 → −0.25 |
| Section-7, March 9–12, 2026 | 15.608 → 16.569 | 5.761 → 7.081 | −0.170 → −0.081 | −0.080 → −0.018 | −0.33 → +0.66 |

At Acqua Alta (45.3139°N, 12.5083°E), neighborhood median westward speed increased from 0.022 to 0.104 m/s. The short-term freshening in the station cell was modest. The local 2025 hydrographic survey describes seasonal cooling, wind mixing and reduced autumn stratification. Sediment transport during northeasterly Bora and southeasterly Sirocco winds is established in the northern Adriatic. At a water depth of about 16 m, resuspension and coastal inputs remain plausible. [CNR 2025 survey, Fig. 8](https://www.ismar.cnr.it/web-content/wp-content/uploads/2026/02/TR-46.pdf) · [CNR station](https://www.ismar.cnr.it/en/what-we-do/oceanographic-infrastructures/acqua-alta-oceanographic-tower/) · [Adriatic sediment dispersal](https://doi.org/10.1029/2006JC003868)

At Section-7 (44.5458°N, 29.4466°E), the neighborhood median salinity rose from 15.61 to 16.46. The interquartile range shifted from 15.44–15.93 to 16.18–16.61. OISST independently records warming. These data support a transition toward warmer, saltier surface water at the time of the optical event. The salinity and current evidence comes from the regional model. A field salinity transect together with optical absorption and backscattering would help separate changing plume influence from mixing, suspended matter and changes in the phytoplankton. River plumes and wave-driven sediment transport on the western Black Sea shelf have prior observational support. [Black Sea observations](https://doi.org/10.3389/fmars.2022.998970)

Run `python regional_hydrography.py` with NumPy and the system Blosc library installed. The script verifies the coordinate arrays, reads the nearest-to-surface model level, records hashes of the source metadata and data chunks, and includes daily OISST values. The CSV covers November 4–9, 2025 and March 7–15, 2026 to show event context.

The Copernicus Black Sea multi-year reanalysis supplies 924 overlapping three-day intervals for March 1993–2025 at the nearest surface cell. Eight intervals warmed by at least 1.320 °C, 21 gained at least 0.961 in salinity and two matched both unrounded March 2026 thresholds. They were March 5–8 and March 13–16, 2015. The joint hydrographic change is unusual and has historical modeled analogues. For comparability, inputs retain the distinction between the historical multi-year reanalysis and the 2026 operational analysis. Their model configurations and resolutions in depth may differ. A historical hydrographic analogue has no accompanying optical classification in this test. Run `python black_sea_reference.py` to replay the same-cell, duration- and season-matched reference. [Black Sea multi-year reanalysis](https://data.marine.copernicus.eu/product/BLKSEA_MULTIYEAR_PHY_007_004/description)
