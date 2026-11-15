# Does the July 13, 2021 far-side eruption belong to NOAA active region 12849?

**Provisional source identification—not confirmed.** Solar Orbiter recorded an estimated M3 flare during the July 13, 2021, fast halo eruption. Its measured X-ray position angle, 104.99° eastward from solar north, favors the region later numbered NOAA 12849. Propagating the July 22 NOAA region centers into Solar Orbiter’s July 13 view gives a 2.92° position-angle residual for AR 12849 and 14.00° for its closest numbered alternative, AR 12847. The proposed identification would place a third fast far-side halo eruption in the sequence that includes the July 15 and July 17 eruptions previously attributed to AR 12849. The published X-ray centroid lies above Solar Orbiter’s apparent limb, leaving an important ambiguity in longitude.

## Archived observations

Times below refer to July 13 UTC. The X-ray entries use the science FITS header’s 75.31 s Sun–spacecraft versus Sun–Earth light-time difference. STEREO and radio times retain their catalog conventions. Minute-level ordering between spacecraft is approximate.

| Observation | UTC | Measurement |
| --- | --- | --- |
| STEREO-A/EUVI | ≈18:25 | DONKI analysts report field opening beyond the east limb. The STEREO COR1 preliminary catalog lists a northeast–east ejection at 18:26. |
| STIX, 15–25 keV | 18:28:29–18:28:53 | Peaks from the four-second quick-look and approximately one-minute science-pixel products. The latter reaches 286.35 counts/s against a 17:40–18:00 spacecraft-time median of 18.31 counts/s. |
| SOHO/LASCO C2+C3 | 18:30:20 | Linear back-extrapolation of nine measured leading-edge heights to one projected solar radius. Its 1,100.2 km/s sky-plane slope reproduces the catalog speed of 1,101 km/s. This extrapolated time is model dependent. |
| STIX, 4–10 keV | 18:35:45–18:38:53 | Quick-look and science-product peaks. The latter reaches 62,906 counts/s, 789 times the baseline. Both telemetry products observe the same photons. |
| Wind/STEREO Waves catalog | 18:39–18:50 | Decametric–hectometric type-II feature drifts from 2.1 to 1.1 MHz. The catalog calls the source behind the east limb and gives no NOAA ID. |
| STEREO/SOHO coronagraphs | 18:53 / 19:00:05 | DONKI and LASCO first appearances, respectively. DONKI’s modeled direction is W178, N1. It measures the ejecta’s propagation, which can differ from its footpoint. |

A second CME changes the interpretation of the coronagraph images. DONKI lists a 16:12 ejection toward W121, S23 and describes overlap in STEREO/COR2. The later 18:53 ejection is the one contemporaneous with the hard X-ray, radio and limb signatures. The evidence does not resolve a unique radio-emitting segment of the shock.

## Geometry and controls

STIX lists helioprojective offsets of (−1103.8″, −295.5″). The measured centroid radius is 1.0312 times the Sun’s apparent radius. The imaging catalog uses a spherical screen for off-disk positions. Its July 13 Carrington coordinates, 52.45°, −14.99°, depend on that assumption. The comparison instead uses the measured position angle, which is insensitive to a radial change in source height. We projected all six sunspot regions and the additional H-alpha plage in the July 22 USAF/NOAA Solar Region Summary back into Solar Orbiter’s view, holding their tabulated Carrington coordinates fixed.

| July 22 region | Carrington longitude, latitude | July 13 predicted position angle | Position-angle residual | Radial height matching observed radius |
| --- | --- | ---: | ---: | ---: |
| AR 12849 | 58°, −18° | 107.90° | 2.92° | 0.036 R☉ |
| AR 12847 | 73°, −28° | 118.99° | 14.00° | 0.089 R☉ |
| Next five regions / plage | See output | | 44.74° or greater | |

For each region we perturbed longitude by ±5° and latitude by ±3° in one-degree steps. AR 12849 spans 0.00–6.01° residual, versus 10.37–17.96° for AR 12847. These offsets illustrate sensitivity to rotation and to displacement of the flare from the tabulated spot. They are not measured uncertainties. For the July 13 centroid alone, radial heights up to 0.06 R☉ admit Carrington longitudes 38.48–66.15°. Allowing 0.10 R☉ expands this to 31.33–73.37°. Off-radial loops, centroid errors, and an unnumbered region at similar latitude all remain viable. The July 22 list samples regions still visible eight days after the eruption. It cannot represent a region that dispersed before reaching the Earth-facing disk.

The search included 29 catalog backside type-II events in 2021–2024. Ten had at least one localized M-equivalent or brighter STIX peak in the exploratory −35 to +25 min window. Some were known, and some matched the wrong hemisphere. All triggers, duplicate peaks and misses are in `measurements.json`. In July 2021 the STIX science list has 59 entries, including three localized M-equivalent or brighter peaks, on July 13, 15 and 17. A ±30 min window about the July 13 peak occupies 60/1,440 one-minute clock phases on the same day. The three July peaks occupy 180/44,640 across the month. STIX science coverage and source brightness vary. The selection was exploratory, so these occupancies cannot supply a discovery significance.

## Prior attribution and limits

The July 13 halo, radio burst and STIX flare were each cataloged before this analysis. The July 2021 pairing of July 15 and 17 with AR 12849 is explicit in the BepiColombo/EPD intercalibration paper. The July 17 gamma-ray detection and its STIX location were already reported. The new proposition here is that July 13 belongs to that same region. The NASA radio and halo lists, DONKI’s July 13 entry, the STEREO COR1 list and the July 15/17 analyses examined here did not supply that numbered July 13 identification. The [CoSEE-Cat DR1 main table](https://coseecat.aip.de/files/dr1/csv/CoSEE_main.csv) contains 303 solar energetic-electron events and has no July 13 entry. Its next July entries are July 16, 17 and 18. This bounded check cannot establish publication priority.

An independent Earth-side record places the region at S18E75 with Carrington longitude 58° in the July 22 NOAA report. DONKI also identifies subsequent July 23/24 activity with AR 12849. The raw STIX timing and the LASCO fit were reproduced here. The X-ray image centroid comes from the published STIX imaging pipeline. The proposed AR 12849 attribution remains **provisional**, not a confirmed association: projected position angle ranks candidate regions but cannot uniquely identify the eruption footpoint. An independent image reconstruction and a photospheric identification for July 13 remain outstanding.

## Reproduction

Run `python reproduce.py` from any directory. It regenerates `measurements.json` from the frozen extracts using the Python standard library. `python extract_fits.py` downloads and checksum-verifies the two original STIX FITS products, then regenerates the compact channel sums and headers. This step uses NumPy and Astropy. Rates are total downlinked counts per accumulation divided by its duration. The quoted STIX GOES-equivalent classes come from the catalog and carry calibration uncertainty. The archive records full payload URLs, sizes and SHA-256 values in `inputs/sources.json`. The shipped STIX list includes the complete month of July 2021 and all catalog peaks within two hours of the 29 radio onsets.

[STIX imaging catalog and coordinate convention](https://github.com/hayesla/stix_flarelist_science/tree/5c988cf41764f8856e3fc018a72e4f48ca43736a) · [NASA type-II catalog](https://cdaw.gsfc.nasa.gov/CME_list/radio/waves_type2.html) · [July 13 DONKI record](https://ccmc.gsfc.nasa.gov/DONKI/view/CME/17292/-1) · [STEREO-A COR1 July catalog](https://stereo-ssc.nascom.nasa.gov/cor1/cor1/catalog/cme/2021/HongXie_COR1_preliminary_event_list_2021-07.html) · [NOAA July 22 region summary](https://www.ngdc.noaa.gov/stp/space-weather/swpc-products/daily_reports/solar_region_summaries/2021/07/20210722SRS.txt) · [Prior July 15 and 17 identification](https://doi.org/10.1051/0004-6361/202660642) · [Published analysis of July 17](https://arpi.unipi.it/handle/11568/1161207)
