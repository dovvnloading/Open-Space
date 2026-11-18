# Public PDS ranges clip Voyager plasma observations

**Noteworthy correction · Correction.** The public NASA Planetary Data System (PDS) pages for the September 2026 Voyager interstellar plasma release display both endpoint hours 12 hours earlier than their source XML labels and version-matched Search API results. A client applying the displayed bundle limits to the source tables would omit Voyager 2's latest record. The displayed collection limits would omit that record and Voyager 1's latest record. Three comparison pages show the same 12-hour clock pattern. [September 2026 release](https://pds.nasa.gov/datasearch/subscription-service/SS-20260930.shtml) · [Public bundle](https://pds.nasa.gov/ds-view/pds/viewBundle.jsp?identifier=urn:nasa:pds:voyager-pws-vlism-density) · [Public collection](https://pds.nasa.gov/ds-view/pds/viewCollection.jsp?identifier=urn:nasa:pds:voyager-pws-vlism-density:data) · [Frozen inputs and receipts](catalog-sources.json)

## Version-matched UTC bounds

The comparison matches logical identifier and version across native XML, the NASA Search API and the public ds-view HTML. All ten start/end fields on the five pages have the same calendar dates, minutes, seconds and fractional seconds as their native labels. The hours are consistent with a 12-hour clock with the AM/PM marker omitted, including midnight rendered as noon. The four available version-matched Search API records agree with the native labels. MAVEN provides an additional XML-to-HTML check. The exact v2.26 API request returned HTTP 404, and the unversioned API response returned v2.24. This sample establishes the display pattern without determining which other pages are affected or how the site produces it. [Native and display comparisons](catalog-summary.json) · [Search API documentation](https://nasa-pds.github.io/pds-api/guides/search/endpoints.html)

| Archive record | Native start → displayed start (UTC) | Native stop → displayed stop (UTC) |
| --- | --- | --- |
| Voyager density bundle v2.0 | 2012-10-22 21:44:19.795 → 09:44:19.795 | 2025-09-01 13:29:00 → 01:29:00 |
| Voyager density collection v2.0 | 2012-10-22 21:44:19.795 → 09:44:19.795 | 2025-04-14 17:27:35.917 → 05:27:35.917 |
| Voyager 1 cruise trajectory v1.0 | 1977-09-09 00:00:00 → 12:00:00 | 2005-01-24 00:00:00 → 12:00:00 |
| Psyche SPICE v5.0 | 2023-10-13 15:21:48 → 03:21:48 | 2026-07-05 23:59:59 → 11:59:59 |
| MAVEN calibrated SEP v2.26 | 2014-09-20 00:00:00 → 12:00:00 | 2025-12-05 23:59:59.999 → 11:59:59.999 |

The Voyager 1 table contains **755 rows at 748 distinct epochs**. The Voyager 2 table contains **5 rows at 5 distinct epochs**. Both products appear in the collection inventory and both source tables agree with their own versioned XML and API endpoints. The native collection end is the Voyager 1 endpoint. It already omits the listed September Voyager 2 epoch, irrespective of the rendering discrepancy. The following counts come from applying each inclusive interval locally to all 760 rows. They do not describe the behavior of a NASA search or download. For reading both products together, the native bundle covers the two tables and the native product bounds cover each individually. [Source bytes](catalog-replay.zip) · [Coverage arithmetic](catalog-summary.json)

| Interval applied | Voyager 1 rows excluded | Voyager 2 rows excluded |
| --- | --- | --- |
| Native bundle | 0 | 0 |
| Displayed bundle | 0 | 1, September 1 at 13:29 UTC |
| Native collection | 0 | 1, September 1 at 13:29 UTC |
| Displayed collection | 1, April 14 at 17:27:35.917 UTC | 1, September 1 at 13:29 UTC |

## One frequency field needs its paper context

The standard cold-plasma conversion is nₑ = (fₚₑ / 8980 Hz)², with density in cm⁻³. Applying it to all 755 Voyager 1 rows reproduces the tabulated density to within **0.059%**. For four of the five Voyager 2 rows, the recomputed value lies inside the tabulated bounds. The exception is June 26, 2020 at 21:16 UTC. Its listed **3110 Hz** gives **0.11994 cm⁻³**, above the reported **0.081–0.094 cm⁻³** range around **0.087 cm⁻³**. [Product arithmetic](catalog-summary.json)

Kurth and Gurnett analyzed the peak near 21:15 UTC using the relative responses of the 1.78 and 3.11 kHz analyzer channels. Their estimated signal frequency of **2.65 kHz** gives 0.0871 cm⁻³ and explains the published density. The table lists 3110 Hz under both the marked frequency and the field the label describes as electron plasma frequency. For this event, a downstream calculation using that field alone would overestimate the published density by about **38%**. This difference concerns the interpretation of the archived frequency field. The source paper supports retaining its density estimate. The public errata discusses coordinates and fill values. [Kurth and Gurnett, 2020, Section 2 and Figure 3](https://doi.org/10.3847/2041-8213/abae58) · [Voyager 2 product label](https://pds-ppi.igpp.ucla.edu/data/voyager-pws-vlism-density/data/vg2-vlism-density-2019-2025.lblx) · [PDS errata](https://pds-ppi.igpp.ucla.edu/data/voyager-pws-vlism-density/document/errata.md)

## A separate check with Voyager 1’s magnetometer

A fixed density-to-field ratio trained on **21** jointly covered 2023 days gave **0.29083 cm⁻³/nT**. In 2025, four of six plasma-wave rows had enough direct 48-second MAG data in the specified ±12-hour window. These fell on four UTC days in three 14-day blocks. The fixed-ratio model gave **0.02164** median absolute log error on those days. A stationary density baseline fitted to the same 2023 dates gave **0.00841**. Two field-model predictions and three stationary predictions fell inside the published density bounds. [Frozen protocol](PROTOCOL.md) · [Native MAG files and digests](mag-sources.json) · [All matched windows](field-inputs.json) · [Scores and dates](field-summary.json)

The shifted windows share only two, two, one and one days with the simultaneous sample at offsets of −14, −7, +7 and +14 days, respectively. The ±6-hour sensitivity has three matched 2025 days in two blocks. In the embedded-field comparison, simultaneous and +12-hour MAG estimates have similar median absolute differences on the shared archived rows. The Voyager 1 table has no usable embedded field in 2025. Available data favor the stationary baseline on the frozen holdout. Sparse and discontinuous coverage prevents a timing-specific inference or an inference about a new pressure front. The 2020 compression and the field history through mid-2025 have already been studied. [Independent MAG archive](https://cdaweb.gsfc.nasa.gov/sp_phys/data/voyager/documents/vgrmag_website/data.html) · [Known 2020 pressure front](https://doi.org/10.3847/1538-4357/abeb6a) · [Park et al. 2026](https://doi.org/10.3847/1538-4357/ae63ba)

## Reproduce

Run these commands from this folder with Python 3.11 or newer. The catalog comparison reads the byte-verified ZIP of source labels, CSVs and API and HTML responses. The field comparison reads the committed matched windows.

~~~text
python audit_catalog.py --output catalog-summary.json
python analyze_field.py --output field-summary.json
~~~

To rebuild the matched windows, download the three [original MAG files](mag-sources.json) to a separate data directory and install NumPy and cdflib. The extraction verifies the original SHA-256 digests and sizes before reading the CDFs. This run staged the original CDFs on D. Separate native-data replays reproduced both the timestamp exclusions and the held-out scoring.

Run `python prepare_mag.py --data-root MAG_DATA_DIRECTORY --output replay-field-inputs.json` followed by `python analyze_field.py --input replay-field-inputs.json --output replay-field.json`.
