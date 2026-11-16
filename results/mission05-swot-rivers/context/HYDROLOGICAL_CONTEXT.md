# Upper Zambezi hydrometric and precipitation evidence

## Seasonal discharge and recovery

The Chavuma gauge recorded a weak and late 2024 flood pulse. The full daily record for February–May shows a March mean of 367 m³/s, 82% below the published March climatology. The deficit narrowed in May. The controls below use the station column of Lourenco et al. table 5. Its caption gives 1959–2002 for Chavuma, while the supplementary figure labels the start as 1960. The historical monthly means have no published year-by-year spread in this table. [Historical station and model columns](https://doi.org/10.5194/hess-29-4557-2025)

| Month | Chavuma 2024 mean (m³/s) | Historical mean (m³/s) | Departure |
|---|---:|---:|---:|
| February | 415 | 1,411 | −71% |
| March | 367 | 2,031 | −82% |
| April | 704 | 1,770 | −60% |
| May | 588 | 684 | −14% |

Matching only April dates reported in both years gives the following station means. The sets of days differ between comparisons. A sustained increase in river flow followed the drought at both sites, with a particularly large upstream response in April 2025. A 2025 annual discharge total or a permanent change in channel behavior cannot be inferred from a record ending on May 2.

| Station | Comparison | Paired days | Earlier mean (m³/s) | Later mean (m³/s) | Change |
|---|---|---:|---:|---:|---:|
| Chavuma | April 2023–2024 | 26 | 2,237 | 690 | −69% |
| Chavuma | April 2024–2025 | 23 | 775 | 3,034 | +291% |
| Nana’s Farm | April 2023–2024 | 29 | 2,010 | 693 | −66% |
| Nana’s Farm | April 2024–2025 | 29 | 692 | 2,173 | +214% |

At Chavuma, the highest reported March–May daily flow in 2023 was 2,614 m³/s on April 15. In 2024 it was 973 m³/s on April 25. ZRA's annual record at Nana’s Farm places the downstream maxima on May 27, 2023 (2,718 m³/s) and June 7, 2024 (938 m³/s). These peak separations were 42 and 43 days. The 2023 upstream series begins in late March and has missing April days. Different tributary inputs and temporary storage in the Barotse floodplain can shift the downstream hydrograph. The observed separations are compatible with earlier dry-year timing estimates, including about 1.5 months in 2014/15. They provide no persuasive evidence of changed flood-wave travel time. ZRA gives 1,109 m³/s as the long-term annual mean at Nana’s Farm and 521 m³/s for 2024. The seasonal drought and its downstream consequences were already in the public record. [ZRA annual report, table 1](https://zambezira.org/sites/default/files/ZAMBEZI%20RIVER%20AUTHORITY%202024%20ANNUAL%20REPORT_FINAL.pdf) · [Zambia Ministry water-resources assessment, section 6.3](https://www.mwds.gov.zm/wp-content/uploads/2025/06/Annexure-C-CZWTS-Water-resources-availability-assessment.pdf)

## Rainfall and downstream water levels

CHIRPS v2 final observations place February 2024 last by rainfall in the 45 Februaries from 1981 to 2025 in all five administrative areas below. The nested Angolan areas and neighboring Zambian provinces share weather and source pixels. The ranks describe this event and should not be multiplied or interpreted as five independent significance tests. October–March rainfall includes the preceding calendar year's October–December. All numbers are millimetres. [WFP Angola rainfall archive](https://data.humdata.org/dataset/ago-rainfall-subnational) · [WFP Zambia rainfall archive](https://data.humdata.org/dataset/zmb-rainfall-subnational)

| Area (historical P-code) | Feb normal, 1991–2020 | Feb 2023 | Feb 2024 | Feb 2025 | Oct–Mar normal | 2023/24 Oct–Mar | 2024/25 Oct–Mar |
|---|---:|---:|---:|---:|---:|---:|---:|
| Alto Zambeze, Angola (AO15128) | 204 | 119 | 111 | 188 | 1,099 | 879 | 1,038 |
| Moxico, Angola (AO15) | 180 | 104 | 95 | 174 | 1,012 | 769 | 916 |
| North-Western, Zambia (ZM108) | 219 | 133 | 73 | 219 | 1,081 | 787 | 930 |
| Western, Zambia (ZM110) | 168 | 91 | 35 | 169 | 822 | 443 | 652 |
| Southern, Zambia (ZM109) | 157 | 121 | 23 | 215 | 692 | 426 | 611 |

Alto Zambeze covers part of the Angolan headwaters upstream of Chavuma. Its October–March total increased 18% from 2023/24 to 2024/25. The contemporaneous increase in headwater rainfall supports a runoff-driven recovery. Rainfall averaged across administrative units differs from rainfall integrated over a contributing basin. The runoff response also depends on antecedent soil moisture, floodplain storage and rainfall timing.

On 12 matched calendar dates between May 8 and September 27, Kariba Reservoir stood 1.67–2.62 m lower in 2024 than in 2023. Both upstream discharge gauges were lower on every one of those dates. The reservoir incorporates local inflow, evaporation and hydropower withdrawals. ZRA reports the spillway gates were closed in 2024. Kariba provides a downstream basin-scale indicator. A separate river-stage gauge at Katima Mulilo read 1.48 m on January 31, 2024 versus 1.69 m in 2023, and 0.49 m on January 16, 2025 versus 1.33 m in 2024. Stage uses a local gauge datum and these observations are upstream of the Chobe confluence. For floodplain context, NASA compared Terra MODIS scenes at the Chobe–Zambezi confluence on April 13, 2023 and April 20, 2024. The published images show less open water in the 2024 Chobe channel and drier surrounding floodplain. These two scenes provide regional context without an independently classified area. [ZRA annual report](https://zambezira.org/sites/default/files/ZAMBEZI%20RIVER%20AUTHORITY%202024%20ANNUAL%20REPORT_FINAL.pdf) · [NASA MODIS comparison](https://modis.gsfc.nasa.gov/gallery/individual.php?db_date=2024-04-24) · [Hydrological Services Namibia, January 2025 bulletin](https://www.namibia-forum.ch/media/kunena/attachments/9603/2025_01_16_FloodBulletin.pdf)

## Independence and source handling

The comparison uses independently published 2023–25 gauge readings and the unconstrained SWOT v3 estimate. V3 uses δHBV2 and GRADES-hydroDL model priors. Both learned from historical gauge data. The published original δHBV2 training ran through 2000 and the final GRADES-hydroDL training through 2020, before every gauge date evaluated here. This establishes temporal separation from the training periods. Historical measurements from particular basins may have informed the model priors and hence the discharge scale. No skill claim here uses a gauge-constrained discharge estimate or treats agreement across SWOT algorithms as independent corroboration. [V3-specific prior description, section 2.1](https://doi.org/10.1029/2026GL122871) · [δHBV2 training](https://mhpi.github.io/datasets/global/) · [GRADES-hydroDL training, section 4](https://doi.org/10.1029/2024WR039764)

A separate field team measured 923.4 m³/s with an acoustic Doppler current profiler downstream of Chavuma on May 17, 2023. The ZRA gauge that day reported 941 m³/s. The difference was 17.6 m³/s (1.87% of the gauge value). The expedition reports a point-in-time transect. Its time of day, exact survey coordinates and uncertainty are absent from the tabulated data. No selected SWOT pass coincides with that date. This tests the gauge at one flow level. [Expedition supplement, table S1, site 24](https://doi.org/10.5194/hess-29-4557-2025-supplement) · [ZRA bulletin, page 1](https://zambezira.org/sites/default/files/Press%20Statement%20-%20Hydrological%20Situation%20at%20Kariba%20Dam_03%20June%202024-fm.pdf)

The daily archive spans March 29, 2023 to May 2, 2025 and contains 707 accepted Chavuma days and 702 Nana’s Farm days from 194 dated ZRA bulletins. Each accepted observation points to the bulletin IDs that report it. Each bulletin has its own URL and SHA-256. The compact source-table archive preserves the actual HTML table cells, date column, both flow columns and headers. Bulletins may share rows. The extraction uses only the current-period flow column, accepts dates inside the bulletin's local reporting window, and excludes 17 station-dates with conflicting published values. Eight date labels outside their reporting periods are retained in the exclusion log. Dates in the bulletin metadata use UTC+02:00. The authority supplies daily values without specifying the reading hour, averaging window or station-specific measurement uncertainty. Treat abrupt single-day jumps cautiously. ZRA identifies Nana’s Farm as 33 km upstream of Victoria Falls and subsequently combined its operational station name with Victoria Falls. Only bulletins that identify Nana’s Farm by name were included. [ZRA station information](https://www.zambezira.org/hydrology/river-flows/victoria-falls-nanas-farm-station)

Run `python analyze_context.py` to reproduce the precipitation and discharge summaries from the small data extracts. `extract_gauges.py RAW_HTML_DIR OUTPUT_DIR` reproduces the daily record from the listed HTML bulletins. The full-page hashes capture retrievals on November 16, 2026. ZRA updates sidebar content when later bulletins arrive, which can change full-page hashes. The archived tables exclude the live sidebar. `context_sources.csv` identifies static source files, their hashes and table locations. The rainfall extract retains all January–March and October–December CHIRPS final records for the five specified units from January 1981 through March 2025. The `rfh` column gives total rainfall in millimetres during days 1–10, 11–20 or 21–end of month. The full administrative-unit identifiers use the HDX boundary definitions contemporaneous with the rainfall data (Angola 2018 and Zambia 2023). Numeric daily comparisons have serial dependence, and only one drought-to-recovery transition is observed. The historical monthly means and annual peak dates support effect-size comparisons without a reliable recurrence probability for a structural change in river hydraulics.
