# PACE matchup screen and historical field comparison

NASA's official remote-sensing reflectance matchup export, generated June 24, 2026, contains 932 OCI/field pairs acquired March 6, 2024 through April 28, 2026. Each OCI record includes all 61 samples from 400 to 700 nm at 5 nm spacing. The seasons DJF, MAM, JJA and SON contribute 90, 268, 343 and 231 matchups respectively. The record is opportunistically sampled, with concentration at coastal platforms. The export includes parent observation identifiers and derived reflectance. Its header omits the ocean-product processing version. The separate validation viewer identifies its own OCI release as V3.2. [Official data and method](https://pace.oceansciences.org/pace_rrs.htm) · [SeaBASS validation](https://seabass.gsfc.nasa.gov/wiki/validation_description) · [Viewer](https://seabass.gsfc.nasa.gov/validation-viewer/)

Full field spectra on the same grid exist for 118 rows across 36 cruise labels. For these rows, the median OCI-to-field spectral angle is 3.27°, with quartiles of 1.84° and 7.55°. The angle uses the entire 61-component vector, including negative retrievals, normalized by its Euclidean length. All 118 records and their coordinates, parent identifiers and angle are in `complete_field_spectral_pairs.csv`. USC_SEAPRISM itself provides multispectral field measurements. Blanks in its spectrum table identify unmeasured bands, and reflectance and uncertainty are in sr⁻¹.

## Fixed-platform screen

Six AERONET-OC platforms contributed 743 matchups or 727 distinct site-days. For repeat overpasses on one UTC day, the screen selects the matchup with the smallest absolute field/satellite time separation. It evaluates consecutive platform visits 1–8 days apart. A cruise name alone cannot identify a repeat observation at a fixed point.

For wavelengths a < b < c, define C as ln(Rrs at b) minus the linear interpolation in wavelength of ln(Rrs) at a and c. A uniform multiplicative brightness change cancels. The blue and red triplets are (443, 490, 560) and (560, 620, 665) nm respectively. A nonpositive input leaves that metric undefined. The screen considered 644 intervals, with 599 usable paired blue deltas. Both sensors returned the same sign in 360 intervals. Ranking those by the smaller absolute change across the two sensors selects March 11–15 at USC_SEAPRISM. The field data participated in selection. Associated scaled median absolute deviation scores and ranks are descriptive. Serial correlation, opportunistic sampling and searching multiple intervals preclude treating them as discovery probabilities.

| California interval, 2024 | OCI blue ΔC | Field blue ΔC | OCI broad ΔC | OCI and field endpoint scores | Other-year comparison days |
| --- | ---: | ---: | ---: | ---: | ---: |
| March 11–15 | −0.3771 | −0.3821 | −0.3424 | −5.33, −7.82 | 21 |
| April 10–15 | −0.3513 | −0.2323 | −0.3178 | −5.22, −10.41 | 21 |

The broad check takes the median of eight OCI log residuals at 475–510 nm, in steps of 5 nm, using a 445/560 nm log-linear continuum. It requires all eight residuals. Endpoint scores compare with the same site's other years within 30 calendar days either side of the endpoint, using (endpoint − median)/(1.4826 × MAD). The scripts align month/day across leap years. The brief March change appears in the broader spectral shape as well as the common sensor bands.

## Historical field data

The supplementary archive uses AERONET-OC Version 3 Level 2.0 daily mean normalized, bidirectionally corrected water-leaving radiances at USC_SEAPRISM. The matchup export labels its field rows V4 Level 1.5. The archival comparison is an alternative processing and time aggregation of the same station, so it does not add a separate physical observer. Fixed per-channel multiplicative normalizations cancel in temporal log contrasts. Data received for 2024 contain the 560 nm channel. Historical radiometers used 551 nm instead. Exact-band historical comparisons therefore use 443/490/667 nm in all years. They measure a related broad-band shape.

| Interval | Daily ΔC with 560 nm | Daily ΔC with 667 nm | Earlier exact-lag changes as negative | Lowest earlier ΔC |
| --- | ---: | ---: | ---: | ---: |
| March | −0.3685 | −0.2920 | 0 of 77 | −0.2683 |
| April | −0.2300 | −0.0467 | 16 of 51 | −0.3739 |

Each reference endpoint falls within ±30 calendar days of the 2024 endpoint, before 2024. The 163 March reference endpoint days span eight years with coverage between 2012 and 2023. The four-day comparison permits overlapping intervals. Daily C443/490/560 shifted from 0.3534 on March 13 to 0.0739 on March 14 and −0.0140 on March 15, then rose to 0.3007 on March 16. Absolute values calculated from normalized radiance and reflectance have different band-dependent offsets. [AERONET-OC data service](https://aeronet.gsfc.nasa.gov/print_web_data_help_v3_seaprism.html)

## Quality and interpretation

- On March 15 only 14 of 25 OCI pixels contributed. The other three California overpasses retained all 25. Absolute QWIP, a published spectral shape quality index, exceeded 0.2 on March 15 (+0.201) and April 10 (−0.203). Both endpoint pairs fail that sensitivity screen. The primary source supplies summaries without per-pixel flags or masks. [QWIP paper](https://doi.org/10.3389/frsen.2022.869611)
- March 15 OCI Rrs at 443, 490 and 560 nm was 0.00374, 0.00344 and 0.00307 sr⁻¹. The field instrument measured 0.00180, 0.00190 and 0.00239 sr⁻¹. OCI Rrs665 was negative at each starting overpass. The local temporal blue-curvature agreement therefore coexists with substantial discrepancies in the absolute spectra and uncertainty about correction errors.
- With at least 20 valid pixels at both ends and absolute QWIP at most 0.2, 228 of 644 intervals remain, including 220 with usable blue deltas in both instruments. The strongest same-sign blue change then occurs at Irbe Lighthouse during June 20–24, 2024 (smaller absolute ΔC of 0.264). Its endpoint scores are −0.78 and −1.91. These are sensitivity checks against sparse spatial support and unusual spectral shape.
- The field record supports a short-lived local water-color change. Its spectral shape alone cannot distinguish runoff, resuspended material, changing colored dissolved matter, advection or phytoplankton optical changes. The individual cruises sample a range of conditions and do not provide dates repeated over a geographically uniform study area.

## Reproduction

Run from this directory with Python 3.12, pandas 2.2, NumPy and curl. The scripts create the output directory if needed. Input fields with fill value −999 become missing. SHA-256 values identify the NASA release and individual annual AERONET responses. The NASA file is 5,896,423 bytes. A new release requires an explicit choice to rerun the exploratory selection using updated input.

```sh
mkdir -p raw
curl --fail --location --retry 2 https://pace.oceansciences.org/data/matchup_Rrs.csv -o raw/matchup_Rrs.csv
python screen_matchups.py raw/matchup_Rrs.csv --out reproduced
python field_history.py raw --out reproduced
```

`california_matchups.csv` and `california_spectra.csv` give the four actual overpasses, band values, uncertainties and geometry. `fixed_platform_intervals.csv` includes all eligible intervals and the quality-screen inputs. `site_season_coverage.csv` uses UTC platform-days. `hyperspectral_coverage.csv` preserves original campaign labels and meteorological seasons. The field archive outputs include annual source URLs and hashes, every historical comparison interval, event-day series and summary statistics.

NASA's PACE mission, Ocean Biology Processing Group, SeaBASS, AERONET-OC and contributing field teams provided the measurements. Site PIs Burton Jones and Matthew Ragan are credited in the USC_SEAPRISM file. The AERONET-OC data service requests offering the site PIs co-authorship for research use. [Data use and PI guidance](https://aeronet.gsfc.nasa.gov/print_web_data_help_v3_seaprism.html)
