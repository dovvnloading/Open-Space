# 4C +49.29 flared again in gamma rays and visible light in February 2026

The radio quasar 4C +49.29 (TXS 1729+491, redshift 1.038) was bright at two independent telescopes during 2026 February 2–9 UTC. Fermi-LAT placed the gamma-ray excess 0.061° from the VLBI position. ZTF recorded positive difference flux in both g and r on two nights at the radio source. The optical median flux during the second flare was 4.52 times the pre-February 2025 baseline in g and 5.67 times in r. **Fermi had already announced the first event** in February 2025, and both gamma-ray flares were already listed in the public FAVA archive. **The new result here is the optical association of the 2026 recurrence** with the two-band ZTF brightening—not the discovery of either LAT flare or evidence for a period. [Original telegram](https://www.astronomerstelegram.org/?read=17048) · [FAVA 2026 flare](https://fermi.gsfc.nasa.gov/ssc/data/access/lat/FAVA/SourceReport.php?week=914&flare=11) · [ALeRCE optical alert record](https://alerce.online/object/ZTF19abddedc)

### Measurements

| FAVA week, UTC | Offset from VLBI position | LAT 95% radius | Likelihood TS, 0.1–0.8 / 0.8–300 GeV | Fixed-position high-energy counts, observed / expected | Same-week median ZTF g / r |
|---|---:|---:|---:|---:|---:|
| 2025 Feb 17–24 | 0.048° | 0.101° | 81.15 / 71.16 | 14 / 1.59 | 18.334 / 18.292 mag |
| 2026 Feb 2–9 | 0.061° | 0.119° | 67.11 / 80.80 | 17 / 1.22 | 17.918 / 17.341 mag |

The 2026 flare and the 2025 flare rank first and second by high-energy aperture significance in the 919 available weekly extractions at the radio position. FAVA localized each independently within the relevant week. The Sun lay over 82° away in both high-energy fits. The nearest gamma-ray catalog positions to the radio core are FL16Y J1728.3+5013 at 1.181° and 4FGL J1728.3+5013 at 1.187°. FL16Y's integration ended in August 2024, before these detections. [2025 FAVA record](https://fermi.gsfc.nasa.gov/ssc/data/access/lat/FAVA/SourceReport.php?week=864&flare=13) · [2026 FAVA record](https://fermi.gsfc.nasa.gov/ssc/data/access/lat/FAVA/SourceReport.php?week=914&flare=11) · [FL16Y release](https://fermi.gsfc.nasa.gov/ssc/data/access/lat/fl16y/)

ZTF's direct photometry supplies 254 g and 260 r pre-trigger weekly medians, centered at 19.556 and 19.226 mag. No pre-trigger weekly median was as bright as the 2026 same-week median in either filter. The baseline gives equal weight to weeks with different visit counts. The two 2026 nights have high-confidence positive ZTF alerts with corrected magnitudes, DRB ≥0.9, no dubious flag, <0.15 mag uncertainty and coordinates within 1 arcsec of the VLBI core. On MJD 61076 the object reached g = 17.796 ± 0.064 and r = 17.176 ± 0.045. The accepted alerts and direct photometry agree to within 0.023 mag on four matched exposures during the 2025 flare. The two optical reductions share ZTF exposures. The independent instrument is LAT. [IRSA ZTF light curves](https://irsa.ipac.caltech.edu/docs/program_interface/ztf_lightcurve_api.html) · [Alert detections](https://api.alerce.online/ztf/v1/objects/ZTF19abddedc/detections) · [Radio positions](https://astrogeo.org/sol/rfc/rfc_2026c/)

### A longer high-energy excess surrounds the two localized peaks

At the fixed radio position the FAVA high-energy aperture also shows an excess outside the two selected weeks. Eight comparison positions, at 3° angular separation and position angles spaced by 45°, use the same weekly extraction and phase boundaries. The displayed ratio is the sum of observed counts divided by the sum of FAVA's expected counts.

| Interval | Covered weeks | Target counts | FAVA expectation | Target ratio | Range of eight comparison ratios |
|---|---:|---:|---:|---:|---:|
| Before week 864 | 840 | 1,006 | 986.29 | 1.02 | 0.95–1.10 |
| Between weeks 864 and 914, both excluded | 42 | 150 | 61.23 | 2.45 | 0.86–1.29 |
| After week 914 through week 951 | 35 | 91 | 50.32 | 1.81 | 0.70–1.24 |

The elevated counts in both intervals support additional gamma-ray activity in this direction. FAVA aperture counts are spatially smoothed. A time-resolved likelihood analysis that models neighboring sources would be needed to apportion the between-flare counts to the quasar. The FAVA sky list has 925 numbered weeks between August 2008 and October 26, 2026. The point extraction returns 919 weeks. Missing weeks are excluded throughout. The two observed epochs establish recurrence. They do not measure a period or continuous emission. [FAVA method](https://fermi.gsfc.nasa.gov/ssc/data/access/lat/FAVA/About.html)

### Selection and prior art

The exploratory survey began with 9,471 FAVA flare rows. A blank `fglassoc` designation (the literal value `None`) and absolute catalog Galactic latitude above 10° retained 1,110 rows. Requiring the fitted flare center to be farther than 0.25° from all DR4 and FL16Y positions left 288. Matching the nearest of 23,866 RFC entries within 0.2° gave 22 flares at 21 distinct radio sources. Only J1730+4906 repeated. The same two events are the only two flare-list centroids within 1° of the quasar, before the selection cuts. The optical light curves were inspected after selecting the radio/gamma-ray repeat. These exploratory cuts, the optical historical ranks and the eight neighboring apertures do not provide a globally calibrated false-association probability. The Fermi discovery telegram independently identifies the quasar at the 2025 epoch. [All 22 accepted survey rows](survey.json)

On 2026 November 15, searches using the aliases `4C +49.29`, `4C 49.29`, `TXS 1729+491`, `J1730+4906` and `ZTF19abddedc` with the February 2025 and February 2026 dates, gamma-ray and optical terms recovered the February 2025 Fermi telegram. It explicitly dates the first LAT flare and invites multiwavelength observations. FAVA already publishes both individual flare records, including the February 2026 one. We found no source-specific publication or telegram identifying the February 2026 recurrence with the corresponding ZTF event in the searches performed. Independent publication priority remains open. The measured addition here is the cross-identification of that second outburst in two wavebands, plus the comparison-field result between and after the peaks. [ATel 17048](https://www.astronomerstelegram.org/?read=17048)

### Reproduction

Run `python replay.py` in this directory to regenerate `measurements.json` offline using the standard library. `inputs` contains the two FAVA reports, relevant RFC row, all 919 weekly point extractions for each of the nine positions, the ZTF cone light curve and the ALeRCE detection series. The compressed photometric and aperture tables preserve the columns used in the analysis. `sources.json` records the original URLs and SHA-256 hashes of the full upstream responses. The URLs for live tables can change contents as observations arrive.

Run `python survey.py` with NumPy and SciPy for the offline all-sky replay. The three compressed `sky-*.csv.gz` tables contain every flare row and every catalog position used in the selection. The output reproduces `survey.json`, including the 22 accepted pairs. To fetch the upstream catalogs and rerun, use `python survey.py --live-cache /path/to/cache` with Astropy installed. The four raw files are cached as `flare_list`, `4fgl`, `fl16y` and `rfc`. Compare their hashes with the manifest when using the live catalogs.
