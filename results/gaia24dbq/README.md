# A recurrent optical variable at the position of 1eRASS J111110.6−434214

Gaia24dbq (AT 2024zci) brightened in 2019, faded, and reached a brighter optical state in 2024. Its position lies 1.44 arcsec from an X-ray source that faded by a factor of 9.67 in six months, then partially recovered. WISE measured higher flux in two infrared bands from June 2019 through June 2020, following the January–May 2019 optical high state. The measurements connect a published eROSITA variable with a subsequently reported Gaia alert. The source’s physical class and distance remain open.

## What the telescopes measured

| Measurement | Dates or visits | Result |
| --- | --- | --- |
| Gaia G band, median of nightly medians | 2015–2017, 11 dates | 20.82 mag |
| Gaia G band | January–May 2019, 8 dates | 19.54 mag |
| Gaia G band | 2020, 13 dates | 20.235 mag |
| Gaia G band | 2023, 14 dates | 20.292 mag |
| Gaia G band | June–October 2024, 5 dates | 19.00 mag, reaching 18.84 on October 22 |
| WISE W1, forced flux in AB nanomaggies | Six visits at MJD 57000–57999 versus three at MJD 58600–59099 | 1.64× and 1.84× higher in the two scan directions |
| WISE W2, same windows | Six versus three visits | 1.84× and 2.23× higher in the two scan directions |
| eROSITA, 0.2–2.3 keV, in units of 10⁻¹³ erg cm⁻² s⁻¹ | Five visits at MJD 58838, 59022, 59206, 59390 and 59575 | 4.31 ± 0.54, 0.45 ± 0.19, 1.75 ± 0.36, 2.42 ± 0.41 and 2.15 ± 0.40 |

Gaia returned numeric measurements in 224 of 364 transits, covering 96 distinct calendar dates. Numeric detections select brighter epochs near Gaia’s faint limit. Null and untrusted values have been retained and excluded from all magnitude summaries. The windows are retrospective descriptions of sampled states. Their spacings do not establish a period. The first and last eROSITA visits fall between the two optical high states. There is no simultaneous X-ray measurement here for the 2024 rise. Grotova et al. measured the X-ray variability and classified its five-visit shape as “other.” [Gaia transit data](https://gsaweb.ast.cam.ac.uk/alerts/alert/Gaia24dbq/lightcurve.csv) · [Gaia measurement documentation](https://gsaweb.ast.cam.ac.uk/alerts/docs/) · [Grotova et al. and the eRO-ExTra catalog](https://doi.org/10.1051/0004-6361/202451253)

## Association and controls

The search compared all 304 Legacy Survey optical counterpart positions in eRO-ExTra with the centers of all 27,076 entries in the Gaia Science Alerts list retrieved on 2026-11-15 UTC. Seven alert centers lie within 8 arcsec. Six have a TNS counterpart in the eRO-ExTra release. Gaia24dbq is the sole match with an empty TNS field in that release. Repeating the search at eight nearby control positions per host (60, 120, 300 and 600 arcsec east and west) yielded one match across 2,432 positions. The control demonstrates the dependence on position. It does not give a per-source false-association probability. The program returns all seven pairs, including the known counterpart of J0744 whose assigned Legacy object is the neighboring galaxy.

Gaia24dbq is 0.056 arcsec from the catalog’s proposed optical counterpart and 1.44 arcsec from the independent X-ray centroid (quoted positional uncertainty 1.12 arcsec). Gaia and the Legacy Surveys share an astrometric reference. The Gaia source is DR3 5388725488667359104. The nearest other Gaia DR3 source is 7.22 arcsec from the optical variable and 8.57 arcsec from the X-ray centroid at the reference epoch. Sky density, unresolved blends, and the broader instrumental X-ray PSF leave a chance alignment possible. [eRO-ExTra data dictionary](https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+A/693/A62) · [Gaia alert discovery](https://www.wis-tns.org/object/2024zci) · [Gaia archive](https://gea.esac.esa.int/archive/)

The infrared series uses the Legacy Surveys DR11 Tractor photometry of time-resolved unWISE coadds, brick 1677m437, object 5969 (DR10 object 1629). Each visit in the comparison has at least five contributing exposures. Nine other sources within 60 arcsec in the same brick passed the same field selection (full-catalog W1 above 30 and W2 above 20 nanomaggies). The two WISE scan directions repeatedly give different absolute fluxes at this position. Comparing each 2019–20 visit with earlier visits in the same season yields the ratios in the table. Three early visits enter each seasonal median, compared with one winter visit or two summer visits in the later window. The largest increase among the nine controls was 1.03× in W1 and 1.49× in W2 across either direction. Infrared magnitudes and the timing of maxima are sensitive to blending and the choice of release. Neighbor contributions to the modeled infrared profiles can be appreciable, especially in the faint state. The packet retains the blend metric, formal inverse variance, profile chi-square and exposure count for each epoch. The last usable DR11 visit for this target was MJD 60312, before the 2024 Gaia high state. At MJD 59011 the infrared flux was still high even though eROSITA measured a substantially lower flux eleven days later. The three surveys do not establish a single synchronous flare. [Tractor quantities and photometric units](https://www.legacysurvey.org/dr11/catalogs/) · [Legacy Surveys data via NOIRLab](https://datalab.noirlab.edu/data/legacy-surveys)

## Physical interpretations and discovery history

The parent catalog gives a model-based photometric redshift of 0.4865 for a point-like source. Its DR10 Legacy measurements have zero inverse variance in r, i and z, leaving g and infrared measurements for this source in that table. Gaia DR3 gives a position-only solution, with no parallax or proper-motion estimate. A distance and an absolute luminosity cannot be established from these data. Variability in an accreting galactic nucleus and accretion in a Galactic binary are both plausible. The X-ray power-law index was 1.90 ± 0.27. The measured X-ray and optical variability alone cannot identify a black hole or a tidal disruption.

Gaia Alerts reported the 2024 brightening. The Transient Name Server has no spectrum or classification report for AT 2024zci. The Rochester list labels the source a cataclysmic variable. The Gaia alert list records its class as unknown. On 2026-11-15 UTC, searches for the exact Gaia, TNS and eROSITA identifiers, alone and with the other identifiers, found the two discovery records and the eRO-ExTra catalog. A 3 arcsec search of VizieR VII/294 (Milliquas v8) and a Gaia DR3 `qso_candidates` lookup by source ID gave no entries. Literature incompleteness and alternate source names limit any priority claim. A spectrum in a bright state could establish a redshift from multiple consistent features or reveal a Galactic accretion disk. High-resolution X-ray imaging during a new bright state could locate the emission relative to the optical source. [TNS discovery record](https://www.wis-tns.org/object/2024zci) · [Gaia alert page](https://gsaweb.ast.cam.ac.uk/alerts/alert/Gaia24dbq/) · [Rochester listing](https://rochesterastronomy.org/sn2024/boneyard.html) · [Milliquas](https://cdsarc.cds.unistra.fr/viz-bin/cat/VII/294)

## Replay

Install NumPy and SciPy. Run from the repository root.

```sh
python results/gaia24dbq/replay.py
python results/gaia24dbq/replay.py --live
```

The first command reads `evidence.json`. The live run fetches the full eRO-ExTra and Gaia alert tables, repeats the sky search, retrieves the full optical curve and rebuilds the target/field infrared series. Inputs carry SHA-256 digests in the packet. VizieR includes a query timestamp in its response header and live catalogs can change. The numerical products are additionally available offline. Each Gaia record gives the reported date, JD(TCB) and G magnitude or quality status. Each infrared record gives the brick-local object ID, band, MJD, flux in AB nanomaggies, inverse variance, number of exposures, fractional neighbor flux and reduced chi-square. Credit for the underlying discoveries and observations belongs to Grotova et al., the eROSITA team, the Gaia Photometric Science Alerts team, ESA Gaia/DPAC, WISE/NEOWISE, unWISE, and the DESI Legacy Imaging Surveys.
