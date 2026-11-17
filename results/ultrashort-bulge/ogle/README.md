# OGLE-IV EWS • 2025–2026 bulge alerts

Four of 2,429 bulge alerts have a preliminary EWS timescale below 0.5 day. None of the four yields a robust isolated finite-source lens or an angular Einstein radius from the public OGLE data. Three illustrate how the catalog threshold can select an asymmetric transient, a one-sided low-amplitude peak, or an event with a longer timescale once blended source flux is allowed.

| OGLE event | EWS τ (days) | Independent local test | Disposition |
|---|---:|---|---|
| 2025-BLG-0867 | 0.237 | 160 points in ±15 d. Binary trial χ² 159.8, asymmetric profile 162.6, FSPL 444.4. The light curve fades by 1.39 mag in 0.041 d. | Symmetric single lens fails. KMT-2025-BLG-1213 is the independent catalog alias; its preliminary KMT timescale is 17.37 d. |
| 2025-BLG-1151 | 0.293 | Eight measurements precede the EWS peak, followed by 0.682 d without data around the peak and decline. A solution with tE < 0.5 d is within Δχ² = 0.43 of the longer FSPL minimum. Imposing ρ ≥ 1 costs 0.39. | A genuine brightening with little leverage on tE or ρ. Two Gaia sources lie within 0.44″. |
| 2025-BLG-1201 | 0.158 | The FSPL minimum has tE = 0.075 d and ρ = 1.32. Restricting ρ ≤ 1 costs Δχ² = 1.90. On the same 42 points χ² is 53.6 (FSPL), 39.4 (exponential flare pulse) and 25.2 (binary lens trial). The last epoch that night is still elevated. | Partial, structured brightening. Missing decay, unresolved source color and alternative fits prevent an isolated-planet interpretation. |
| 2026-BLG-0843 | 0.371 | 132 local measurements, 0.640 d gap across the catalog peak. A blended FSPL fit yields tE = 3.80 d. Forcing tE ≤ 0.5 d costs Δχ² = 180.4. | The KMT-2026-BLG-1569 catalog fit (3.40 d) independently favors days. An eclipsing variable 2.97″ away is a possible contaminating neighbor. |

![Measured I-band observations and event-centered trial fits for the four EWS alerts](short-alerts.svg)

## Inventory and previous work

The ASCII tables in [tables.json](tables.json) contain 1,528 alerts in 2025 and 952 in 2026. Selecting IDs containing `-BLG-` yields 1,500 and 929. The alert positions occupy 174 and 152 distinct field numbers, respectively. The survey lists 208 observed fields. [catalog.csv](catalog.csv) supplies event IDs, original field/CCD, star number, coordinates, EWS peak time, impact parameter and timescale. [photometry.csv](photometry.csv) records the raw SHA-256, URL, retrieval or validation time, point count, HJD extremes, and the number of measurements around the catalog peak and in either baseline. Its median near-peak gap includes day/night gaps if they fall in the window. The two source tables were fetched on 2026-11-17 at 10:29 UTC. The server reported last modifications on 2025-11-15 14:07:37 UTC and 2026-11-04 18:47:22 UTC, respectively.

Both raw table hashes exactly match the November 15 `results/ogle-dark-lenses/sources.json` freeze. Its deposited `data/screen.csv.gz` has SHA-256 `471e675d35bc07776d8a4b4ccfc56c2976dc4a47ecc01b45d9b967ce6a68069c`. Re-filtering its existing 2,429 bulge rows gives 2,427 successful fits and two “few points” records. Its sole sub-half-day PSPL fit is 2025-BLG-1151 (0.45081 d). That pipeline seeded tE above 0.5 d, removed outliers after a parallax fit and targeted long events. The local fits here allow hour-scale seeds and retain bright peak points. No dark-lens search was repeated. The fetcher, I-band reader, flux/error conversion, PSPL magnification and conditional linear flux solve come from that project. [Earlier analysis](../../ogle-dark-lenses/README.md)

## Local fits and contaminant checks

[vet_short.py](vet_short.py) fits all valid OGLE points within a fixed ±15 d of the EWS peak in linear flux, adding a 0.004 mag error floor and rescaling once using the full outside-±30 d baseline. Source and blend flux must each be nonnegative. No peak residual clipping is applied. PSPL and uniform-disk FSPL are optimized with differential evolution and local multi-start refinement. Each restricted fit reoptimizes the other parameters. The asymmetric model has different rise and fall timescales with the PSPL shape; the flare trial has an exponential rise and decay with the decay constrained to last at least as long as the rise. The static finite-source binary search samples 0.35 ≤ s ≤ 3.5 and 0.0003 ≤ q ≤ 1. These fits test plausible competitors. They do not uniquely identify a binary, constrain an arbitrary wide host, or incorporate a source-color measurement. The underlying χ², timings, fluxes, bounds and checks are in [short-fits.json](short-fits.json). Source/blend zero points follow the earlier lens module. The plotted error bars include the fitted noise inflation.

| OGLE event | Valid epochs | Time span (yr) | Points more than 30 d before / after EWS peak | Gaia DR3 matches within 2″ | VVV DR2 matches within 2″ |
|---|---:|---:|---:|---:|---:|
| 2025-BLG-0867 | 4,558 | 8.70 | 3,913 / 310 | 1 | 1 |
| 2025-BLG-1151 | 4,564 | 8.70 | 4,210 / 89 | 2 | 1 |
| 2025-BLG-1201 | 345 | 8.31 | 278 / 15 | 2 | 2 |
| 2026-BLG-0843 | 5,046 | 9.73 | 4,862 / 8 | 2 | 1 |

No outside-±30 d observing night contains three measurements more than five scaled σ brighter than its curve's median baseline. This recurrence check loses sensitivity in gaps and to eclipses, slow variables, single flares, and infrequent dwarf-nova eruptions. Blended map V−I for the matching OGLE catalog star is 1.684, 2.776 and 1.942 for 0867, 1151 and 0843; 1201's downloaded map is empty. These unresolved baseline colors cannot establish achromatic magnification. Gaia/VVV neighbors and OGLE mean colors are positional and photometric diagnostics. A public CADC ObsCore HST/JWST image-footprint query and ESA Euclid ObsCore image-footprint query returned zero rows at all four positions. The records in [counterparts.json](counterparts.json) include exact query URLs, response hashes and returned rows. Archive coverage and catalog flags cannot replace an imaging contrast curve.

[SkyBoT queries](moving_objects.json) use Las Campanas topocentric ephemerides, a 120″ radius, no orbit-error filter, and up to three distinct epochs in the observed brightening (first, highest-significance, last). The known body 2011 WJ136 reached 4.66″ from 0867 late in its bright night; it was 117″ away during an earlier bright exposure, with V ≈ 22.2. For 1201, 1978 UD5 was 12.02″ away at one bright point and 83.15″ away at a later bright point. A nearby catalogued asteroid merits image-level checking. Centroid tracking and original difference images are pending, as are a source-color measurement, quantitative high-resolution imaging limits, and an end-to-end survey efficiency.

## Reproduce

From the repo root with Python, NumPy, SciPy, pyerfa, Matplotlib and VBMicrolensing 5.6.1 installed:

```sh
python results/ultrashort-bulge/ogle/freeze.py /path/to/ogle-cache --jobs 6
python results/ultrashort-bulge/ogle/vet_short.py /path/to/ogle-cache
python results/ultrashort-bulge/ogle/counterparts.py /path/to/counterpart-cache
python results/ultrashort-bulge/ogle/moving_objects.py /path/to/ogle-cache
python results/ultrashort-bulge/ogle/window.py /path/to/ogle-cache
python results/ultrashort-bulge/ogle/plot_short.py /path/to/ogle-cache
```

Public EWS photometry is preliminary. OGLE requests contact before publication using the photometry. [2025 EWS](https://ogle.astrouw.edu.pl/ogle4/ews/2025/ews.html) · [current EWS](https://ogle.astrouw.edu.pl/ogle4/ews/ews.html) · [Udalski et al. 2015](https://acta.astrouw.edu.pl/Vol65/n1/a_65_1_1.html) · [SkyBoT](https://ssp.imcce.fr/webservices/skybot/api/conesearch/)
