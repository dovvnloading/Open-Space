# JPL catalogue screen

**Exactly three of 2,290 selected JPL records passed both nominal osculating speed ≥5 km/s and formal heliocentric eccentricity excess ≥5σ. Their persistent IDs identify 3I/ATLAS, 2I/Borisov and 1I/‘Oumuamua.** Five further records reached the speed threshold. All five had eccentricity excess below 0.36σ. The thresholds define a finite catalogue screen. They do not measure how many incoming objects the surveys could have missed.

## Snapshot and replay

JPL returned every record satisfying `e|GE|1` in a single full-precision request without a row limit. The API returned 2,290 unique SPK IDs, including fragments. There were 1,763 exactly parabolic solutions and 527 nominally hyperbolic solutions. Positive finite formal eccentricity uncertainties were available for 537 records. JPL lists nonzero A1, A2, A3, DT or S0 model terms for 51 records. Absence of a fitted term leaves unmodelled activity possible. In the selected population 1,638 records had reported arcs ≤30 days, five had missing arc lengths and 49 lacked an observation count. The most recent last observation in this response is dated 2026-11-09. See [JPL field definitions](https://ssd-api.jpl.nasa.gov/doc/sbdb_query.html).

From the parent directory, `python catalog/audit.py catalog` verifies the bytes and regenerates the complete ranked `catalogue.csv` and `summary.json`. The raw HTTP response is compressed in `sbdb.json.gz`. The SHA-256 in `receipt.json` refers to the decompressed bytes. To fetch a later independent snapshot into an empty directory, run `python catalog/audit.py /path/to/new-snapshot --fetch`.

The score uses `v_osc = sqrt[(GM_sun/AU) (e−1)/q]` with q in au, GM in km³/s² and AU in km. The second column of scores, `z_e = (e−1)/sigma_e`, uses only the reported marginal formal error. These scores share an eccentricity measurement and are correlated. Encounter geometry, barycentric frame transformations, nongravitational forces, correlated residuals and short-arc degeneracies all affect an incoming-energy test. We do not convert this screen into a population upper limit.

## Five other records above 5 km/s

| SPK ID | Object | v_osc (km/s) | z_e | Independent evidence |
|---|---|---:|---:|---|
| 1003934 | C/2023 F2 (SOHO) | 17.556 | 0.024 | Ye et al. identify a Meyer-group near-Sun comet. Their attempted ground recovery covered 5–10% of the predicted sky probability to R ≈ 20. [Paper](https://doi.org/10.3847/2515-5172/acda9b) |
| 1001481 | C/1999 U2 (SOHO) | 17.381 | 0.052 | The Sungrazer Project associates this Marsden-group object with the 2005 return C/2005 W5. [Project cross-identification](https://sungrazer.nrl.navy.mil/periodic_soho) |
| 1001088 | C/1954 O1 (Vozárová) | 9.145 | 0.179 | The JPL fit used four observations. Hasegawa's gravitational fit used nearly 60 positions across 144 days and reported heliocentric e = 1.0000686 ± 0.0000274 and q = 0.6771121 au, equivalent to 0.300 km/s at that epoch. [Hasegawa (1966), pp. 283–284](https://doi.org/10.1093/pasj/18.3.281) |
| 1003178 | C/2012 E2 (SWAN) | 7.913 | 0.354 | The SOHO team classifies it as a Kreutz sungrazer. The JPL solution spans two days. [Project discovery log](https://sungrazer.nrl.navy.mil/sites/sungrazer.local/files/html/legacy_reports/comets201205_arch.html) |
| 1003589 | C/2011 N3 (SOHO) | 6.732 | 0.140 | This sungrazer was directly observed in the solar atmosphere. The JPL arc spans one day. [Schrijver et al. (2012)](https://doi.org/10.1126/science.1211688) |

## Slower objects and prior identification

C/1980 E1 (Bowell) and C/2024 L5 (ATLAS) reached 3.902 and 3.101 km/s in the screen. Published reconstructions find planetary ejection and eventual outbound speeds of about 3.8 and 2.8 km/s respectively. C/1997 P2 is an additional example. Its 2.411 km/s osculating speed follows a 0.652 au approach to Jupiter six months before perihelion. The preferred CODE solution describes a likely Oort-spike comet with substantial planetary perturbation. [Bowell and ATLAS reconstruction](https://doi.org/10.1051/0004-6361/202451920) · [CODE solution for Spacewatch](https://apollo.astro.amu.edu.pl/CODE/orbit.php?int=1997p2a5&orb=osculating)

C/2025 V1 (Borisov) reached 3.850 km/s with z_e = 5.27 on a ten-day arc of 116 JPL observations. It is an identified active comet, distinct from 2I/Borisov, and already had a published hyperbolic solution in its discovery circular. Our catalogue screen prioritizes a test of its incoming energy. A survey-image footprint would establish an opportunity to look for it. A precovery requires an independently identified moving source with a documented astrometric and photometric measurement. [Discovery circular, MPEC 2025-V40](https://minorplanetcenter.net/mpec/K25/K25V40.html)

The JPL subset contains seven asteroid-class rows, including 1I. The newest two, A/2026 T2 and A/2026 T3, have 25- and 35-day arcs respectively with formal osculating eccentricity excesses under 0.4σ. A/2024 U2 has a 721-day arc and e = 1.00266633 ± 0.00001138 at q = 8.1791 au. Its 0.538 km/s osculating speed illustrates why a precise heliocentric eccentricity can still require a barycentric reconstruction.

The C/2019 Y4 and C/2025 K1 components appear as separate JPL records. In particular, the fits to C/2025 K1-B and C/2025 K1-D begin in November after the parent's October 8 perihelion and during or after its observed breakup. These are documented comet fragments. Their measured trajectories do not independently count incoming bodies. [Hubble observations of C/2025 K1](https://science.nasa.gov/asset/hubble/fragmenting-comet-c-2025-k1/)

## Published incoming orbits

CODE listed a preferred osculating model for 301 of the 2,290 JPL entries. For each matched designation we retrieved the corresponding *original barycentric* page and uncertainty. Of those preferred models, 286 have positive reciprocal semimajor axis and 15 have negative values. This comparison uses published orbits. Independent fits to the astrometry are reported separately. Every row and source URL appears in [code-original.csv](code-original.csv). The [retrieval receipts](code-receipts.json) identify each downloaded page by SHA-256 and UTC response date.

| Comet | Preferred original 1/a ± σ (10⁻⁶ au⁻¹) | Formal excess | Checks on the interpretation |
| --- | ---: | ---: | --- |
| C/1995 Y1 (Hyakutake) | −48.77 ± 9.04 | 5.39σ | CODE model n5 gives the negative value. Alternative model m5 gives +55.50 ± 7.44 from the same astrometry. [n5](https://apollo.astro.amu.edu.pl/CODE/orbit.php?int=1995y1n5&orb=original) · [m5](https://apollo.astro.amu.edu.pl/CODE/orbit.php?int=1995y1m5&orb=original) |
| C/2007 W1 (Boattini) | −29.84 ± 1.29 | 23.13σ | Catalogue link bt is strongly negative. The authors favor a distant-data solution d5, which gives −0.09 ± 1.22. [bt](https://apollo.astro.amu.edu.pl/CODE/orbit.php?int=2007w1bt&orb=original) · [d5](https://apollo.astro.amu.edu.pl/CODE/orbit.php?int=2007w1d5&orb=original) |
| C/2020 W5 (Lemmon) | −202.96 ± 59.39 | 3.42σ | CODE’s 56 observations cover a heliocentric distance range of only 3.36–3.38 au. The authors use a gravity-only model of quality class 2b and identify planetary perturbations during passage. [CODE solution description](https://apollo.astro.amu.edu.pl/CODE/orbit.php?int=2020w5a2&orb=osculating) |
| C/2011 M1 (LINEAR) | −1,546.55 ± 506.15 | 3.06σ | CODE used 122 pre-perihelion observations over 41 days. Its authors could not determine nongravitational effects, assign orbit-quality class 3a, and discuss compatibility with an Oort-cloud origin at about three sigma. [CODE solution description](https://apollo.astro.amu.edu.pl/CODE/orbit.php?int=2011m1a5&orb=osculating) |
| C/2012 C2 (Bruenjes) | −2,592.99 ± 883.18 | 2.94σ | Its 264 pre-perihelion positions span 16 days. CODE assigns quality class 3b and reports that nongravitational forces are unmeasurable. [CODE solution description](https://apollo.astro.amu.edu.pl/CODE/orbit.php?int=2012c2a5&orb=osculating) |

The remaining 10 negative nominal values are all within 2.2σ of zero. The matched set has incomplete coverage of the JPL list. A catalogue cross-match cannot establish how many distinct incoming bodies populate the rest of that list.

The four largest formal high-speed SOHO alternatives in the blind batch have prior independently fitted paths and documented LASCO observations. [Spacecraft astrometry, photometry and historical orbit fits](soho/README.md).
