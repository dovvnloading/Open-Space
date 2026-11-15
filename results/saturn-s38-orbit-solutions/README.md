# S/2023 S 38 and the extreme-eccentricity orbit

The nine published positions of Saturn's irregular moon S/2023 S 38 favor a moderate-eccentricity orbit over the extreme orbit in JPL's SAT455 ephemeris. Independent short-arc refits give two-dimensional astrometric RMS errors of **0.192″** and **1.595″**, respectively. The moderate-eccentricity solution predicts the excluded measurements better in **all five leave-one-night-out tests**. This finding bears directly on the claim that this moon reaches the orbital distance of Titan.

Ashton, Gladman, Petit and Alexandersen discovered the moon. Alexandersen's discovery orbit already had moderate eccentricity. The contribution here is to adjudicate the disagreement with SAT455 using the measurements, local orbit refits, seven control moons and held-out observing nights. A literature and source search as of 2026-11-15 found both orbit descriptions in circulation and did not locate this residual comparison. The discovery team explicitly emphasized multiyear arcs in their population analysis. [Discovery circular](https://minorplanetcenter.net/mpec/K25/K25EF4.html) · [Discovery-team study](https://phas.ubc.ca/sites/default/files/RNAAS_preprint.pdf) · [JPL mean elements](https://ssd.jpl.nasa.gov/sats/elem/) · [Published example quoting eccentricity 0.909](https://skyandtelescope.org/astronomy-news/solar-system/a-guide-to-planetary-satellites/)

## What the two solutions imply

| Orbit description | Epoch and elements | Eccentricity | Periapsis (million km) | Discovery-position RMS |
|---|---|---:|---:|---:|
| JPL SAT455 published summary | Mean elements, J2000 reference | 0.909 | 1.17, from mean elements | 1.686″ from the full ephemeris |
| MPC discovery circular | Osculating, 2025-05-05 TT | 0.317 | 14.92 | Circular supplies rounded individual residuals |
| Moderate-eccentricity refit | Osculating, JD 2460200.5 TDB | 0.346 | 14.38 | 0.192″ |
| Extreme-eccentricity refit | Same epoch and force model | 0.915 | 1.03 | 1.595″ |

Titan's orbital semimajor axis is 1.22 million km. At the fit epoch, the two candidate periapses differ by a factor of 14. Comparing predictions on the sky at common times avoids treating mean and osculating elements as interchangeable. Our astrometric comparison tests these two solutions. Another apparition is needed to establish the moon's long-term trajectory. [NASA Saturnian satellite fact sheet](https://nssdc.gsfc.nasa.gov/planetary/factsheet/saturniansatfact.html)

## Measurement and controls

The discovery circular and the MPC's bulk natural-satellite observations each supply the same nine T14 positions over 62.85 days. They are the same measurements. Magnitudes reported on three exposures span 26.6–26.9 in the w band. Each observing night contributes one or two measurements. Archived topocentric positions come from Horizons target 65273, observer T14 (CFHT, Maunakea), quantity 1, astrometric ICRF RA and Dec, sampled at the actual UTC observing times. The frozen Horizons solution is labeled `sat455_merged_DE440`. [MPC bulk archive](https://minorplanetcenter.net/iau/ECS/MPCAT-OBS/SatObs.txt.gz) · [JPL satellite solutions](https://ssd.jpl.nasa.gov/sats/ephem/sep.html) · [Horizons conventions](https://ssd.jpl.nasa.gov/horizons/manual.html)

The seven control satellites are every other S/2023 S 30 through S 41 with discovery measurements in MPEC 2025-E154. Across these controls, the public ephemeris gives 0.553–0.658″ RMS. S 38 gives 1.686″, and its residuals change direction with date. This separates a target-specific discrepancy from a shared time-scale, site or coordinate-frame offset. Each residual is predicted minus observed, in (ΔRA cos Dec, ΔDec). RMS means the square root of the average squared two-dimensional separation.

The independent fit integrates a test moon about the Saturn-system barycenter with Saturn-system monopole gravity and differential accelerations from the Sun and Jupiter system. Positions and velocities from DE440/Horizons drive cubic Hermite interpolation. The model iterates down-leg light time to the CFHT observer, uses TDB for dynamics, and fits six Cartesian parameters with equal weight per coordinate. Four starts cover both signs of depth along the line of sight and two radial velocities per depth. Seeding the same model with the public SAT455 state reproduces the archived nine predicted positions to within 0.00006″. Individual Saturn moons, oblateness and relativistic effects are omitted. The last check validates the short-arc implementation against the official ephemeris, and the long-baseline checks below quantify the increasing difference.

| CFHT night omitted, 2023 UTC start | Positions | Moderate fit, held-out RMS | Extreme fit, held-out RMS |
|---|---:|---:|---:|
| Aug 15 | 2 | 7.32″ | 22.51″ |
| Sep 15 | 2 | 0.316″ | 2.74″ |
| Sep 18 | 2 | 0.286″ | 2.63″ |
| Oct 15 | 1 | 0.278″ | 2.32″ |
| Oct 17 | 2 | 0.344″ | 3.06″ |

Each test drops the entire night, then reoptimizes all six parameters separately around each solution. Dropping August shortens the arc by a month. Its large errors in both branches measure how strongly the earliest night constrains the result. Night exclusions reuse the remaining observations and do not supply independent discovery or recovery epochs.

## Where another observation would discriminate

Illustrative CFHT astrometric ICRF directions at 07:00 UTC are shown below. The last column is the maximum shift of the moderate solution among the five single-night exclusions. It measures sensitivity to omissions and supplies no statistical coverage probability.

| UTC date | Moderate orbit RA, Dec (degrees) | Official SAT455 RA, Dec (degrees) | Separation | Maximum moderate-orbit shift |
|---|---|---|---:|---:|
| 2023-11-16 | 332.41656, −13.19905 | 332.40893, −13.20382 | 0.53′ | 0.03′ |
| 2023-12-15 | 333.78673, −12.65711 | 333.77406, −12.66396 | 0.85′ | 0.12′ |
| 2024-08-15 | 349.94399, −6.87000 | 348.36774, −7.17220 | **95.60′** | 24.62′ |
| 2026-12-01 | 7.86594, +0.59315 | 7.65273, +0.63922 | 13.09′ | 41.50′ |

The August 2024 positions give the clearest archival discriminator among these sampled dates. For December 2026 the omission sensitivity exceeds the separation of the nominal directions. The simplified model propagated from the official starting state differs from Horizons by 1.24″ in August 2024 and 12.31″ in December 2026. No independent image detection or recovered moon is claimed.

## Reproduce

With Python 3, NumPy and SciPy installed, run `python reproduce.py` in this directory. The program reads the frozen CSV inputs, refits both orbit families, performs all five night exclusions and emits the quantities in `results.json`. Small optimizer-dependent differences in the least significant digits are expected. The coordinates in each CSV are kilometers and km/s for vectors, and degrees for astrometry. The state epoch is JD 2460200.5 TDB.

To refresh from the primary archives in a separate copy, run `python acquire.py`. This checks for new MPC observations and for a change of Horizons solution. Quantities computed under an updated ephemeris need a new comparison. The seven controls and all fit residuals are in the machine-readable results.

The two fitted states and the five moderate-orbit night-deletion states are included in `results.json`. Each gives Cartesian (x, y, z, vx, vy, vz) in ICRF, relative to the Saturn-system barycenter, at JD 2460200.5 TDB. For other dates inside the frozen 2023–2026 integration grid, run, for example, `python predict.py 2460537.7916667`. This fetches the T14 observer vector and prints both nominal directions and the five night-deletion directions.

## Independent CFHT image check

A CADC search for CFHT imaging within an 8° radius of (349°, −6°) from MJD 60492 to 60672, with exposures of at least 50 s, returned 452 distinct MegaPrime or WIRCam observations after removing repeated products. Intersecting the moving orbit predictions with the reported sky polygons found two frames covering either nominal solution. A third intersects one moderate-orbit night-deletion trajectory. Examination of processed detector cutouts confirmed coverage of both nominal directions on 2024-07-31 UTC.

| CADC product | Orbit and detector | Approximate exposure start (UTC) | Mid-exposure RA, Dec (degrees) | Aperture S/N | 5σ aperture depth, CaHK |
|---|---|---|---|---:|---:|
| 3030533p | Moderate, CCD06 | 13:53:08 | 350.68064, −6.50421 | −0.21 | 22.77 |
| 3030535p | Extreme refit, CCD12 | 14:01:17 | 349.17024, −6.78118 | −0.02 | 22.75 |

Each image is a single 200 s CaHK.MP9303 narrowband exposure. Circular apertures of radius 5 pixels (0.93″) use a local sky annulus of 10–16 pixels. A robust MAD estimator over 1,256 and 1,263 nearby apertures gives the local noise. Magnitudes use the Elixir camera-run zero point and airmass term without aperture correction. These depth estimates and the discovery w-band magnitudes use different passbands. Neither nominal position gives a significant source, and the shallow narrowband images supply no independent orbit discrimination. [3030533 product](https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/caom2/?Observation.collection=CFHT&Observation.observationID=3030533) · [3030535 product](https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/caom2/?Observation.collection=CFHT&Observation.observationID=3030535)

Run `python check_cfht_images.py` with Astropy, NumPy and SciPy to request the original detector cutouts and repeat these measurements. The numerical results and SHA-256 digests of the fetched products are in `cfht_image_checks.json`.
