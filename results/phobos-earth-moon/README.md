# Phobos can hide Earth and its Moon together from Gale crater

**At a measured position on Curiosity’s route on Mars, two orbital solutions predict that both Earth and its Moon will disappear behind Phobos for about 9.4 seconds on 5 February 2035.** The Moon disappears first. Earth follows about 2.7 seconds later. When the Moon emerges, Earth stays behind Phobos for roughly another 1.6 seconds. The location is Curiosity’s historical sol-1399 waypoint (RMC 55_2098, 4.68457822° S, 137.3562752° E). This is a forecast for that location. The rover’s 2035 position and operational status are unknown. [Mission waypoint data](https://mars.nasa.gov/mmgis-maps/MSL/Layers/json/MSL_waypoints.json)

| Martian satellite orbit | Both fully hidden, observer UTC on Mars | Duration | Full Moon occultation starts | Full Earth occultation ends |
| --- | --- | ---: | --- | --- |
| JPL MAR099 | 06:24:55.21–06:25:04.64 | 9.43 s | 06:24:52.48 | 06:25:06.25 |
| IMCCE NOE-4-2020 | 06:24:54.00–06:25:03.39 | 9.40 s | 06:24:51.25 | 06:25:04.98 |

The satellite orbits come from different teams and fits. They share historical astrometry. Both calculations also use the same planetary ephemeris, Mars orientation, camera site and Phobos shape. Their 1.2-second difference in predicted onset measures model sensitivity. It does not establish an absolute timing uncertainty. Times use the leap-second table identified below. Future UTC leap seconds can shift the clock labels. [JPL satellite orbit](https://doi.org/10.3847/1538-3881/addab4) · [IMCCE satellite orbit](https://doi.org/10.1051/0004-6361/202039406)

## What an observer would see

Earth and the Moon are separated by 0.03593° (129 arcseconds). They stand about 7.2° above the west horizon at azimuth 273.6°. The Sun is 23.45° below the horizon. The highest terrain on the modeled sightline rises 3.39° at 38.25 km, leaving 3.80° of clearance. The ray follows a great circle out to 100 km over the 128-pixels-per-degree Mars Orbiter Laser Altimeter radius grid. Nearby rocks, dust, optical sensitivity and exposure cadence affect detectability. At 270.6 million km from Mars, Earth subtends about 9.7 arcseconds and the Moon about 2.7 arcseconds. [MOLA grid documentation](https://pds-geosciences.wustl.edu/missions/mgs/megdr.html) · [Measurements for every site](contacts.csv)

Moving the modeled observer 2 or 5 km north, south, east or west preserves the double occultation in all eight sampled directions for both fits. The joint durations range from 6.16 to 9.93 seconds across those samples, and all eight sightlines clear the gridded terrain. At the sol-5073 position, 8.59 km away, the same fits predict full coverage of the Moon alone. Each location uses its own interpolated MOLA radius and a 2 m camera height.

## Search, rejected candidates and prior work

The screen used 52 Gale waypoints and 51 Jezero waypoints from the NASA rover routes. Reading each route in archive order, a waypoint was retained if it lay more than 0.4 km from every previously retained point. The last point of each route was included. From 15 November 2026 to 1 January 2037, a coarse Mars-fixed search sampled the rays every 120 s, kept forward-facing Phobos approaches within 300 km at the route reference point, required Earth above about 4° and the Sun below about −2°, and linearly estimated each site's two impact parameters. Three epochs had a smaller than 13 km minimax separation, all on the Gale route. Linearized timestamps in the screen can differ by tens of seconds from the direct solution. This selection is bounded by its route, angular cuts, and temporal and spatial sampling. [Search code](search.py) · [All sampled sites](search-sites.csv) · [All three leads](screen.csv)

| Geometric candidate at Gale sol 5073 | JPL joint coverage | Elevation of Earth | MOLA terrain | Result |
| --- | ---: | ---: | ---: | --- |
| 27 October 2029 | 8.85 s | 4.66° | 10.63° | Terrain obscures both bodies |
| 26 May 2034 | 7.58 s | 5.80° | 10.12° | Terrain obscures both bodies |

As an observational control, Perseverance’s team recorded Earth passing behind Phobos at Jezero on 2 July 2026. The same calculation at the measured sol-1907 waypoint predicts a roughly four-second full Earth occultation at 10:14 UTC during that sequence. Neither orbit predicts full lunar coverage there. An exploratory stack of the mission’s browse images gave a weak off-Earth peak near the expected lunar position. Temporal subsets varied and a secure lunar identification remains open. The published observation and its discovery belong to the mission team. [NASA image sequence and attribution](https://science.nasa.gov/photojournal/nasas-perseverance-captures-phobos-and-earth/)

Prior searches have cataloged simultaneous solar transits by **Phobos and Deimos**. Curiosity has separately imaged Earth and the Moon and observed Phobos occult Deimos. Searches for the specific combination of Phobos, Earth, Earth’s Moon, Gale crater and the 2035 date found no published prediction of this double occultation through 14 November 2026. Broader unpublished ephemeris work could contain it. [Simultaneous Martian-moon solar transit catalogue](https://arxiv.org/abs/2602.18513) · [Curiosity’s historical imaging](https://www.nasa.gov/history/curiosity-celebrates-10-years-on-mars/)

## Calculation and replay

The final calculation uses converged light times from a rotating surface observer. It tests forward rays from the observer to the centers and 64 points on the apparent circular limbs of Earth and the Moon against the Willner et al. triangular Phobos shape at the light emission epoch. Using Earth’s equatorial radius conservatively encloses its oblate disc. Ray tests omit stellar aberration consistently for every direction. Azimuths and elevations include stellar aberration. Numerical contacts are bracketed on a quarter-second grid and refined to two milliseconds. Doubling the number of limb rays to 128 changes the sol-1399 contact times by at most 0.002 s. Substituting Gaskell’s 3.15-million-triangle shape gives 9.51 seconds of simultaneous coverage for JPL and 9.48 seconds for IMCCE at the reference site. The Gaskell and Willner reconstructions use different methods and share some Viking imagery. The pair measures sensitivity to shape reconstruction. [Shape comparison and source images](https://doi.org/10.1186/s40623-023-01814-7) · [Gaskell shape model](https://sbn.psi.edu/pds/resource/phobosshape.html) · [Alternative shape timings](contacts-gaskell.csv)

Accurate local terrain and Phobos' future position and orientation set the practical limits. [Shape model and methods](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/dsk/satellites/willner_etal_phobos.pdf) · [128-ray check](contacts-128.csv)

Run the commands from this directory with Python 3.11 or newer. The first command below installs the three scientific dependencies into your active environment. Replay downloads approximately 485 MB to the named cache and verifies source fingerprints. The optional Gaskell calculation downloads another 157 MB. JPL and MOLA products have full SHA-256 hashes. The IMCCE BSP has a byte count and SHA-256 hashes for the original header/summary and exact Chebyshev data records and segment trailers used here. Our IMCCE run read those records at their original offsets. The published replay uses the original full file from the IMCCE server. [Input URLs, hashes and site coordinates](sources.json) · [Contact calculation](replay.py)

```sh
python -m pip install numpy scipy spiceypy
python replay.py /tmp/phobos-inputs --download --spatial > /tmp/phobos-contacts.csv
python replay.py /tmp/phobos-inputs --limb-rays 128 > /tmp/phobos-128.csv
python replay.py /tmp/phobos-inputs --download --shape Gaskell > /tmp/phobos-gaskell.csv
OPENBLAS_NUM_THREADS=1 python search.py /tmp/phobos-inputs > /tmp/phobos-screen.csv
diff -u contacts.csv /tmp/phobos-contacts.csv
diff -u contacts-128.csv /tmp/phobos-128.csv
diff -u contacts-gaskell.csv /tmp/phobos-gaskell.csv
diff -u screen.csv /tmp/phobos-screen.csv
```

Offsets in the contact files are seconds from the row’s UTC reference. Empty ingress and egress fields indicate no sampled interval of full-disc coverage within ±30 s. The calculations used SpiceyPy 8.4.0 and CSPICE N0067.
