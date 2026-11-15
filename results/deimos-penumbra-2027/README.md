# Phobos casts three penumbral shadows on Deimos in October 2027

Phobos will cross the Sun as seen from Deimos on October 8, 9 and 11, 2027. At each event its entire silhouette fits inside the solar disk simultaneously across the sunward face of Deimos. The first two shadow passages occur about 30 hours 20 minutes apart. Deimos stands high in a dark sky above Gale crater during the first and above Jezero crater during the second. The third occurs below both horizons. These are predictions for the illumination of one moon by the other. No measurement of a brightness dip is presented here.

| UTC at Deimos, JPL closest approach | Solar disk area blocked at Deimos center | Conservative whole-surface interval, JPL / IMCCE | Deimos and Sun altitude at Gale | Deimos and Sun altitude at Jezero |
| --- | ---: | ---: | ---: | ---: |
| Oct 8, 17:29:01 | 3.10% | 91.8 / 91.4 s | +72.1°, −52.2° | +5.3°, +8.2° |
| Oct 9, 23:48:33 | 2.48% | 63.8 / 63.4 s | +13.9°, −45.1° | +65.3°, −64.8° |
| Oct 11, 06:13:26 | 2.13% | 26.0 / 24.8 s | −64.3°, +38.6° | −9.7°, −19.6° |

The closest-approach times from JPL MAR099 and IMCCE NOE-4-2020 agree to 0.42 seconds across all three dates. A sphere with Phobos's approximate mean radius gives first-to-last contact windows of 170, 113 and 80 seconds for the center of Deimos. Mars's limb stays more than 14° from the Sun at each closest approach. The occulting body here is Phobos. The times are UTC at Deimos. Light reaches an observer on Mars about 0.1 seconds later. UTC conversions assume the last leap-second update in the archived kernel.

## Geometry and precision

The screen covers November 15, 2026 through January 1, 2033 (exclusive). It searches both directions of shadowing between the two moons. A 300-second grid locates local minima of the transverse moon-to-moon distance relative to the Sun. Each candidate within 750 km is refined using SPICE converged light-time and stellar-aberration corrections, finite solar and satellite radii, and a clear solar disk above the Mars limb at the receiving moon. The resulting table contains 30 mean-sphere overlaps. Halving the grid step recovers the same 30 events with a maximum 0.001-second change in the rounded closest-approach times. These results do not establish completeness for other periods, eclipses superposed on Mars's shadow, or very shallow grazes of an irregular limb. [Screen](screen.csv)

The detailed replay projects every triangular plate of the Willner Phobos shape in perspective, unions the projected triangles in raster space, and clips them against the Sun. At closest approach, changing the image width from 512 to 1024 pixels changes the blocked area by less than 0.009 percentage points. The whole-surface interval uses the maximum Phobos vertex radius plus 1 meter, a 7.8-km bounding sphere for the IAU Deimos reference ellipsoid, and the combined angular parallax of the Sun and Phobos over those positions. In that interval Phobos's full outline is inside the solar disk even at the least favorable location in the bounding sphere. Local topography at sunrise and sunset is outside the model. [Shape model](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/dsk/satellites/phobos_2014_09_22.cmt) · [SPICE planetary constants](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00011.tpc)

The area fraction is a uniform-solar-disk illumination estimate. Solar limb darkening, bandpass, surface scattering, terrain, sky background, and camera signal-to-noise affect the measured light curve. Elevations use a geometric horizon and two fixed mission waypoint coordinates, MSL sol 5073 (4.82602331° S, 137.388886° E) and Mars 2020 sol 1906 (18.41996953° N, 77.22012875° E). Future rover positions and operations are unknown. The two satellite ephemerides use overlapping historical observations, shared planetary and rotation inputs, and the same shape. Their time spread measures sensitivity to the orbit choice. [Detailed timings and shape measurements](contacts.csv) · [MSL waypoints](https://mars.nasa.gov/mmgis-maps/MSL/Layers/json/MSL_waypoints.json) · [Mars 2020 waypoints](https://mars.nasa.gov/mmgis-maps/M20/Layers/json/M20_waypoints.json)

Curiosity photographed Phobos passing in front of Deimos in 2013. Cody's 2026 calculation catalogs both moons transiting the Sun from a location on Mars. Those geometries answer different observing questions. The literature search did not locate an earlier prediction of this October 2027 sequence of changes to Deimos's illumination. Publication priority remains unestablished. [Curiosity observation](https://science.nasa.gov/resource/before-and-after-occultation-of-deimos-by-phobos/) · [Cody 2026](https://arxiv.org/abs/2602.18513) · [Mars moon ephemeris from Mars Express astrometry](https://doi.org/10.1051/0004-6361/202039406)

## Replay

From this directory, with Python 3.11+, NumPy, SciPy, Pillow and SpiceyPy installed, run the following command. The download verifies SHA-256 hashes and writes a sparse copy of the independently computed French ephemeris containing only the required records. Replay the 150-second screen by adding `--step 150` with a separate output directory. [Kernel URLs, hashes and byte ranges](kernels.json)

```sh
OPENBLAS_NUM_THREADS=1 python replay.py .cache reproduced --download --screen
```
