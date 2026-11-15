# Umbriel crosses Ariel and its sunlight three times in 22 hours

**Testable prediction for 2050 March 12–13.** Two JPL orbit solutions predict three separate passages in which Umbriel both occults Ariel as seen from Earth and casts a shadow on Ariel. The second and third passages last about 110 and 93 minutes in the current solution. Between them, both disks separate and Ariel clears Umbriel’s shadow. Timing the recovery and the return to overlap would discriminate between orbital solutions. These are future model predictions. A targeted literature search did not establish priority.

## Earth-reception times

Dates below are **geocentric reception TDB**, rounded to a minute, on a 24-hour clock. The delay from Uranus to Earth is about 144 minutes. A 2050 UTC observing schedule will depend on the future TDB−UTC offset. Eclipse contacts are where the penumbra first or last touches Ariel. Projected disk coverage is the fraction of Ariel hidden by Umbriel at minimum disk separation. It is independent of albedo.

| Date in March 2050 | Occultation interval | Occultation minimum | Fraction of Ariel's disk behind Umbriel | Eclipse interval | Eclipse minimum |
|---|---|---|---:|---|---|
| 12 | 04:02–04:18 | 04:10 | 45% | 04:03–04:15 | 04:09 |
| 12 | 21:38–23:27 | 22:20 | 66% | 21:44–23:25 | 22:23 |
| 13 | 00:51–02:25 | 01:50 | 53% | 00:57–02:14 | 01:44 |

The center of Ariel reaches 527, 320 and 450 km from the center of Umbriel’s apparent disk in the three passages. The corresponding distances from the center of its solar shadow are 950, 621 and 820 km. The sum of the satellite radii is about 1,164 km. All three passages produce partial occultations and place part of Ariel inside the umbra. Eclipse and occultation affect overlapping parts of Ariel’s disk. Their depths cannot be added when predicting the combined light curve. Uranus is about 176° from the Sun and Ariel is approximately 14.8, 7.8 and 11.9 arcseconds from the planet's center at the three occultation minima.

## Discovery screen and checks

- **Selection.** With the five major moons, screen every ordered pair for mutual solar eclipses and Earth occultations. The 2048–2051 coarse screen found 114 eclipse and 85 occultation intervals. Only three same-pair overlaps occurred in this screen. All were Umbriel over Ariel. The coarse calculation used three-hour Hermite-interpolated states and ten-minute minimum bracketing. Repeating all searches with five-minute and two-minute brackets recovered all 199 intervals and the same three overlaps. These counts describe the implemented screen and are sensitive to grazing events. [Provisional overlapping intervals](screen.csv) · [Online screen](screen.py)
- **Refinement.** Ten-minute barycentric position and velocity vectors are interpolated to the epoch when light leaves Ariel. The calculation solves separately for the Sun at emission, Umbriel when the ray crosses its plane, and Earth at reception. The finite solar radius sets the penumbra width. Exact circle intersections give occulted disk areas in the tangent plane. A 60-second bracket search and scalar root refinement resolve separate intervals. [All contacts, durations and separations](events.json) · [Frozen vector queries and hashes](horizons.json)
- **Second dynamical solution.** Replace both moon trajectories from Horizons (URA184) with the JPL URA111xl kernels. Reuse the same Sun, Earth and Uranian barycenter trajectories. URA111 reproduces all three combined passages. Occultation minima differ by 46 seconds, 4.6 minutes and 3.7 minutes. The late *occultation-only clear interval* changes from 83.8 minutes to 51.6 minutes. Counting the eclipse as well, the interval entirely clear changes from 83.8 to 50.7 minutes. The two JPL solutions share observations. Their spread is a model comparison, not a calibrated confidence interval. [JPL kernels](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/) · [JPL ephemeris descriptions](https://ssd.jpl.nasa.gov/sats/ephem/)
- **Voyager-era countermodel.** The analytical GUST86 theory, fitted to observations through 1986, predicts the first combined event and a single continuous late occultation from 22:14 to 01:53 TDB, spanning the recovery interval predicted by both modern solutions. Its late eclipse is also continuous (22:27–01:34). This older model does not support three separate passages. Continuous photometry through the gap would directly distinguish the predicted event patterns. [GUST86 replay](gust86_check.py) · [All GUST86 intervals](gust86.json) · [Laskar & Jacobson (1987)](https://ntrs.nasa.gov/citations/19880037204)
- **Historical timing control.** The same method recovers the observed 2007 May 4 occultation of Umbriel by Oberon at 19:09:42.6 UTC, 9.4 seconds early relative to the published 19:09:52 ± 4 s. The predicted impact parameter is 535 km versus 500 ± 80 km measured. The duration is 349 s versus 352 ± 10 s. Historical mutual-event data may contribute to the fit of the modern ephemeris. This checks geometry and time conversion rather than independent predictive accuracy. [Hidas, Christou & Brown (2008), Table 1](https://arxiv.org/abs/0711.2095)
- **Public literature scope.** Searches for the month, both satellite names, and Uranus mutual-event predictions found published tables for 2006–2010 and papers anticipating the 2050 season. There is precedent for simultaneous mutual eclipses and occultations in other systems. No date-specific identification of this three-passage sequence was located in the sources reviewed. Additional ephemerides or unpublished predictions could contain it. [2006–2010 Uranian mutual-event service](https://www.sai.msu.ru/neb/nss/EveUra/EveUra.htm) · [2019 Uranian satellite astrometry paper](https://academic.oup.com/mnras/article/490/3/3464/5610645) · [Europa–Ganymede observed simultaneous event](https://doi.org/10.1016/j.icarus.2022.115348)

The model uses spherical mean radii of 578.9 km for Ariel, 584.7 km for Umbriel, 761.4 km for Oberon, and 695,700 km for the Sun. Limb topography, rotational figure, surface brightness patterns, atmospheric seeing and topocentric reception geometry are outside the model. Exact millisecond timestamps in the JSON record the numerical replay. Minute-scale differences between recent solutions are material.

## Replay

Run from this directory with Python 3, NumPy, SciPy and jplephem. The contact calculation uses only the frozen inputs. It verifies the decompressed SHA-256 of each Horizons vector table. The screen can be rerun separately and requires Horizons access.

```sh
python replay.py --output /tmp/uranian-events.json
diff events.json /tmp/uranian-events.json
python screen.py > /tmp/uranian-screen.csv
```

To refresh the archived Horizons inputs, run `python fetch_vectors.py`. The URA111 files are five-day excerpts of the official JPL binary kernels. They can be regenerated using the commands below. A later jplephem release can change the excerpt's embedded comments while preserving numerical data.

```sh
for body in 701 702; do
    python -m jplephem excerpt --targets "$body" 2050/3/10 2050/3/15 \
        "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/ura111xl-$body.bsp" \
        "inputs/ura111xl-$body.bsp"
done
```

For the optional GUST86 check, compile the specified Stellarium implementation and run the same event finder. The wrapper rotates its planetocentric ecliptic coordinates into the Horizons frame, samples the analytical model every five minutes, and interpolates position and velocity onto the light path. The other bodies keep their original Horizons vectors.

```sh
revision=03ce048ef38f5afa0b9af4a6af674724d5a7713f
mkdir -p gust86-src
for file in gust86.c gust86.h elliptic_to_rectangular.c elliptic_to_rectangular.h calc_interpolated_elements.c calc_interpolated_elements.h; do
    curl -fsSL "https://raw.githubusercontent.com/Stellarium/stellarium/$revision/src/core/planetsephems/$file" -o "gust86-src/$file"
done
cc -O2 -fPIC -shared gust86-src/*.c -lm -o gust86-src/gust86.so
python gust86_check.py gust86-src/gust86.so
```
