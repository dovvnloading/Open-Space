# A nine-hour total eclipse of Sinope by Jupiter in October 2037

Jupiter will block the entire solar disc as seen from its distant retrograde moon **Sinope** for about **9 h 27 min on 17–18 October 2037**. At mid-eclipse the moon will lie **21.5 arcminutes from Jupiter** in Earth’s sky, with a nominal uneclipsed visual magnitude of **18.89**. Its large apparent separation from Jupiter during totality makes a ground-based light curve possible. At mid-eclipse the geocentric ICRF position is approximately **07h 58m 33.1s, +20° 47′ 41″**.

| Predicted event | Earth reception, UTC |
|---|---|
| Penumbra begins | Oct 17, 09:20 |
| Totality begins | Oct 17, 15:02 |
| Deepest passage | Oct 17, 19:52 |
| Totality ends | Oct 18, 00:28 |
| Penumbra ends | Oct 18, 06:01 |

The times include light travel from Sinope to Earth. The outer phases last several hours because the Sun has a finite apparent diameter at Jupiter. The light curve gives the direct solar illumination of Sinope. Scattered and refracted light from Jupiter’s atmosphere can contribute residual flux. At 14:12 UTC before totality the geometric model gives 9.1% of ordinary sunlight. By 03:32 UTC after totality it has recovered to 58.9%. [Ten-minute curve](observing-curve.csv) · [Event table](events.csv)

## Sky visibility and tests

At the start of totality Sinope stands **68° above Maunakea**, with the Sun 19° below the horizon. At the end of totality it stands **40° above Kottamia**, with the Sun 46° below the horizon. The solar elongation is 86.7° at mid-eclipse. The relevant observing test is differential photometry at Sinope’s moving position through twilight-safe portions of ingress and egress, with unsaturated comparison stars in the same images. A constant source at the predicted coordinate through both transitions would challenge the identification. [Nine-site geometry](checks.json)

The scan covered eight preselected Jovian irregular moons (Himalia, Elara, Pasiphae, Sinope, Lysithea, Carme, Ananke and Leda) from 2020 through 2040 inclusive. All eight had full ephemeris coverage. **Sinope was the only planetary eclipse in this cohort and interval**, including partial eclipses. A three-hour grid located minima of the angular separation of the Sun and Jupiter; continuous contact solving used Jupiter’s oblate limb and a finite solar disc. Halving the search step reproduced all four contact times to the reported second.

Three earlier JPL orbit fits also put Sinope deep inside the shadow. JUP300 (2013) places the moon **847 km** from the JUP347 prediction at 20:00 TDB on October 17. Substituting that orbit shifts the four eclipse contacts by **5 min 12 s to 5 min 44 s**. JUP343 (2020) and JUP344 (2021) differ by about 190 km and keep every contact within **76 seconds** of the nominal solution. All four solutions predict the full eclipse. They share some underlying astrometry. This comparison isolates the effect of changing Sinope’s orbit under a common illumination model.

The moon lies about **21,200 km** inside the total-eclipse boundary at deepest eclipse. Moving the assumed shadow edge by ±1,000 km changes either totality contact by at most **8.6 minutes**. JPL’s separate observer-table algorithm returns continuous full-eclipse flags from **14:43 to 00:31 UTC** on a one-minute grid. That algorithm uses an equatorial radius, and shares the JPL orbit used by the nominal calculation. A large position error would need to displace the moon beyond the full shadow to eliminate totality. [JPL visibility-code definition](https://ssd.jpl.nasa.gov/horizons/manual.html) · [Observer-query digest](observer-sources.json)

## Data and reproduction

JPL Horizons supplied daily ICRF geometric positions and velocities relative to the solar-system barycenter, in km and km/s (JUP347 for the moons, JUP365 for Jupiter, DE441 for Earth and Sun). Cubic Hermite interpolation preserves both quantities at each daily epoch. At each trial time the calculation solves for the apparent positions of Jupiter and the Sun as seen from Sinope, with distinct light times. Earth reception is computed separately. The one-bar Jovian equatorial and polar radii are 71,492 and 66,854 km. The satellite is treated as a point, the solar surface as uniformly bright and the planet’s limb as opaque. The forecast holds TAI−UTC at 37 seconds after 2017. Source digests cover the numerical rows and exclude Horizons’ changing header timestamp. [Exact archive queries and hashes](sources.json) · [JPL satellite ephemerides](https://ssd.jpl.nasa.gov/sats/ephem/)

Run in this directory with Python, NumPy and SciPy. The first command regenerates the observing curve entirely offline from the included, unrounded nine-day state excerpt. The other commands reproduce the parent search or rebuild the excerpt directly from Horizons.

```sh
python observing_curve.py
python search.py --cache /tmp/sinope-data
python check_events.py --cache /tmp/sinope-data
python observing_curve.py --export-from-cache /tmp/sinope-data
python search.py --cache /tmp/sinope-data --step-hours 1.5 --outdir /tmp/sinope-half-step
```

## Scope and prior work

The eight Jovian targets were in the initial nine-moon search specification committed before ranking. Another thread independently submitted the Saturn result from that interval in [PR #41](https://github.com/dovvnloading/Open-Space/pull/41). This result directory covers Jupiter and its irregular moons.

A November 2026 literature and repository search for “Sinope October 2037 eclipse”, “Jupiter IX eclipse by Jupiter”, “Sinope planetary shadow” and irregular satellite eclipse tables found related stellar-occultation predictions and mutual-eclipse programs for the regular moons. We did not locate a dated Sinope observing report or prediction for this event. Horizons users can already generate its spherical visibility flag. Priority remains open. The new result is a dated, observable event, an oblate-limb contact solution and a reproducible light curve. [Irregular-moon orbit and stellar-occultation study](https://academic.oup.com/mnras/article/462/2/1351/2589643) · [Galilean mutual eclipse program](https://doi.org/10.1093/mnras/stac2538)
