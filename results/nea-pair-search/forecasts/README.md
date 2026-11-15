# Topocentric recovery pointings

The five gravity-only epochs in `followup.json` use the full two-apparition state for each asteroid and the station coordinates in `stations.txt`. The observations establishing the link end before these epochs. The corresponding old- and new-designation JPL Horizons positions are frozen public single-designation predictions retrieved on 2026 November 15. The offsets use spherical great-circle distances. The H–G estimates assume H=26.0 for VV15, H=27.8 for RG39 and G=0.15 for both. Rotation, photometric calibration and lightcurve amplitude affect actual brightness.

Run these from the results directory, replacing the build path with the directory containing Find_Orb and DE440. The orbit packet's `build.sh` obtains the pinned source and checks the DE440 hash. After that `replay.py` needs Python 3, NumPy and the local fit records.

```bash
bash orbits/build.sh /path/to/orbit-build
python forecasts/replay.py /path/to/orbit-build
```

All times are UTC. Pointings give topocentric J2000 astrometric RA and Dec. Observatory 950 is La Palma, 568 is Maunakea, and G96 is the University of Arizona Mt. Lemmon Survey. Site codes specify a location and do not imply telescope availability. The site constants come from the MPC list. [Observatory codes](https://minorplanetcenter.net/iau/lists/ObsCodes.html)

Each model evaluation propagates Mercury through Pluto and the Moon with DE440 using Find_Orb. For the gravity-only forecasts the uncertainty calculation propagates plus and minus one standard deviation along each of the six eigenvectors of the station-and-UTC-date clustered fitted covariance. Projecting these differences onto the east/north tangent plane produces the local linear two-dimensional covariance. The two quoted axes are square roots of its eigenvalues. For a bivariate Gaussian its one-standard-deviation contour encloses 39.3% probability. A 95% contour would use 2.448 times both axes under the same assumptions. Timing bias, phase-model errors and unmodelled nongravitational acceleration can exceed these statistical errors, especially in 2038. The code rechecks all five nominal positions and both axes against the frozen values.

## VE138 force-model sensitivity

The alternative force fits allow radial acceleration A1 and transverse acceleration A2, both proportional to inverse heliocentric distance squared. The designation linkage holds across these fits. The fitted outward acceleration has several possible causes, including astrometric systematics. The predicted positions under four selected fits appear in `force_forecasts.json`. `A1` fits only radial acceleration, `A2` only transverse acceleration, `A12` fits both, and `A12_site` omits the J95 and 595 stations from the joint fit. The listed range is the spread of those specific models. It does not carry a coverage probability.

At La Palma on 2026 November 15 at 19:30 UTC the full A1+A2 model gives RA 283.165928°, Dec +42.393725°, 2.94 arcsec from gravity only. The estimated brightness is V=20.8 and the target moves at 48 arcsec/minute. An astrometric time bias of 1 s corresponds to 0.8 arcsec along the track. Independent, accurately midpoint-timed subarcsecond astrometry over multiple sites/nights would test these predictions. The models share historical observations, and their centroids differ by about 1.73 gravity-model error-ellipse units at this epoch. A single measurement has limited power to distinguish them. By November 18 the nominal separation falls to 2.10 arcsec and the target fades to V≈24.

At Maunakea on 2038 April 28 at 08:00 UTC the full A1+A2 and gravity-only centers differ by 274.19 arcsec (4.57 arcmin). The fitted A1-only and withheld-site A1+A2 offsets are 270.65 and 282.21 arcsec. Under the full A1+A2 model the linearized 1-sigma semiaxes are 18.9 and 1.15 arcsec. They condition on that force law persisting to 2038. Any planning use must cover model-dependent positions or update the fit with subsequent astrometry.

Reproduce these positions, displacements and full-model covariance using the orbit packet's force-enabled build and the frozen states in `force_states.json`. `replay_force.py` propagates the eight-parameter clustered covariance in addition to the four point fits.

```bash
bash pressure/build_acceleration.sh /path/to/orbit-build
python forecasts/replay_force.py /path/to/orbit-build
```

## Apparition search

A geocentric two-day scan from November 15, 2026 to February 2047 with a solar elongation cut of 60° identified an additional VV15 opportunity in late April 2038. The cataloged Maunakea pointing gives one visible night within this apparition. Sparse sampling can miss short bright intervals. The November 2026 local pointings provide nearer-term opportunities for both bodies. At La Palma, the fitted instantaneous speed for VV15 is 48 arcsec per minute.

To regenerate a Horizons comparison, use its observer mode with the designation `DES=2026 VV15;` or its listed partner, the corresponding observatory code, topocentric `TLIST` Julian date (UTC), `QUANTITIES=1,9,20,23,36`, `ANG_FORMAT=DEG`, `EXTRA_PREC=YES` and `CSV_FORMAT=YES`. Live positions and uncertainties can change when the archive adds observations or accepts identifications. [Horizons API](https://ssd-api.jpl.nasa.gov/doc/horizons.html)

## Hypothetical Goldstone radar discrimination

JPL's automated detectability calculation for 2026 VV15 (solution 2) lists a DSS-14 visibility interval on November 15, 16:54–23:11 UTC and an estimated integrated SNR of 30 for its assumed transmitter, antenna, spin and diameter. The operations schedule records DSS-14's shutdown for modernization from March 2026 until October 2028. The following is a model comparison at the site. [JPL automated detectability table](https://echo.jpl.nasa.gov/data/NEA_SNR/tracks) · [Goldstone operations schedule](https://echo.jpl.nasa.gov/data/schedule.json)

At 20:02:30 UTC the A1+A2 prediction gives a round-trip light time 2.249 ms shorter than gravity only (337.1 km in delay-equivalent one-way range) and an echo frequency 99.52 Hz higher at 8.56 GHz. The model-conditional one-standard-deviation errors in delay and Doppler are 168 μs and 6.86 Hz for A1+A2, and 151 μs and 1.67 Hz for gravity. The A1-only prediction differs from gravity by −2.168 ms and +95.78 Hz at that same epoch. The model centers are well separated in radar space. These are hypothetical center-of-mass returns. No radar measurements entered the fits.

The native Find_Orb radar routine takes the UTC reception epoch, models the light-time-corrected asteroid position at the bounce, and iterates the transmitting station position at the retarded transmit time. Both legs use MPC 253, Goldstone DSS-14. Reported Doppler is modeled received minus transmitted frequency, with positive differences for the force model's higher echo frequency. The calculation excludes instrumental delays, propagation through terrestrial media and the uncertainty of inferring a center of mass from a rotating body's echo. The two-minute centered finite difference of two-way delay verifies the Doppler sign and scale. Values and projected uncertainties at the start, midpoint and end are in `radar.json`.

Starting with the force-enabled build above, run:

```bash
python forecasts/build_radar.py /path/to/orbit-build /path/to/radar-build
python forecasts/replay_radar.py /path/to/orbit-build /path/to/radar-build
```

The radar builder links a second executable in the supplied output directory. It calls the pinned Find_Orb `compute_radar_info` and uses two-line MPC radar metadata for station and frequency. The dummy observation record supplies the projection geometry only.
