# Independent orbit check with ESA Aegis

The ESA NEO Coordination Centre's mid-arc equinoctial elements for 2024 XS16 predict the same December 2024 shadow transit. Frozen ESA inputs supply the asteroid orbit and its full six-dimensional covariance. JPL enters this calculation only through the positions of Earth, Moon and other planetary systems. The asteroid's JPL trajectory is not an input. ESA and JPL can share underlying optical astrometry. [Aegis orbit formats](https://neo.ssa.esa.int/objects/help) · [ESA programmatic access](https://neo.ssa.esa.int/en/computer-access)

| Geocentric UTC reception time on 11 Dec | Nominal Aegis calculation |
|---|---:|
| First partial ingress | 18:27:24.8 |
| Start of totality | 18:37:09.8 |
| End of totality | 19:20:47.9 |
| Last partial egress | 19:28:20.2 |

We propagated from the mid-arc epoch (12:00:42.42 TT on the same date) with the heliocentric restricted N-body equations, explicitly separating Earth and Moon. The inferred geocentric direction was also checked against 121 direct Aegis sky predictions at one-minute spacing from 18:00 to 20:00 UTC. The median difference is 0.081 arcsec and the maximum is 0.307 arcsec. Comparison with JPL asteroid vectors over the same interval gives 6.8–23.0 km nominal positional differences. The ESA service reports sky position and a conventional unshadowed magnitude. The shadow fractions below are calculations with its orbit.

| L34 epoch (UTC) | Blocked at midpoint | Mean blocked in 219 s | Approximate 95% orbit-only range in mean blocked |
|---|---:|---:|---:|
| 18:27:41.1 | 0.68% | 3.28% | 3.28–3.57% |
| 18:32:58.9 | 55.84% | 55.84% | 55.84–56.81% |

The two exposure averages predict a **+0.851 mag** differential from shadow alone. This uses a uniformly bright solar disk, opaque spherical Earth of radius 6,371 km and uniform exposure weighting across the nominal 219-second interval. The time convention includes the geocentric asteroid light time. Using the L34 site distance changes that delay by at most 0.022 s. Finite atmosphere, Earth oblateness, limb darkening, stacking cadence, actual exposure weights and asteroid rotation are outside the orbit-only intervals. The interval estimates use 256 scrambled Sobol draws in the local multivariate Gaussian defined by the published element covariance. The point solution need not be the median after transforming the thin, correlated uncertainty ellipsoid to Cartesian coordinates. The intervals do not encompass model or calibration error.

Run `OPENBLAS_NUM_THREADS=1 python shadow.py` here with Python, NumPy 2 and SciPy. It uses the frozen source bytes and prints the four contacts, sky comparison and photometric sensitivity. Planetary states are in au and au/day and are interpolated by cubic Hermite splines. The solver uses the Earth and Moon separately and point-mass direct and indirect terms; relativistic and figure effects are omitted. Source URLs and SHA-256 values are recorded in [sources.json](sources.json).
