# Orbit fit to the three observed apparitions of 2025 UQ146

The Minor Planet Center publishes 15 optical positions on five observing nights. DECam at Cerro Tololo measured five positions in 2014, Pan-STARRS 1 measured four in 2021 and WFST at Lenghu measured six in 2025. A fit initialized from the sky positions alone gives a radial RMS of 0.407 arcsec and a maximum residual of 0.828 arcsec. Propagating that state gives negative Jupiter-relative two-body energy on both targeted passages. These dates are model outputs.

| Nominal solution | Closest approach in 1985 (au) | Energy there (km²/s²) | Closest approach in 2070 (au) | Energy there (km²/s²) |
| --- | ---: | ---: | ---: | ---: |
| Refit from positions, with DE440s | 0.237873 | −0.92259 | 0.155322 | −1.00929 |
| Published MPC state, propagated with DE440s | 0.237873 | −0.92260 | 0.155322 | −1.00929 |

Between the zero crossings in Jupiter-relative energy, the refit gives 1,156.90 days from January 1984 to March 1987 and 881.65 days from April 2069 to September 2071. These durations include times outside Jupiter’s Hill radius. The propagator uses the barycenter and total mass of the Jovian system.

## Prediction checks

For each WFST test, two positions from one October 2025 night were withheld before orbit initialization. The fit used the other 13 sky positions across four nights. The largest error among the six withheld positions was 0.417 arcsec.

| Withheld WFST night (UTC) | Radial error for first position (arcsec) | Radial error for second position (arcsec) |
| --- | ---: | ---: |
| 2025-10-27 | 0.292 | 0.315 |
| 2025-10-28 | 0.123 | 0.160 |
| 2025-10-29 | 0.277 | 0.417 |

The nominal refit and single-night tests each construct their own preliminary orbit using the two latest WFST nights available to the fit. The sky positions on those nights give the direction and angular velocity. A 300 × 211 grid spans topocentric range (0.5–9 au) and range rate (−0.015 to +0.015 au/day). Bound heliocentric trial orbits propagate under solar gravity to the first retained position on each earlier night. The grid ranks the orbits by the summed squared chord distance between predicted and measured unit directions. The best trial initializes the N-body fit. DE440s supplies solar and observer positions. The 13 retained positions determine each single-night prediction.

The whole-apparition checks start at the published MPC Cartesian state, which uses all three apparitions. Withholding 2014 gives an RMS sky-position error of 1,051.5 arcsec. Withholding 2021 gives 547.2 arcsec. A fit to just 2014 and 2021, each a single night, reached the 400-evaluation limit. These checks measure sensitivity within the catalog solution branch. Both nominal orbits in the table above use the same 15 observations.

## Replay

Use Python with NumPy, Pandas, SciPy, Astropy and jplephem. Download the 32 MB JPL planetary kernel and run the calculation in this directory.

```sh
curl -fL https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/de440s.bsp -o de440s.bsp
python fit.py de440s.bsp replayed
```

The script verifies the kernel digest. The refit propagates six barycentric ICRF state components and their derivatives with DOP853. The acceleration includes the Sun and the eight planetary systems, with Earth and Moon represented individually. DE440s position and velocity samples on a two-day grid define cubic Hermite splines. UTC reception times, Earth rotation, observer offsets and iterated geometric light travel time give topocentric astrometric directions. The fit uses tangent-plane angular residuals with a soft-L1 loss and 0.6 arcsec transition scale. The long-span integrations then apply the same force model to the two nominal solutions. Broad time windows bracket each passage, a daily grid brackets the energy sign changes, and scalar roots refine them. This numerical orbit check uses independently written replay logic and sky-fit methods adapted from the repository’s asteroid-pair orbit fit at commit `6ba182a81dd7fa010bff05f7d32a2508b51c68f7`.

Input URLs, request bodies and hashes are in [sources.json](data/sources.json). [Per-observation errors](measurements/fit-residuals.csv) and [states and crossings](measurements/fit-summary.json) are included.
