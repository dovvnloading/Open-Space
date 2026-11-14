# Orbit fit to the three observed apparitions of 2025 UQ146

The Minor Planet Center publishes 15 optical positions on five observing nights. DECam at Cerro Tololo measured five positions in 2014, Pan-STARRS 1 measured four in 2021 and WFST at Lenghu measured six in 2025. An independent fit to these positions gives a radial RMS of 0.407 arcsec and a maximum residual of 0.828 arcsec. The resulting state gives negative Jupiter-relative two-body energy on both targeted passages. These dates are model outputs.

| Nominal solution | Closest approach in 1985 (au) | Energy there (km²/s²) | Closest approach in 2070 (au) | Energy there (km²/s²) |
| --- | ---: | ---: | ---: | ---: |
| Refit from positions, with DE440s | 0.237874 | −0.92259 | 0.155322 | −1.00929 |
| Published MPC state, propagated with DE440s | 0.237873 | −0.92260 | 0.155322 | −1.00929 |

Between the zero crossings in Jupiter-relative energy, the refit gives 1,156.90 days from January 1984 to March 1987 and 881.65 days from April 2069 to September 2071. These durations include times outside Jupiter’s Hill radius. The propagator uses the barycenter and total mass of the Jovian system.

## Prediction checks

The three single-night exclusions in the WFST apparition give maximum absolute sky-position errors of 0.315, 0.160 and 0.417 arcsec. Whole-apparition exclusions are much less predictive. Withholding 2014 gives an RMS sky-position error of 1,051.5 arcsec. Withholding 2021 gives 547.2 arcsec. A fit to just 2014 and 2021, each a single night, reached the 400-evaluation limit. All fits start from the MPC Cartesian solution using all three apparitions. These tests measure sensitivity to the excluded observations within that solution branch. The two nominal orbits above share the same 15 observations.

## Replay

Use Python with NumPy, Pandas, SciPy, Astropy and jplephem. Download the 32 MB JPL planetary kernel and run the calculation in this directory.

```sh
curl -fL https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/de440s.bsp -o de440s.bsp
python fit.py de440s.bsp replayed
```

The script verifies the kernel digest. The refit propagates six barycentric ICRF state components and their derivatives with DOP853. The acceleration includes the Sun and the eight planetary systems, with Earth and Moon represented individually. DE440s position and velocity samples on a two-day grid define cubic Hermite splines. UTC reception times, Earth rotation, observer offsets and iterated geometric light travel time give topocentric astrometric directions. The fit uses tangent-plane angular residuals with a soft-L1 loss and 0.6 arcsec transition scale. The long-span integrations then apply the same force model to the two nominal solutions. Broad time windows bracket each passage, a daily grid brackets the energy sign changes, and scalar roots refine them. This numerical orbit check uses independently written replay logic and sky-fit methods adapted from the repository’s asteroid-pair orbit fit at commit `6ba182a81dd7fa010bff05f7d32a2508b51c68f7`.

Input URLs, request bodies and hashes are in [sources.json](data/sources.json). [Per-observation errors](measurements/fit-residuals.csv) and [states and crossings](measurements/fit-summary.json) are included.
