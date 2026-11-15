# Independent optical observations and later TESS visits

Ground observatories recovered the primary rotation of both asteroids during the 2018 discovery sector. The weaker dimmings remain uncertain in ground data. Later TESS visits give smaller amplitudes or shifting phases when tested at the recurrence measured in Sector 3.

## Ground test

The Minor Planet Center supplied public optical measurements for each object. We omitted station C57 (TESS) and I41 (Palomar ZTF) from the ground test. The ground and spacecraft observations span overlapping dates. JPL Horizons supplied geocentric distances, heliocentric distances, and solar phase angles at one-day intervals. The test linearly interpolates these quantities to the ground observation times. [MPC Observations API](https://docs.minorplanetcenter.net/mpc-ops-docs/apis/get-obs/) · [JPL Horizons API](https://ssd-api.jpl.nasa.gov/doc/horizons.html)

We froze the primary rotation, recurrence, epoch, and folded dimming shape from Sector 3. The primary template uses five rotational harmonics fit to the original TESS light curve along with the dimming template and a quadratic solar phase correction. At the ground observation times we fit both template amplitudes, an offset for each station/filter combination, and a continuous piecewise-linear solar phase correction with a break at seven degrees. Apparent magnitude also includes the standard distance correction of five times the base-ten logarithm of the product of the two distances in au. Times use reported observation JD without an additional target light-time correction. The Earth-spacecraft light-time difference is small compared with the two-hour events. [TSSYS-DR1](https://archive.konkoly.hu/pub/tssys/dr1/)

| Asteroid | Window | Ground measurements | Observing nights | Rotation scale | Dimming scale |
|---|---|---:|---:|---:|---:|
| (29489) | Sector dates | 49 | 9 | 1.09 ± 0.08 | 0.38 ± 0.26 |
| (29489) | Sector dates plus 14 days on each side | 111 | 21 | 1.02 ± 0.11 | 0.55 ± 0.27 |
| (3177) | Sector dates | 57 | 12 | 0.84 ± 0.43 | 1.52 ± 2.17 |
| (3177) | Sector dates plus 14 days on each side | 112 | 21 | 1.14 ± 0.22 | 0.69 ± 0.76 |

A scale of one is the Sector 3 amplitude. Uncertainties are approximate night-clustered sandwich standard errors with a finite-sample correction. The small number of nights limits their accuracy. The estimates for 29489 are suggestive and depend on the window. The 0, 7, 14, and 25-day extensions in ground.csv share data. MPC magnitudes come from heterogeneous filters and often have 0.1 mag precision. The spacecraft recurrence and viewing geometry can evolve between apparitions.

## Later TESS test

The C57 observations are additional TESS visits reported to the MPC. The original sector checks the identification and analysis procedure. Its exposures overlap TSSYS-DR1. Each subsequent visit in the table was obtained at a new apparition. [MPC observatory codes](https://docs.minorplanetcenter.net/mpc-ops-docs/apis/obscodes/)

| Asteroid | Visit begins (JD) | Measurements used | First-half template scale | Second-half scale at the first-half phase | Difference between independently fitted phases (cycles) |
|---|---:|---:|---:|---:|---:|
| (29489) | 2458383.74 (Sector 3) | 1,102 | 0.82 | 0.90 | 0.00 |
| (29489) | 2458842.53 | 720 | 0.55 | 0.22 | 0.04 |
| (29489) | 2460288.96 | 241 | 1.45 | 0.57 | 0.19 |
| (3177) | 2458383.72 (Sector 3) | 1,132 | 0.71 | 1.20 | 0.01 |
| (3177) | 2458899.38 | 946 | 0.25 | 0.19 | 0.09 |
| (3177) | 2459361.78 | 2,975 | 0.80 | −0.06 | 0.34 |

The test groups C57 measurement dates into runs separated by more than ten days and keeps runs with at least 100 observations. Each run removes a quadratic time trend and five harmonics of the primary rotation. We refine the spin over ±0.5% of the Sector 3 value and exclude large photometric outliers. The dimming shape and recurrence stay fixed. The first half selects the phase from a grid spaced by 0.01 cycle. The last column in later_tess.csv measures the second-half amplitude at that fixed phase. The primary spin and trend used the full run, and the first-half phase search biases the first-half amplitude upward. A change in mutual-event visibility, viewing geometry, photometric quality, or synodic period could reduce the fitted amplitude.

## Replay

Run from this directory with Python and NumPy. This reproduces ground.csv and later_tess.csv offline using the compressed measurement snapshots and the parent directory's measurements.csv and fold.csv.

```bash
OPENBLAS_NUM_THREADS=1 python check.py
```

The input snapshot contains the MPC measurements used, the daily Horizons predictions, and the selected Sector 3 rows used to fix the spin template. sources.json gives each original URL, MPC request body, and SHA-256 of the original response. The archived data were retrieved on 2026-11-15. Running `python freeze.py` explicitly refreshes the snapshots from the live sources. MPC can add observations and Horizons recomputes the solution and includes a generation time in its response.
