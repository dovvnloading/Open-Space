# CEXJ3W2 — a 77.6-minute orbit ambiguity

**The 19 observations admit both a heliocentrically bound orbit and a strongly unbound orbit.** A bound solution with eccentricity 0.917 has 0.873 arcsec radial RMS. An unbound solution with eccentricity 1.431 and nominal excess speed 27.75 km/s has 0.871 arcsec radial RMS. Their sums of squared residuals differ by only 0.0443 arcsec² across all 19 positions. The astrometry does not establish an interstellar object or a securely hyperbolic heliocentric orbit.

| Angular positions fitted | Example solution | Topocentric range (au) | Topocentric range rate (km/s) | Heliocentric eccentricity | Radial RMS (arcsec) | Per-coordinate RMS (arcsec) |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| All 19 | Bound | 0.008382 | −24.86 | 0.91670 | 0.87268 | 0.61708 |
| All 19 | Unbound, 27.75 km/s excess | 0.014308 | −38.15 | 1.43099 | 0.87134 | 0.61613 |
| Published 15 | Bound | 0.008800 | −25.23 | 0.92898 | 0.38944 | 0.27537 |
| Published 15 | Unbound, 48.37 km/s excess | 0.022160 | −54.36 | 2.10234 | 0.38388 | 0.27144 |

All 110 grid distances yielded converged fits for each mask. On all 19 positions the solutions within one additional arcsec² of the minimum cover 0.003329 to 0.1000 au and eccentricities 0.503 to 14.877. The interval reaches the upper edge of the grid. For the published 15 positions, the bound and unbound examples differ by 0.0645 arcsec² total. The different best-fitting excess speeds after changing four observations illustrate the range–range-rate degeneracy.

## Observational basis

The archived pseudo-MPEC reports 19 optical positions on 2026-10-29, 11 from G96 and eight from I52. Four rows carry the published exclusion mark. The published nominal heliocentric eccentricity is 1.53585 ± 1.66, obtained from 15 of the 19 positions. The endpoint angular speed is approximately 7.88°/day. The two Mount Lemmon stations' geocentric MPC constants give a separation of 26 m, equivalent to 0.0036 arcsec of simultaneous parallax for a source at 0.01 au. During the observations Earth's rotation moved the observer by about 1,820 km. The orbit calculation includes that changing observer position. [Archived astrometry and prior solution](https://www.projectpluto.com/neocp2/mpecs/CEXJ3W2.htm) · [Raw retrieval receipt](../evidence/neocp_history_initial.json.gz)

A second association check fits a quadratic to tangent-plane coordinates at just one station and extrapolates across the 8.110-minute gap to the nearest observation from the other. Using the 15 retained measurements gives a 1.70 arcsec forward miss (G96 to I52) and a 0.51 arcsec backward miss. Using all 19 gives 3.52 and 0.76 arcsec. These checks support continuity of the observed track. Converting them into an accidental-association probability would require a measured local background and the number of association trials. For the retained MPC band-code G photometry, G96 reported a median 20.07 mag (19.65–20.51, nine rows), and I52 a median 19.94 mag (19.67–20.24, six rows). These magnitudes come from the archived submissions. [Station holdout calculation](continuity.py)

## Independent range profile

`profile.py` fits the original angular positions directly, using neither Find_Orb nor the published orbit as an initial state. At each of 110 logarithmically spaced reference-epoch topocentric ranges from 0.0005 to 0.10 au, it optimizes two angular offsets, two angular rates and the topocentric range rate. Each observation receives equal weight per tangent-plane coordinate. The first fit uses all observations; a second run uses the exact 15 rows retained by the archived pseudo-MPEC. The range-rate search is bounded at ±400 km/s and tightens at small distances to keep the initial line-of-sight range positive over the arc. This is a profile over the specified grid, with no prior over distances, velocities or object type.

The short-arc calculation integrates Earth-centered inertial motion with Earth gravity plus solar and lunar differential acceleration. ERFA epv00 and moon98 supply the Earth and Moon state, and the model uses each station's geocentric parallax constants, Earth rotation, precession–nutation, two-day IERS polar motion and UT1−UTC values, and iterated light travel time. Differential rotation, acceleration and light travel are retained in every fit. Equatorial residuals are projected onto east and north at the arc midpoint. Step refinement checks motion integration over this short interval. Heliocentric energy and eccentricity are calculated at the mean observation epoch using the fitted state. They describe the local osculating elements during an Earth encounter.

The profile's total squared residual is the sum of squared east and north angular residuals in arcsec². Radial RMS uses one denominator per observation. A per-coordinate RMS uses two and therefore equals radial RMS divided by √2. Fixed differences in summed squared residuals are descriptive; the sources do not provide a calibrated independent error or intra-station covariance for this profile. These fits establish admissible orbital alternatives and do not measure a prior probability of interstellar origin or a false-tracklet probability.

## Reproduction

Use Python 3.12 with NumPy 2.2.6, SciPy 1.15.3 and pyerfa 2.0.1.5. Run `python profile.py` from this directory. The input hash covers the 19 frozen 80-column lines with one LF per record. The pinned IERS rows used by the calculation are archived in `evidence/eop.json`. The grid table records convergence, range rate and RMS for every attempted distance. The representative state vectors and every per-observation residual are in the JSON result. [Profile code](profile.py) · [All range fits](evidence/ranges.csv) · [State vectors and residuals](evidence/profile.json) · [Read-only execution](https://github.com/dovvnloading/Open-Space/actions/runs/44030941200)
