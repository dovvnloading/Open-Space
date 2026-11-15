# Topocentric astrometric orbit fits

This folder contains a replay of four proposed designation linkages. The observations are time-stamped ground-based optical positions. The 2007 and 2017 links connect separate apparitions. Two additional links connect tracklets within the same apparition.

## Reproduce the results

A Linux build requires Bash, curl, tar, patch, make and a C++17 compiler. The replay uses Python 3, NumPy and SciPy. Run these commands in this directory.

```sh
bash build.sh /tmp/nea-find-orb-build
python replay.py /tmp/nea-find-orb-build /tmp/nea-orbit-check
```

Use `--pair VE138`, `--pair RB15`, `--pair SA12` or `--pair VO` for a single test. The builder pins Find_Orb, lunar, sat_code and jpl_eph to immutable commits and verifies the SHA-256 digest of JPL's 102 MB DE440 planetary ephemeris. The patch exposes the existing Find_Orb topocentric residual calculation as a streaming forward model. It also permits an isolated temporary working directory. All differential corrections are explicit in `replay.py`. No external orbit solution for a linked pair is an input.

## Inputs and model

- `astrometry.csv` freezes all 260 MPC ADES positions returned for the eight designations. It retains MPC observation IDs, UTC, observatory code, catalog, observing mode, precision fields and photometry. Photometry and reported individual error estimates are archived for context. Equal weights in both sky coordinates make the residuals directly interpretable across eras with different uncertainty reporting. The fit uses every listed position. Counts are 110, 62, 48 and 40.
- `seeds.json` contains the eight distinct SBDB catalog records returned on 2026-11-15. Each fit starts from a six-dimensional heliocentric orbit. For full and site-withheld fits it uses the newer designation's catalog elements. In a site-withheld fit those catalog elements are an optimizer start and may already reflect the site's astrometry. The optimizer freely refits all six coordinates using only the retained rows. Residuals on held sites are a site-withheld refit diagnostic, distinct from a completely prospective observation.
- `stations.txt` contains MPC's fixed-width longitude and geocentric parallax constants for the stations. It includes S11, Serra La Nave, which postdates the pinned software's bundled station list. The coordinates, Earth rotation and light propagation matter for these nearby objects.
- The forward model integrates a heliocentric massless test particle with Find_Orb's gravity model, planetary and lunar perturbations and DE440. Observed epochs are UTC. The seed and fitted states are at JD 2461200.5 TDB, in heliocentric J2000 ecliptic coordinates. Positions are AU and velocities are AU/day. No nongravitational acceleration is fitted. Independent RA and declination uncertainties, star-catalog biases, or small timing offsets can change the precision achievable with the same positions.
- Results include a full common fit, a common fit withholding complete observatories, and two fits using just one designation. The label-specific fits extrapolate the refitted published orbit to the other label's timestamps. Over a long time gap, a single short-arc nominal can miss by degrees even if the two labels identify the same body. The program propagates the fitted state to every observation epoch with Find_Orb’s light-time and observer-geometry calculation. For multi-degree, single-label extrapolations the component norm is only a tangent-coordinate diagnostic.

## Output conventions

Each JSON file gives the state vector and osculating elements for each fit, together with the withheld residual summary where applicable. Each CSV gives the computed-minus-observed east and north offset in arcseconds, explicitly including cos(declination) in the east coordinate. Radial RMS means the square root of the mean squared two-dimensional offset. The `in_fit` column specifies which positions influenced each solution.

Formal covariance uses the residual variance per degree of freedom and the pseudoinverse of the full local Jacobian. An alternative finite-sample sandwich covariance clusters on observatory and UTC date. Both are expressed in the state axes and units above. They condition on the selected linkage, local solution, observations and adopted force model. They do not span alternative revolution counts, radiation pressure, or remaining catalog and timing systematics. The station-withheld fits quantify some of these limitations empirically.

MPC observations were requested with a GET to `https://data.minorplanetcenter.net/api/get-obs` with JSON body `{"desigs":["2007 VE138"],"output_format":["ADES_DF"]}`, replacing the designation for each request. Sources include the [MPC observation API](https://docs.minorplanetcenter.net/mpc-ops-docs/apis/get-obs/), [station list](https://www.minorplanetcenter.net/iau/lists/ObsCodes.html), [JPL SBDB Query API](https://ssd-api.jpl.nasa.gov/doc/sbdb_query.html), [JPL DE440](https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/de440/), and [Project Pluto Find_Orb](https://github.com/Bill-Gray/find_orb). Retrieval date was 2026-11-15 UTC.
