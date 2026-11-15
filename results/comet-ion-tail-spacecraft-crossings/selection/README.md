# Search records and ephemerides

The dated selection is preserved in `screen.csv` (125 approximate leads), `refinement-status.csv` (all 80 refinement requests) and `refined.csv` (one row per retained trial speed for each successfully modeled lead). Dates in these geometry tables are **TDB**, referenced to the J2000 ecliptic with the Sun at the origin. The in situ instrument tables in the parent investigation use UTC. Historical/future positions for the rows called `vex` and `maven` use the planet centers, not the orbiting spacecraft. `earth` is also a planet center. There are no corresponding instrument detections attached to those rows.

`net.py` specifies JPL SBDB selection fields, trajectory identifiers and date bounds. `screen.py` propagates the rounded two-body elements locally and calculates an advected position. Trial emission true anomalies are spaced at 0.009 radians; this can miss narrow encounters. `refine.py` uses 6-hourly Horizons comet vectors, daily observer vectors interpolated with their velocities, and jointly consistent emission and receipt radii for assumed radial speeds between 250 and 900 km/s. The objective is the transverse distance on the spacecraft's heliocentric sphere. It first searches a +/-42-day neighborhood of the coarse guess, then refines the grid minimum. Very small or boundary radial distances are geometrical approximations. There is no plasma tracer in these three tables. Conditions such as comet activity, nongravitational-orbit uncertainty, and wind changes along the trajectory need individual investigation. UCL's Tailcatcher had already established this type of search. [Tailcatcher](https://www.ucl.ac.uk/mathematical-physical-sciences/tailcatcher-service)

Run the following in this directory with Python, NumPy and SciPy. Network access to JPL is needed for a fresh selection. Numerical results for the primary Wind event have a separate replay. The downloads are cached, and new survey tables appear under `cache/`. `survey-inputs.csv` records hashes of responses from this run. JPL may update solutions and response headers. Failed requests or changed orbital solutions can change a fresh search. The signed positional agreement of two Wind instruments and the 2011 rank are evaluated independently in the event folder.

```sh
python net.py
python screen.py
python refine.py
python provenance.py
```
