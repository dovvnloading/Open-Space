# Tracklet search and recovery

The latest complete 29-day UTC window in the 17 November 2026 ITF snapshot contained zero Rubin (X05) rows and 7,170 Pan-STARRS 1 (F51) rows. Rubin's latest available ITF observations in this snapshot were from July. The snapshot and the separate Rubin submission feed use different publication/ingestion rules. A forced track supplied on three observed nights measures linkage recovery conditional on having the detections.

## Observations and selection

The main window is 2026-10-19 through 2026-11-16 inclusive. Read all MPC1992 rows from the [MPC Isolated Tracklet File](https://www.minorplanetcenter.net/iau/ITF/itf.txt.gz), select stations from columns 78–80, and use the actual UTC observation date. The primary F51 rectangle has 40° ≤ RA < 50° and 10° ≤ Dec < 20°. The all-sky expansion followed inspection of that rectangle's single-night coverage. The Rubin fallback is the 29-day window ending on its last archived observation date. `rubin_feed.py` separately inventories the [Asteroid Institute's Rubin-only MPC export](https://b612.ai/rubin-mpc-downloads/) and reads the most recent completed UTC ingestion-date partition, retaining its observation dates and publication statuses.

Keep optical CCD/CMOS rows with a nonblank submitted name. Form a tracklet by submitted name and UTC day within each station. Deduplicate identical timestamps and fit the unit direction as a linear function of time, using east/north angular motion implicitly on the unit sphere. Require ≥2 distinct observations, a 5-minute to 6-hour span, 0.08–5 deg/day speed, and ≤2 arcsec maximum linear residual. A UTC date split, trail endpoints, nonlinear motion and very short pairs can lower acceptance.

## Linkage and synthetic controls

For every qualifying seed track, combine its topocentric position and angular velocity with 13 logarithmic ranges (0.4–3 au) and 13 radial velocities (−90 to +90 km/s). Add the observing station's heliocentric position and velocity. The ERFA built-in Earth ephemeris and terrestrial rotation supply the observer, using MPC station parallax constants. UT1 is approximated by UTC, TDB by TT, and polar motion is omitted. Compare a solar-gravity Taylor propagation of these Cartesian trial states with all eligible observations on later UTC nights. Require 0.7–3.25 days between consecutive tracklet midpoints and a common trial within 600 arcsec on the second and third nights. The initial angular prefilter and fit criteria are explicit in `link.py`.

Fit all directions in a proposed triple to a six-coordinate solar two-body orbit with numerical DOP853 propagation and light-time correction. Retain fits with ≤1 arcsec angular RMS, ≤2.5 arcsec maximum separation, and nominal heliocentric excess speed >10 km/s. Distinct triples are independent proposals and may share detections. A screen-positive fit would require an independent multi-solution, perturbed original-barycentric refit and identity checks. This short-arc screen alone cannot establish binding energy.

Draw synthetic objects at an actual seed direction and an actual three-night cadence from the same station and interval. The geocentric range is log-uniform over 0.45–2.8 au, the true heliocentric excess speed uniform over 15–60 km/s, and the heliocentric velocity direction isotropic. Use numerical two-body propagation, iterative light time and independent 0.2 arcsec Gaussian tangent-plane noise per coordinate. Force these astrometric detections into the full real background one object at a time. Apply the same angular, timing, residual and fitted-speed cuts as for real triples. A recovery requires all three injected tracklet identities. Count accepted mixed associations separately. In the rectangular subset, reapply its boundaries to every synthetic observation.

Injections assume usable astrometry already exists. Later positions are placed at actual observing times but are not checked against exposed pixels, camera footprints, weather, flux limits or trailing losses. Recovery is conditional on the injection distribution and excludes single-night and two-night objects. Five controls independently rotate the real sky per night by a random RA offset in ±3°, preserving each tracklet's internal geometry while perturbing cross-night alignments. They measure chance-link behavior under this specific control. No footprint-complete or space-density limit follows.

## Reproduction

From this directory, install `numpy==2.2.6 scipy==1.15.3 pyerfa==2.0.1.5` under Python 3.12. Run `python extract.py --output run` and compare hashes with the archived receipts. A future live download can change. With the archived observation slices in `run/`, execute:

```sh
python link.py --data run --survey X05 --trials 100
python link.py --data run --survey F51 --trials 100
python link.py --data run --survey F51-all --trials 100
python rubin_feed.py run
```
