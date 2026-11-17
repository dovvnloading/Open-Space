# Tracklet search and recovery

The all-sky Pan-STARRS 1 search returned zero candidate triples from 7,170 real observations on 17 nights. It recovered 99 of 100 synthetic hyperbolic objects whose three tracklets were supplied at real observing times (exact two-sided 95% interval 94.55–99.97%). Of the 52 objects moving inward at the seed heliocentric epoch, 51 were recovered. The trial distribution covers 0.45–2.8 au in geocentric distance and 15–60 km/s in heliocentric excess speed. Coverage of exposed pixels and object brightness differs from this forced-astrometry experiment.

| F51 or Rubin sample | Real rows | Submitted tracklet-days | Motion-qualified tracklet-days | Real passing triples | Forced recovery |
| --- | ---: | ---: | ---: | ---: | --- |
| All-sky F51, Oct 19–Nov 16 | 7,170 | 2,131 | 2,130 | 0 | 99/100, including 51/52 inbound |
| Fixed F51 rectangle in the same window | 57 | 18 | 18, all on Nov 3 | 0 | No three-night opportunity |
| All-sky X05, Oct 19–Nov 16 | 0 | 0 | 0 | 0 | No observations |
| All-sky X05, Jun 16–Jul 14 | 36 | 8 | 5 on Jul 9 and 12 | 0 | 97/100, including 53/55 inbound |

For Rubin's July control, usable timestamps existed on July 9, 10, 12 and 14 before the motion and within-night linearity cuts. Only two nights retained real tracklets. The 97/100 recovery (91.48–99.38%) describes forced objects sampled at the four pre-cut nights. The F51 all-sky and July Rubin experiments each yielded zero accepted mixed-source triples and zero accepted links in all five rotated-night controls. A matched injection counts as recovered when the three true groups are linked and the fitted nominal excess exceeds 10 km/s. One true 51.59 km/s F51 object fitted at 134.19 km/s, despite recovering its identity.

The separate Rubin MPC replica check retrieved all 29 daily ingestion partitions for Oct 19–Nov 16, totaling 45,440 distinct rows. The statuses were 45,352 accepted/waiting (`p`), 59 published (`P`), and 29 ITF (`I`). Their actual observing dates ranged from May 6, 2025 through July 11, 2026. Zero ITF observations in this export were taken during the primary window. The Nov 16 partition contained a header only. The [query receipt](run/rubin-feed.json) records generation-pinned file URLs, byte counts, SHA-256 digests, observation dates and statuses.

Under an independent Poisson model for the mean number of objects **already supplying three timestamped astrometric groups at the modeled cadence and noise** drawn like the forced F51 population, the conservative 95% upper bound is 3.90 such objects in this slice. This uses 3.6889 = −ln(0.025) for zero detections divided by the one-sided 97.5% lower recovery bound 0.94554. The two error probabilities sum to 0.05. The assumptions exclude scheduling, exposure footprints, photometric completeness, source detection, tracklet formation and the probability of entry in the ITF. A sky density, encounter rate or Solar System population limit would require those factors.

### Exposure-time proximity

For each forced, noisy direction, find the nearest archived real direction from the same station within one second of the supplied observation timestamp. At each of 0.1°, 0.5° and 1.5° require either every sample, or at least two samples per tracklet, on each of the three nights. In **both** 100-object simulation sets, 100 first-night tracklets met every radius by construction. Zero second-night or third-night tracklets had two near-time positions within even 1.5°. Thus zero forced three-night paths in the 99/100 F51 experiment or 97/100 X05 experiment passed this observing-support proxy. A missing ITF witness cannot establish that no telescope imaged the location. Actual chip coverage, an exposure identifier and photometric detectability were not measured. The per-trial distances and exact-time trial rerun are archived. [Read-only run](https://github.com/dovvnloading/Open-Space/actions/runs/44027015100)

### Search using observed sky positions and times

A separate replay drew 20,000 solar-hyperbolic paths with independent seed 317801. Each started at the sky position of one of the 2,131 timing-valid real F51 tracklet-days chosen uniformly with replacement. The distance, excess-speed and velocity-direction distributions matched the forced experiment. For each path, search **all** later real F51 groups in the window. A 2° preliminary midpoint gate proposes tracks before precise propagation. A three-night sequence passes the 1.5° proxy only if every modeled exposure direction is within 1.5° of the **same real tracklet's** submitted direction at that timestamp on each night. Where multiple sequences pass, choose the one with the smallest maximum separation. This choice uses only geometry. Add independent 0.2 arcsec noise, inject the three groups into the original background, and apply the frozen linker. The prespecified fit cap of 100 was never reached. [Frozen replay code](pointing/replay.py) · [Execution log](https://github.com/dovvnloading/Open-Space/actions/runs/44027435511)

| Radius | Paths supported among 20,000 | Exact 95% interval for supported fraction | Supported paths recovered | Exact 95% interval for conditional recovery |
| --- | ---: | ---: | ---: | --- |
| 0.1° | 0 | 0–0.0184% | 0 tested | Undefined |
| 0.5° | 1 | 0.00013–0.0279% | 1/1 | 2.50–100% |
| 1.5° | 37 | 0.1303–0.2549% | **35/37** | **81.81–99.34%** |

At 1.5°, 37/20,000 paths (0.185%) found a three-night geometric witness. Of those, 35 (94.59%) passed all per-track quality cuts and recovered the true association, including 14 of 15 inbound cases. No accepted link mixed real and synthetic identities. Both losses occurred at within-track quality cuts. The model-specific yield from all seed paths was 35/20,000 (0.175%). Thirty-five supported cases used Nov 1–3. The other two used observations during Nov 13–16. Only 24 distinct anchor groups contributed to the 37 cases.

These radii are proxy choices. The first position is forced onto a real track, later positions are checked against recorded directions, and both the scheduling distribution and the first-position distribution are sampled from *isolated tracklets*. A nearby detection does not identify the chip occupied by the synthetic source or establish adequate depth, photometric recovery or inclusion in the ITF. The 35/37 conditional probability and the 99/100 forced-time probability have separate denominators. [Geometry summary](run/pointing-result.json) · [All 37 paths and witness groups](run/pointing-trials.json)

## Observations and selection

The main window is 2026-10-19 through 2026-11-16 inclusive. Read all MPC1992 rows from the [MPC Isolated Tracklet File](https://www.minorplanetcenter.net/iau/ITF/itf.txt.gz), select stations from columns 78–80, and use the actual UTC observation date. The primary F51 rectangle has 40° ≤ RA < 50° and 10° ≤ Dec < 20°. The all-sky expansion followed inspection of that rectangle's single-night coverage. The Rubin fallback is the 29-day window ending on its last archived observation date. `rubin_feed.py` separately inventories the [Asteroid Institute's Rubin-only MPC export](https://b612.ai/rubin-mpc-downloads/) and reads each ingestion-date partition within the same 29 completed UTC days, retaining its observation dates and publication statuses.

Keep optical CCD/CMOS rows with a nonblank submitted name. Form a tracklet by submitted name and UTC day within each station. Deduplicate identical timestamps and fit the unit direction as a linear function of time, using east/north angular motion implicitly on the unit sphere. Require ≥2 distinct observations, a 5-minute to 6-hour span, 0.08–5 deg/day speed, and ≤2 arcsec maximum linear residual. A UTC date split, trail endpoints, nonlinear motion and very short pairs can lower acceptance.

## Linkage and synthetic controls

For every qualifying seed track, combine its topocentric position and angular velocity with 13 logarithmic ranges (0.4–3 au) and 13 radial velocities (−90 to +90 km/s). Add the observing station's heliocentric position and velocity. The ERFA built-in Earth ephemeris and terrestrial rotation supply the observer, using MPC station parallax constants. UT1 is approximated by UTC, TDB by TT, and polar motion is omitted. Compare a solar-gravity Taylor propagation of these Cartesian trial states with all eligible observations on later UTC nights. Require 0.7–3.25 days between consecutive tracklet midpoints and a common trial within 600 arcsec on the second and third nights. The initial angular prefilter and fit criteria are explicit in `link.py`.

Fit all directions in a proposed triple to a six-coordinate solar two-body orbit with numerical DOP853 propagation and light-time correction. Retain fits with ≤1 arcsec angular RMS, ≤2.5 arcsec maximum separation, and nominal heliocentric excess speed >10 km/s. Distinct triples are independent proposals and may share detections. A screen-positive fit would require an independent multi-solution, perturbed original-barycentric refit and identity checks. This short-arc screen alone cannot establish binding energy.

Draw synthetic objects at an actual seed direction and an actual three-night cadence, before applying cuts on the observed rate or linearity, from the same station and interval. The geocentric range is log-uniform over 0.45–2.8 au, the true heliocentric excess speed uniform over 15–60 km/s, and the heliocentric velocity direction isotropic. Use numerical two-body propagation, iterative light time and independent 0.2 arcsec Gaussian tangent-plane noise per coordinate. Force these astrometric detections into the full real background one object at a time. Apply the same angular, timing, residual and fitted-speed cuts as for real triples. A recovery requires all three injected tracklet identities. Count accepted mixed associations separately. In the rectangular subset, reapply its boundaries to every synthetic observation.

Injections assume usable astrometry already exists. Later positions are placed at actual observing times but are not checked against exposed pixels, camera footprints, weather, flux limits or trailing losses. Recovery is conditional on the injection distribution and excludes single-night and two-night objects. Five controls independently rotate the real sky per night by a random RA offset in ±3°, preserving each tracklet's internal geometry while perturbing cross-night alignments. They measure chance-link behavior under this specific control. The calibrated denominator is the supplied-astrometry input class.

## Reproduction

From this directory, install `numpy==2.2.6 scipy==1.15.3 pyerfa==2.0.1.5` under Python 3.12. Run `python extract.py --output run` and compare hashes with the archived receipts. A future live download can change. The original compressed ITF SHA-256 is `7624f8ebbfe37117c353998626b7d2f3782e55bbc02e5e21671b9199072262b1` (9,060,422 lines). The F51 all-sky slice SHA-256 is `5c273cd20fc7feecda82d627bc7feed3a2f01904ebc419899a90c23decaeaef0`. Earlier experiments used already motion-qualified tracklets as cadence templates and produced a 97/100 F51 result. The final control draws from all timing-valid groups before motion cuts; all three slices, input hashes, trial rows and control counts are archived. [ITF receipt](run/input.json) · [F51 result](run/f51-all-result.json) · [F51 trials](run/f51-all-injections.json) · [Rubin result](run/x05-result.json) · [Rubin trials](run/x05-injections.json)

With the archived observation slices in `run/`, execute:

```sh
python link.py --data run --survey X05 --trials 100
python link.py --data run --survey F51 --trials 100
python link.py --data run --survey F51-all --trials 100
python pointing/replay.py --data run --states 20000 --fit-cap 100
python rubin_feed.py run
```
