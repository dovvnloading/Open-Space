# Search for inbound hyperbolic small bodies

**No previously unrecognised interstellar object has been validated.** The November 17 JPL snapshot contains 2,290 records with nominal heliocentric eccentricity at least one. Eight have an osculating excess speed of at least 5 km/s. Three are the known interstellar controls. For each of the other five, the eccentricity exceeds one by less than 0.36 of its quoted standard error. This is a census of catalogue solutions meeting those thresholds. Objects with short arcs, missing uncertainties or less extreme nominal speeds can escape the screen.

## Completed

- The [JPL extraction and row-level audit](catalog/) include 527 strictly hyperbolic solutions and 1,763 set exactly to a parabola. The latter supply no positive nominal excess speed. The newest observation represented in this selected JPL snapshot is November 9, 2026.
- In the [bound-orbit comparison](comparison/evidence/screen_summary.json), 91 of the 902 entries with `0.95 ≤ e < 1` have an MPC single-opposition source or a measured arc of at most 30 days. Forty-two of the 76 with usable eccentricity errors reach the parabolic boundary within three formal standard deviations. Two of those 42 have observations in 2025 or later. This comparison includes comets and asteroids. The other 15 of 91 lack a usable error.
- The [Get-CP rolling NEOCP/PCCP orbit file](comparison/evidence/neocp_rolling.orb.gz) contained 275 distinct temporary IDs and zero nominal eccentricities of one or greater. The provider retains current objects plus the last nominal orbit of recent removals. These short-arc nominal solutions have no demonstrated completeness for interstellar visitors.
- The full November 17 isolated-tracklet file contained no Rubin X05 rows from October 19 through November 16 UTC. The latest available 29-day X05 slice ended July 14, 2026 and contained 36 observation lines. [Extract and query receipts](tracklets/)

The speed used for triage is `sqrt(GM_sun (e−1)/q)`. It describes a heliocentric two-body orbit at the catalogue epoch. An incoming speed needs a fit to astrometry and propagation backwards past the planets. The single-element significance omits systematic errors and sensitivity to non-gravitational acceleration.

## Pending at this checkpoint

Independent gravitational and activity-sensitive fits, a bulk refit against the MPC observation file, tracklet linking and recovery tests, and archival follow-up of any surviving lead. No population bound has been derived. C/2025 V1 (Borisov), a previously designated comet distinct from 2I, has a nominal 3.85 km/s excess on a ten-day JPL arc and is under review.

## Prior work and reproduction

The [NEA identity study](../nea-pair-search/README.md) established station holdouts and archival checks. The [Borisov encounter study](../borisov-704021/README.md) supplies orbit constants and an example of testing predicted image positions. The [Rubin binary study](../rubin-jr80-mutual-motion/README.md) documents public Rubin alert cutouts. This search covers different observations and hypotheses. See [the repository's inclusion criteria](../../docs/REMARKABLE_FINDINGS.md).

The result subdirectories contain the scripts, input hashes and machine-readable evidence. Archive inputs change over time. Compare their hashes and dates before comparing separate runs.
