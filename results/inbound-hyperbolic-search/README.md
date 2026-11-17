# Search for inbound hyperbolic small bodies

**No previously unrecognised interstellar object has been validated.** Searching 2,130 eligible Pan-STARRS 1 tracklets on 17 nights produced no three-night orbit passing the high-speed selection. Testing at locations close to actual detections recovered 35 of 37 simulated objects for which three nights of nearby observations existed. Only 37 of 20,000 simulated paths met that sky-and-time condition. These numbers measure the search on a particular sample and cadence.

## Completed measurements

| Test | Result | Scope |
| --- | --- | --- |
| [Full-precision JPL census](catalog/README.md) | 2,290 unique records with heliocentric `e ≥ 1`, including 527 strictly above one. Eight have nominal excess speed ≥5 km/s. Only 1I, 2I and 3I also exceed five quoted eccentricity errors. | For each of the five other fast nominal solutions, `(e−1)/σₑ < 0.36`. Only 537 of all 2,290 had a positive quoted eccentricity error. |
| [Published original barycentric comparison](catalog/code-original.csv) | 301 catalogue matches. Their preferred original reciprocal semimajor axes were positive in 286 cases and negative in 15. | Two negative preferred solutions exceeded five formal errors. Both have published alternative force models or arc selections compatible with binding. See the [object-by-object discussion](catalog/README.md). |
| [Independent astrometric fits](orbits/README.md) | All three interstellar controls recovered at 250 au, with incoming speeds of 26.34, 32.28 and 57.97 km/s. C/1980 E1 and C/2024 L5 were bound before their Jupiter and Saturn encounters, respectively. | Distinguishes incoming energy from the heliocentric osculating energy used to select the catalogue. |
| [Bound-side comparison](comparison/evidence/screen_summary.json) | Of 902 JPL orbits with `0.95 ≤ e < 1`, 91 were from an MPC single-opposition source or had arcs ≤30 days. Of those, 76 had usable eccentricity errors and 42 reached the parabolic boundary within three formal errors. | Thirteen of the 91 had observations in 2025 or later, including two of those 42. |
| [NEOCP/PCCP rolling orbit file](comparison/evidence/neocp_rolling.orb.gz) | 275 distinct current or recently removed temporary IDs, zero nominal `e ≥ 1`. | File published November 15, 2026. Nominal solutions on short arcs do not exclude unbound alternatives. |

### C/2025 V1 (Borisov)

This previously identified active comet is distinct from 2I/Borisov. Its catalogue orbit has a nominal osculating excess speed of 3.85 km/s. The current MPC record contains 117 observations on ten days. A gravity-only fit gave an incoming `1/a = −0.01598 ± 0.00316 au⁻¹` at 250 au (5.05 formal standard deviations). Allowing radial and transverse water-ice outgassing gave `−0.02682 ± 0.02618 au⁻¹` (1.02). Dropping all seven ZTF observations reduced the gravity-only result to 1.71. The known ZTF measurements include reported r magnitudes of 17.49 and 17.70 on October 29, 2025, and 15.41 on November 1. [Fits and holdouts](orbits/README.md) · [ZTF observation matches](comparison/photometry/evidence/mpc_ztf_matches.json)

In a direct fit constrained to zero incoming energy at 50 au, the cost was Δχ² = 25.37 under gravity. It fell to 1.06 with freely signed radial and transverse Marsden accelerations and 3.39 when enforcing outward radial acceleration. These are differences in fit quality. They do not quantify the probability of an interstellar origin. [Independent model definitions and evidence](orbits/)

### What the tracklet search measures

The 29-day window ran from October 19 through November 16, 2026 UTC. The full MPC Isolated Tracklet File contained 7,170 F51 observations in that window. Eight real triples survived the broad angular prefilter, and zero had a common trial orbit passing the next gate. The final orbit-fit selection required a nominal heliocentric excess speed above 10 km/s. Five shuffled-sky controls produced zero accepted links. [Tracklet selection, dynamics and receipts](tracklets/README.md)

| Injection experiment | Result | Interpretation |
| --- | --- | --- |
| Supply three detections groups at real observing times | 99/100 recovered, exact 95% binomial interval 94.55–99.97% | Measures linkage given the three tracklets. The later visits were chosen without a sky-coverage requirement. |
| Require nearby real astrometry at each exposure on three nights | 37/20,000 synthetic paths qualified within 1.5°. Of those, 35/37 were recovered, 94.6% (exact 95% interval 81.8–99.3%). Fourteen of 15 inbound examples were recovered. | The proxy uses detections near each predicted direction. Thirty-five of the 37 supported paths relied on November 1–3. Within 0.5° only one path qualified, and within 0.1° none did. |

The simulated geocentric distances span 0.45–2.8 au and the injected heliocentric excess speeds span 15–60 km/s. Real-pointing proximity does not establish exposure on a working detector, flux sensitivity, trailing loss or survey completeness. No bound on interstellar space density follows. The latest 29-day ITF window had zero Rubin X05 rows. The last populated X05 window ended July 14, 2026, with only two nights of eligible observations. The [Rubin submission-feed audit](tracklets/) separately checks observation dates against ingestion dates.

## Checks in progress at this checkpoint

The bulk independent Find_Orb fit to the published MPC observation file is in progress. Its row-level output will distinguish far-inbound solutions from insufficient observations, failed fits and model-dependent solutions. Direct MPC queries supplement the dated bulk file with recent arcs. Archived observations of recently unconfirmed NEOCP objects are also being checked.

A public ZTF footprint search for C/2025 V1 identified ten exposures before the known October 29 arc. A difference-image stack searches the nominal ephemeris and its uncertainty region; it yielded no validated precovery. The search did not apply image mask products, and one frame had a problematic-quadrant quality flag. Its empirical ranks and approximate aperture depths are exploratory, not validated sensitivity limits. [Footprints, epochs and archive receipts](catalog/image-search/)

## Prior work and reproduction

The [Borisov encounter study](../borisov-704021/README.md) supplied the starting propagation conventions, the [NEA identity study](../nea-pair-search/README.md) established station holdouts, and the [Rubin binary study](../rubin-jr80-mutual-motion/README.md) documented the status of public Rubin alert cutouts. The results here use different observations or test different hypotheses. The subdirectories provide the scripts, frozen inputs or their hashes, query receipts and controls.

This research uses data and services from the IAU Minor Planet Center, NASA/JPL, CODE, Pan-STARRS, ZTF, the Asteroid Institute and other archives credited in each analysis. [Repository inclusion criteria](../../docs/REMARKABLE_FINDINGS.md)

Frozen input bytes can be checked from the repository root with `python results/inbound-hyperbolic-search/verify_inputs.py`. The [historical 38eef5a workflow](orbits/historical-38eef5a/README.md) preserves the earlier execution separately; its fit and the newer models use the same astrometry and are not independent observations.
