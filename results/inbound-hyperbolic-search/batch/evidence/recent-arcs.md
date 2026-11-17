# Recent astrometric arcs

A direct MPC Observations API refresh matched astrometry for **43 of 44** selected JPL objects. The 44 comprise 31 from the `e ≥ 1` census and 13 recent bound-side comparison objects. Some arcs contain observations through November 16, 2026. [Query, identifier alternatives, HTTP receipts, astrometry hashes and fit outputs](latest_mpc_refits.json.gz)

The second pass returned 28 numerical fits. An earlier pass also returned a solution for C/2022 QE78 on the identical input SHA-256, using 3,967 of 4,225 observations. The joined table retains that result when the retry times out. The first and second outputs are preserved separately. Line counts include observation-continuation records where present. [Initial request receipts and fits](initial_mpc_queries.json.gz)

All six non-control asteroid-class rows in the `e ≥ 1` census yielded nominally bound far-incoming orbits in the live pass. A/2024 U2 provides a clear example. Its 159 MPC lines extend to October 17, 2026; Find_Orb gives reciprocal semimajor axis +0.0000254 au⁻¹ inbound at 244 au. A/2026 T2 and A/2026 T3 also give positive nominal reciprocal semimajor axes at 242 and 246 au. These values describe the particular gravity solution and do not prove every admissible orbit bound. Fits to the complete observations of the three interstellar controls are reported [separately](../../orbits/README.md); the original live packed-ID matching omitted records with only a permanent I designation.

One of the thirteen recent bound-side objects, 2026 UX3, crossed to a nominally hyperbolic inbound orbit in the live fit. It used all nine positions over 7.11 days, including three days beyond the saved JPL last-observation date. The solution at 273 au has `e = 1.00898 ± 0.0185`, only 0.49 formal standard deviations above one. Its 1.80 km/s nominal incoming excess cannot establish its origin. For another comparison object, 2026 RY158, the MPC designation service rejected both its unpacked and packed designation. It remains unresolved here.

See [batch methods](../README.md) and the [full single-opposition or short-arc comparison](../../comparison/evidence/short_arc_screen.csv).
