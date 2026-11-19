# Preliminary extreme-halo scale — numerical bug-repair receipt

**Status: provisional computed geometry/HMF baseline; not a catalogue test. No measured null or discovery.**

A GitHub Actions run ([#44445915491](https://github.com/dovvnloading/Open-Space/actions/runs/44445915491), source commit `93e32dc575b2c99eb1d61efb2d5b1b2a8671aeb2`, artifact ID `14495125270`) computed redshift-integrated Colossus ST99/Tinker08 halo abundances on a nominal 0.54 deg² COSMOS-Web area. The original `halo_bounds.json` has SHA-256 `ae10b31f2f96c4f0d49f37a75e6df53c594973b7ac48f65fba4e317041b0f7ed`. A **coding error** applied base-10 exponentiation to a natural-log mass interpolation and printed impossible 10^28–10^29-solar-mass maxima. Fixed in commit [4204c26](https://github.com/dovvnloading/Open-Space/commit/4204c265832d437958262ba1472a30c82fdc15df); separate referee Poisson-tail checks added in [73f91b6](https://github.com/dovvnloading/Open-Space/commit/73f91b611295277b88a0f4df72825ae4d26f2935).

The table is an **exact algebraic inverse of the exponent bug**, not an independently verified post-fix rerun: `M_correct = exp(log10(M_archived))`, since the archived `M_archived = 10**(ln(M_correct))`. The corrected HMF code and referee await execution on unmodified new outputs. All values below are **log10(maximally converted stellar baryon mass in solar masses)**, using `0.16 M_halo`, not absolute cutoffs.

| z bin | ST99 FOF expected-one | ST99 FOF 95% max-halo CDF | Tinker08 M200m expected-one | Tinker08 M200m 95% max-halo CDF |
| --- | ---: | ---: | ---: | ---: |
| (6,7] | 11.548 | 11.870 | 11.394 | 11.714 |
| (7,8] | 11.261 | 11.582 | 11.094 | 11.414 |
| (8,9] | 10.994 | 11.315 | 10.815 | 11.134 |
| (9,10] | 10.743 | 11.063 | 10.553 | 10.872 |
| (10,12] | 10.521 | 10.834 | 10.318 | 10.632 |
| (12,15] | 10.086 | 10.398 | 9.862 | 10.176 |

**Interpretive cautions:** The 95% CDF is for the random most-massive halo, **not** a 5% ΛCDM significance threshold for an observed galaxy. Counts are Poisson in an unmasked nominal area; exclude completeness, selection, cosmic variance, mass-systematic spread and lensing. ST99 FOF and Tinker08 M200m are different halo mass definitions; high-z HMF extrapolation is uncertain. **No observed COSMOS2025 galaxy has been compared to these figures**, because the CDS archive download timed out. No own joint posterior, independent SED refit or spectroscopic/interloper vetting exists.

Scientific disposition: **NO INFERENCE about a COSMOS2025 violation or absence.** Methodological baseline only. See [referee review](INDEPENDENT_REFEREE.md).
