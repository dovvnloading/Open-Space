# Does measured outgassing supply 3I/ATLAS's nongravitational acceleration?

**Initial claim (before replay): No unexplained acceleration is established.** This is a new comparison of full-arc astrometry and measured gas production, *not* a new interstellar-object discovery. Quantitative claims remain pending the independent GitHub Actions orbit fit and physical recoil grid.

## Frozen design before accessing new observations

- Astrometry: obtain MPC 80-column 3I observations and retain their raw bytes, endpoint dates, selection counts and SHA-256. Fit gravity plus the RTN Marsden nongravitational law using a documented N-body integrator, with JPL DE440/DE440s planetary positions, barycentric ICRF coordinates, finite observer position and iterative light time. Freeze quality cuts **before** inspection: supported ground-based optical stations only, remove identical (station, epoch, RA, Dec) duplicates; one arcsecond per coordinate floor and station/night total inverse-variance cap of four; three fit passes (soft-L1, then two linear), 4-sigma two-coordinate clipping after each of the first two passes. Compare independently to JPL SBDB, separating formal errors and astrometric model uncertainty. Preserve gravity-only, CO2-like r^-2 RTN, and water Marsden-law sensitivity tests and pre/post-perihelion time splits.
- Gas: independently compile time-stamped measured **production rates** from primary papers, distinguishing nuclear gas from distributed coma grains, and upper limits from detections; do not turn CN/Ni directly into total gas mass.
- Recoil: `|a| <= eta * sum_i(Q_i * m_i * v_i) / ((4/3) pi rho R^3)`, with `0 <= eta <= 1`. Report matched temporal and RTN projections; flux and speed assumptions are **bounds**, not measured jet axes. An acceleration-magnitude agreement alone does not identify a jet direction.
- Falsification: an excess would require *lower-bound* fitted acceleration beyond the *upper-bound* physically plausible recoil at overlapping dates, stable under selection/force-law and size/density uncertainty. No special label without independent validation. Validate held-out astrometry and output provenance before any substantive conclusion.

Earlier work that this result **reuses, not repeats**: [3I gravity and original-energy positive control](../inbound-hyperbolic-search/orbits/README.md) and [3I–383P radial geometry](../3i-atlas-383p-radial-crossing/README.md). The independent fit and momentum-budget question are new. [Remarkable findings criteria](../../docs/REMARKABLE_FINDINGS.md).

## Reproduction

The dedicated workflow on branch `research/3i-atlas-outgassing-recoil` is push-triggered only by changes to this results folder's *source* files, explicitly excluding `generated/**`. Its single job performs all network retrieval and computation on a GitHub-hosted runner and commits output files in one batch. Full receipts, numerical outputs and referee assessment will be added here after the run.
