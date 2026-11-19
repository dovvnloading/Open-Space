# Does measured outgassing supply 3I/ATLAS's nongravitational acceleration?

## Result from the independent 19 November 2026 GitHub Actions replay

**Conclusion — no validated excess beyond natural outgassing, but also no demonstrated three-dimensional gas-jet fit.** The public MPC orbit shows statistically detectable nongravitational acceleration. The measured CO2/CO/H2O budget supplies its *magnitude* for some low-mass nuclei; for the HST pV=0.04 radius of 1.3±0.2 km, specific observed gas epochs cannot supply it under fiducial recoil parameters. This is a model-dependent mass/activity tension, **not** a validated requirement for another physical force.

**Astrometry:** The replay downloaded **8,497** raw MPC 80-column positions from **2025 May 8 to 2026 April 14**, of which **8,345** pass the frozen ground-optical and station rules. The fit retained **8,301** and obtained symmetric inverse-square radial/transverse/normal coefficients, in units of 10^-8 au/day²:

| Component | Independent coefficient ± formal 1σ | JPL SBDB #55 coefficient |
| --- | ---: | ---: |
| Radial A1 | +1.856 ± 0.315 | +4.618 |
| Transverse A2 | +3.865 ± 0.240 | −0.956 |
| Normal A3 | −0.745 ± 0.023 | −0.829 |

The gravity-only fit has **0.753 arcsec** per-coordinate RMS, improving to **0.637 arcsec** with inverse-square RTN recoil. A standard water-law 3-component fit has A ≈ (+2.157,+3.820,−1.157) in the same units, showing dependence on the assumed distance law. The time-ordered 10%-of-data holdout has **0.875 arcsec** RMS versus **0.620 arcsec** on its training subset.

**Not like-for-like with JPL:** Retrieved JPL SBDB orbit **#55** (solved 2026 Oct 21) is a *temporary JWST-support* orbit using **117 observations**, **2025 Oct 3–2026 Oct 16**, with fitted outgassing delay **DT=+28.465 days**, and the DE441 force ephemeris. The present fit uses 8,345 filtered MPC ground observations to **2026 Apr 14**, a zero-delay model, and DE440s. JPL and this fit have a **74.7°** acceleration-vector separation, driven in particular by differing A2; this is **not an extra-physics residual** when the data selections and acceleration timing differ. The *preperihelion-only* 9-parameter optimization failed to converge, and the *postperihelion-only* A2 is very poorly constrained (formal σ≈50 in the stated units). These failed sensitivity checks are retained in [orbit.json](generated/orbit.json), not suppressed.

**Measured gas versus required momentum.** At fixed density **500 kg/m³** and gas exhaust speed **500 m/s**, the following are required momentum efficiencies `η = nucleus mass × |fitted acceleration| / (gas mass production × exhaust speed)`. `η > 1` cannot be supplied by the *listed* gas for that particular size/speed/density scenario:

| Measured epoch; independent alternative gas extractions kept separate | Radius 0.42 km | Radius 1.30 km |
| --- | ---: | ---: |
| JWST 2025 Aug 6, *inner/nucleus-centered* CO2+CO+H2O | 0.305 | **9.04** |
| JWST 2025 Aug 6, *terminal/coma-integrated* | 0.166 | **4.92** |
| JUICE/MAJIS 2025 Nov 2 H2O+CO2 (CO not measured by this instrument) | 0.0487 | **1.44** |
| JWST 2025 Dec 22–23 CO+CO2+H2O | 0.138 | **4.10** |

The small 0.42 km model that easily passes the scalar recoil bound may be difficult to reconcile with the HST photometric nucleus at assumed low albedo; a real reconciliation requires the albedo, density, jet momentum geometry and source-coupled gas flux, not merely a successful small-radius row. The full grid spans **0.22–2.8 km**, **300–1000 kg/m³** and **200–1000 m/s**. CO2/H2O and CO/H2O ratios, OH-water proxy assumptions, distributed ice-grain sources, upper limits, and large-aperture water can strongly bias momentum assigned to the nucleus. The measurements do **not** determine radial/transverse/normal jet orientations or a continuous Q(t). No simultaneous full 3-D, same-arc physical gas-recoil fit exists here; a robust residual requiring extra physics is not demonstrated.

The complete public inventory with scholarly links and original-versus-revised measurements is in [SOURCES.md](SOURCES.md) and [production.csv](production.csv). See [computed orbit and coefficient comparison](generated/REPORT.md), [raw MPC positions](generated/mpc_3i.obs.gz), [recoil grid](generated/recoil_grid.json), [JPL JSON snapshot](generated/sbdb.json), [fetch receipts](generated/retrieval_receipts.json), [pipeline hashes](generated/pipeline.json), and the [independent referee assessment](REFEREE.md). These are archived in one path-filtered push-only GitHub Actions job. **Classification: routine, conditional measured replication; no [REMARKABLE] label or log entry.**

**Pre-registered claim (written before any new replay): No unexplained acceleration is established.** This is a new comparison of full-arc astrometry and measured gas production, *not* a new interstellar-object discovery. Quantitative claims remain pending the independent GitHub Actions orbit fit and physical recoil grid.

## Frozen design before accessing new observations

- Astrometry: obtain MPC 80-column 3I observations and retain their raw bytes, endpoint dates, selection counts and SHA-256. Fit gravity plus the RTN Marsden nongravitational law using a documented N-body integrator, with JPL DE440/DE440s planetary positions, barycentric ICRF coordinates, finite observer position and iterative light time. Freeze quality cuts **before** inspection: supported ground-based optical stations only, remove identical (station, epoch, RA, Dec) duplicates; one arcsecond per coordinate floor and station/night total inverse-variance cap of four; three fit passes (soft-L1, then two linear), 4-sigma two-coordinate clipping after each of the first two passes. Compare independently to JPL SBDB, separating formal errors and astrometric model uncertainty. Preserve gravity-only, CO2-like r^-2 RTN, and water Marsden-law sensitivity tests and pre/post-perihelion time splits.
- Gas: independently compile time-stamped measured **production rates** from primary papers, distinguishing nuclear gas from distributed coma grains, and upper limits from detections; do not turn CN/Ni directly into total gas mass.
- Recoil: `|a| <= eta * sum_i(Q_i * m_i * v_i) / ((4/3) pi rho R^3)`, with `0 <= eta <= 1`. Report matched temporal and RTN projections; flux and speed assumptions are **bounds**, not measured jet axes. An acceleration-magnitude agreement alone does not identify a jet direction.
- Falsification: an excess would require *lower-bound* fitted acceleration beyond the *upper-bound* physically plausible recoil at overlapping dates, stable under selection/force-law and size/density uncertainty. No special label without independent validation. Validate held-out astrometry and output provenance before any substantive conclusion.

Earlier work that this result **reuses, not repeats**: [3I gravity and original-energy positive control](../inbound-hyperbolic-search/orbits/README.md) and [3I–383P radial geometry](../3i-atlas-383p-radial-crossing/README.md). The independent fit and momentum-budget question here are new. [Remarkable findings criteria](../../docs/REMARKABLE_FINDINGS.md).

## Reproduction

The dedicated workflow on branch `research/3i-atlas-outgassing-recoil` is push-triggered only by changes to this results folder's *source* files, explicitly excluding `generated/**`. Its single job performs all network retrieval and computation on a GitHub-hosted runner and commits output files in one batch. Full receipts, numerical output and referee assessment have been archived after the replay. Output-only commits are excluded from retriggering and document-only changes are also excluded.
