# Known-planet timing controls

We measured 35 transits across 15 JWST/NIRSpec PRISM visits from programs 6456, 9256 and 12492. The measurements use the integration barycentric time column and test the sign, magnitude and planet dependence of the 2024 JWST-based forecast errors. Planets b and e in the December 2025 DDT exposure arrive **5 seconds late** and **198 seconds early**, respectively. The eight-minute mask padding in the additional-transit search covers both displacements.

The paired observations show a clear contrast across the same visits. From June–July 2025 planet e arrives 116–131 seconds early in six visits. The five visits from 2025-10-30 through 2025-12-11 give 160–198 seconds early. On 2026-06-06, it arrives 277 seconds early, whereas b arrives 16 seconds late. During the 2025-11-23 visit d arrives 379 seconds early and h 550 seconds late, with b offset by less than one second. The proposal team had already modeled orbital timing anomalies and discussed an outer-planet hypothesis. Our measurements supply a reproducible check on a fixed, older forecast. They do not identify the cause of the deviations.

| Visit | Planet | Measured minus forecast | Measurement scatter, 3- and 9-minute blocks | Published forecast uncertainty |
| --- | --- | ---: | ---: | ---: |
| 2025-10-30, 9256-007 | d | −367 s | 39, 23 s | 124 s |
| 2025-11-23, 9256-012 | d | −379 s | 3, 1 s | 133 s |
| 2025-11-23, 9256-012 | h | +550 s | 5, 3 s | 231 s |
| 2025-12-11, 12492-001 | b | +5 s | 3, 1 s | 30 s |
| 2025-12-11, 12492-001 | e | −198 s | 3, 1 s | 87 s |
| 2026-06-06, 9256-013 | e | −277 s | 4, 1 s | 110 s |

All 35 rows are in [times.csv](times.csv). Observations after the forecast's 2024 release test it prospectively with respect to that fit. We selected and inspected these data retrospectively. The forecast errors have shared orbital covariances; combining the rows as independent significance tests would overstate their evidence for any dynamical hypothesis. Individual offsets are within about three times their published forecast errors.

## Measuring the times

We compared normalized 20-second median flux bins to a symmetric transit template, with independent midpoint, depth, duration and ingress time, robust soft-L1 loss, and a local quadratic baseline. The shape uses smooth cubic ingress/egress and a flat bottom. An alternate linear baseline is fitted separately. Other forecast planet transits are masked with three minutes of padding. Symmetric-shape approximations, flares, starspots and reductions can produce systematic differences larger than the statistical errors.

The published Allen et al. seven midtimes from three visits furnish an external timing check. Our first two b times agree to within 2 seconds and our other quiet-visit midtimes agree to within approximately 8 seconds. The strongly flaring 2024-07-05 e transit disagrees by about 77 seconds. The team explicitly treated this flare in its reduction. Our fit's formal errors underestimate systematics for that event, which we exclude from the trend quoted above. Values labeled block errors use half the 16–84 percentile spread of 60 deterministic wild-bootstrap fits per elapsed-time block size. Each block receives one random residual sign. These estimates describe local fitting noise, not a guarantee of absolute accuracy.

Agol et al.'s measured-transit reference is [arXiv:2409.11620](https://arxiv.org/abs/2409.11620). The forecast table, `T1_forecast_JWST_chi_19Aug2024.csv`, is in [the author's source archive](https://arxiv.org/src/2409.11620). Its SHA-256 is `55c3534cbaaf75f428ef7ce05341c70c62c76d33a67a673e0c532549d656d345`. [forecast_subset.csv](forecast_subset.csv) preserves the visited rows' planet, epoch, forecast midpoint (column 3) and inflated uncertainty (column 6). The published five July 2024 post-fit timing checks include the visits used here. Allen et al. published the first three 6456 visits' seven midtimes. [Allen et al., Table 2 and discussion](https://doi.org/10.3847/1538-3881/ae28cb)

## Replay

Run `OPENBLAS_NUM_THREADS=1 python replay.py /path/to/cache` from this directory with NumPy and SciPy. The script retrieves or verifies all 15 exact MAST products listed in [products.csv](products.csv). It uses `BJD_TDB_NRS1` plus 2,400,000.5 for barycentric JD. Existing files with a different checksum stop the replay, which protects the comparison against later pipeline revisions.
