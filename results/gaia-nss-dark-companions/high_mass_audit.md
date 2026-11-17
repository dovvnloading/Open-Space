# What happened to the 383 stars above the three-solar-mass screen

Twelve had clean RVs on at least two nights with the same independent instrument. Twelve more had clean velocities with insufficient repeat coverage under that rule. The other 359 had no quality-passing RV in the selected DESI, SDSS and LAMOST products. Gaia catalogue summaries have no RV epoch information for the predictive test. [Source-by-source coverage](coverage.csv.gz) · [Coverage totals](coverage.json)

| Gaia DR3 ID | Maximum screen mass (solar masses) | Clean nights | Evidence affecting the compact-object interpretation |
|---|---:|---|---|
| 182029303939094784 | 3.54 | 2 LAMOST | Two nights supply one velocity difference. The star occurs in the Jayasinghe and El-Badry SB1 screens. |
| 275582835635329024 | 3.61 | 7 LAMOST | UCAC4 720-033289 has a reproducible 470-day orbit. Primary-mass models straddle the threshold. The held-out NUV flux and optical Hα wings favor substantial stellar light. |
| 422991301559494528 | 3.09 | 3 LAMOST | Over 52 days the nightly velocities span 0.42 km/s. The nominal 467-day orbit predicts only 0.59 km/s at those phases. |
| 460498563758448512 | 30.45 | 4 LAMOST | The nominal Gaia orbit predicts a 247.35 km/s span. The measured span is 81.30 km/s. None of 4,073 physical Gaia covariance draws achieved a shape-test p-value above 0.01. The star occurs in Jayasinghe et al. |
| 878555832642451968 | 7.82 | 3 LAMOST | The observed 52.80 km/s span disagrees with the nominal 14.84 km/s span and no physical draw passes the shape test. Already screened by Jayasinghe and El-Badry. |
| 1594082579405063168 | 3.08 | 3 APOGEE | The 17-day span covers about 2% of the astrometric period. Observed range 1.67 km/s, and a constant model has χ² 1.83 for 2 degrees of freedom. Present in prior tables. |
| 2224573745414713984 | 3.02 | 3 APOGEE | Thirty days cover 4% of the 751-day period. The measured span is 0.16 km/s. Lowering the flagged giant primary mass by 30% yields a 2.49-solar-mass screening estimate. |
| 2951241529971223808 | 3.51 | 4 DESI | A 17.25 km/s change matches the SB1 prediction over 46 days. The source is a long-period variable, with chromatic variability near half the 276-day period and strongly variable spectroscopic parameters. |
| 3372818683473351936 | 3.18 | 4 LAMOST | The nominal RV model gives χ² 193.1 for 3 degrees of freedom. Fewer than 1% of physical orbit draws pass the shape test. Already screened by Jayasinghe and El-Badry. |
| 3664684869697065984 | 3.24 | 3 APOGEE | The observed change spans 0.90 km/s over 27 days of a 1,220-day orbit. Li et al. explicitly print this source in their 2026 sample. |
| 3998734261904478848 | 3.57 | 2 LAMOST | A two-night baseline supplies one velocity difference. Uncertainty in the 692-day Gaia ephemeris allows compatible draws. Already in two published SB1 lists. |
| 5713565792393224960 | 3.24 | 3 APOGEE | A 6.92 km/s trend matches the predicted 6.43 km/s over 26 days of a 348-day orbit. A 30% reduction of the model giant primary mass yields 2.64 solar masses. Companion light is unconstrained. |

The entries use the larger selected mass where the two input passes differ. The mass column includes edge-on SB1 estimates, astrometric estimates under a zero-companion-light assumption, and mass-function floors when no primary estimate was available. A row in a published selection records prior examination. Several of these studies already identify stellar or spurious solutions. [Initial LAMOST predictions](lamost/shape_checks.csv) · [Extension LAMOST predictions](extension/lamost/shape_checks.csv) · [SDSS and DESI predictions](orbit/survey_union_shapes.csv) · [Prior-table audit](novelty/README.md)
