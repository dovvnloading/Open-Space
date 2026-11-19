# PISN/CCSN/Type-Ia comparison: specified but not executed

**Status: NOT EXECUTED.** No retrieved, hashed, compatible Heger–Woosley/Takahashi/CCSN/Ia yield arrays have been fitted to these stars; no Bayes factor, likelihood ratio, PISN probability or preference is claimed in this result. This is the explicitly pending inference layer, not evidence for either progenitor class.

## Immutable model inputs required before evaluation

| Family | Upstream provenance | Requirements |
| --- | --- | --- |
| PISN, nonrotating | Heger & Woosley (2002), *ApJ* 567:532–543; [StarFit model database](https://2sn.org/starfit/data/) | Document grid version and model masses, explosion energy/metallicity assumptions, full isotope yields and SHA-256 |
| PISN, rotating | Takahashi et al. (2018), *ApJ* 857:111; [VizieR J/ApJ/857/111/yields](https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=J/ApJ/857/111/yields) | Pin all 1602 published yield rows, rotation, initial mass, model IDs and source checksum |
| Core collapse | Explicit CCSN grid (e.g. Heger & Woosley 2010, or Nomoto et al. 2013) | Comparable ejected-mass normalization and mixing/fallback nuisance; avoid cherry-picking |
| Type Ia / mixed | Published Iwamoto et al. (1999) W7/WDD models (or comparable modern yield grids) plus CCSN mixing | Fit Fe-contributing mixture fraction; same penalties and priors as PISN/CCSN |

No published tables are bundled or adapted in the present deposit. A generic toy abundance vector is **not** a substitute for a physical yield grid.

## Pre-analysis fitting specification

1. Convert all ejected isotope masses into elemental number abundances (sum isotopes; divide by atomic weights); compute \`[X/Fe] = log10[(N_X/N_Fe)_yield] - log10[(N_X/N_Fe)_solar]\` with **one documented, frozen solar composition** across families. Treat ejecta mixing and dilution consistently; [Fe/H] may constrain dilution but do not freely adjust each [X/Fe] coordinate.
2. For each measured [X/Fe] \`y_i\` use \`log N(y_i; mu_i(theta), sigma_i² + tau_i²)\`, with documented observational covariance when accessible and shared nuisance model-systematic scatter \`tau\`. For cross-element ratios [Na/Mg], retain their shared Mg/Fe covariance rather than treating ratios and original abundances as independent data points.
3. For any upper limit \`U_i\` use a one-sided censored log likelihood \`log Φ((U_i - mu_i)/s_i)\`. No data for absent/unreliable elements; no inference that a missing Sr/Al/Co detection is a low abundance. If [X/Fe] supplied as absolute abundance, convert consistently before fitting.
4. Optimize or integrate within each PISN and CCSN/Ia model **family** under precommitted prior grids and nuisance parameters. Report \`2(log L_PISN - log L_alt)\` and a comparable penalty/evidence for grid multiplicity and extra mixture freedoms, including look-elsewhere calibration. Specify the favored *non-PISN* family if found.
5. Test J1010's original table as a censored-data control and *independently* attempt updated high-resolution abundance measurements from Skúladóttir et al.; report how the model preference changes and explicitly compare CCSN and CCSN+Ia.
6. Freeze a mixed-population injection experiment on the observed errors, flag completeness and [Fe/H], including a distribution of PISN mixture fractions and controls under CCSN/Ia, **before** claiming a search-wide selection function or physical null.

## Why skipped here

The primary hard-threshold screen yielded **zero** eligible five-element survivors, and the important extremely metal-poor ([Fe/H] <= −2.5) regime has **zero** such complete-five stars. Most of the most diagnostic abundance measurements are undetectable or missing for the majority of low-metallicity stars (e.g. Sr 17/5189 and Al 217/5189). Any uncalibrated model ratio on that nonrandom subset would not support a population claim. The code deposit confines itself to a quantitative observational **selection-level null**, not an unsubstantiated PISN-model comparison.

## Related cross-survey and optical validations

- APOGEE DR17: compare only genuinely common detected abundances (e.g. Mg, Na, Al, Co); identify by Gaia DR3 source_id or carefully checked coordinates, enforce ASPCAPFLAG and STARFLAG, and document DR17-specific zero points; **not run**.
- LAMOST: match by Gaia source_id/sky position with epochs, separate low- and high-resolution pipeline patterns, do not mix published disputed J1010 Xing Table 1 with updated values; **not run**.
- RV scatter, spectral blends, GALAH line residuals, SIMBAD and SAGA cross-identification: necessary for any future survivors, but none passed the detected-element screen; **not run**.

A complete physical null requires reanalysis of censored objects and confirmed modeling/selection function beyond this first frozen five-detection cut.
