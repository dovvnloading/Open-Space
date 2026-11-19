# GALAH DR4 PISN-fingerprint search — evidence ledger

**Claim (as of opening): No validated new pair-instability-supernova (PISN) descendant has been identified.** This is an *open* search, not a null result. An absent or unusable measurement must never be interpreted as a low abundance.

## Question and precommitted tests

Test whether high-quality metal-poor GALAH DR4 stars exhibit the predicted strong odd–even pattern, suppressed cobalt and zinc, and no neutron-capture enhancement, and whether that chemistry prefers a PISN-enrichment model over core-collapse and Type Ia explanations. Look up APOGEE DR17 and LAMOST measurements independently; use LAMOST J1010+2358 only as a published, disputed control, not as a discovery.

1. **Catalogue receipt:** pin exact DR4 input/table/release, UTC retrieval date, query text, returned row counts, data SHA-256 and column metadata. Keep all query and filtering code.
2. **Element-level QC:** reject failed spectral fits and flagged abundances individually; never equate missing Na/Al/Co/Zn/Sr/Ba (or upper limits) with measured depletions; record effective counts and limits per element.
3. **Frozen screening:** predeclare metallicity, uncertainty, odd-even and Zn/Co conditions *before* inspecting candidate ranks; keep all tested stars and selection denominators.
4. **Model comparison:** use documented, versioned yield tables for Heger & Woosley (2002) and Takahashi et al. (2018), compare against independently specified core-collapse and Ia mixtures using identical likelihoods, observational covariance, limits and model-selection penalties. Report a sensitivity/false-positive calibration, not a naked best-fit score.
5. **Controls:** recover published patterns of the contested J1010+2358 with correct survey cross-identification; perform APOGEE DR17/LAMOST consistency and cross-survey zero-point checks.
6. **Vetting:** inspect blend/neighbor contamination, multi-epoch RV scatter, line-fit residuals and primary literature (SIMBAD, SAGA), using source-specific receipts.
7. **Independent referee:** audit selection, fitting, withheld controls, completeness and potential look-elsewhere effects. A validated discovery alone qualifies for the repository [REMARKABLE] policy.

## Status

- Completed: repository policy read; baseline scientific question and conservative evidence bar frozen here.
- Pending: all astronomical catalogue retrieval, full-sample scan, model fitting, cross-survey checks, control recovery, vetting, referee audit.
- Evidence level: **not yet determined; no [REMARKABLE] claim**.

## Files and reproduction

Queries, scripts, response receipts, counts and hashes will be stored alongside this README. Report an access failure as a blocked search rather than as a physical null.
