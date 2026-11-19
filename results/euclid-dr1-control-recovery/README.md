# Euclid DR1-Foundation: control-first high-z reselection

**Claim (pre-reselection): the prior 462,500-row strict VIS-dropout catalogue is incomplete for known high-z quasars; the 478-object morphology cohort is not a defensible discovery list.** The prior selection recovered **3/24 (12.5%)** of catalogue-matched z>6.5 controls, and its borrowed NISP shape band retained **1/24**. This follow-up will diagnose every matched control from its exact source Parquet bytes before freezing revised cuts. There is **no new candidate, calibrated full-footprint completeness or spectroscopic discovery** at this stage.

The original full HATS census, source/ETag receipts and controls are frozen in [PR #208](../euclid-dr1-wide-highz-vetting/README.md). The 96 unmatched literature controls cannot be counted as failures of a pixel-coverage-complete survey. The z>7 luminosity-function expectation of approximately 8–19 for a hypothetical complete 1,900 deg² is not a candidate yield.

## Precommitted analysis sequence

1. Re-read only the 24 matched quasars and 14 matched brown dwarfs from the **exact ETag-verified original S3 HATS partitions**; audit original Boolean conjuncts, flags, signed/absolute VIS S/N, NISP colours, ellipse and MUMAX_MINUS_MAG object by object. Output: [raw control audit](output/control_audit.json), [diagnostic summary](output/audit_summary.json).
2. Freeze an explicitly specified high-recall *eligibility* rule and a separate high-priority/contaminant rejection rule. Measure numerator, denominator, Wilson uncertainty and dwarf contamination on matched controls; do not train and assess a threshold on the same labels without caveat.
3. Reselect from existing frozen strict-row batch registers where applicable and stream additional HATS source rows whenever changed cuts admit previously excluded sources. Checkpoint in bounded batches through **one** path-filtered, push-triggered Actions job. Never turn 462,500 strict rows into a claim of completeness for a changed cut.
4. Investigate top surviving objects with actual Euclid VIS/Y/J/H mosaic pixels using IRSA SIA metadata and **public S3 FITS byte-range reads** (Astropy WCS/Cutout2D); review RMS/FLAG where accessible, edges, blends and artefacts; record unWISE/CatWISE matching and motion, public spectrum/literature receipts. An unreturned service query is not a nondetection.
5. Independent referee audit; publish a ranked photometric list with measured conditional recovery and open validation gates. **Only a newly confirmed spectrum supports a discovery claim.**

All computation runs on GitHub-hosted runners, not the network-disabled chat sandbox. Generated \`output/**\` and narrative edits do not retrigger the workflow. Receipt snapshots and limitations take priority over premature candidate counts.

This work has made use of the Data Release (DR1) data from the Euclid mission of the European Space Agency (ESA).
