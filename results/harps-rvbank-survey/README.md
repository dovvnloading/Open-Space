# HARPS radial-velocity companion survey

**No HARPS observations have been searched in this branch yet.** There are no measured periods, companion masses, false-alarm probabilities, injection recoveries or archive checksums to report. The exposure count, eligible-target count and candidate yield are unknown. Zero completed fits cannot constrain planet occurrence.

## Scope

The intended population is the public HARPS-RVBank release, including stars with too little time coverage for an orbital fit. The inventory must record those exclusions. An ESO Phase 3 fallback must identify its own footprint, velocity reduction and diagnostics. The two reductions of the same exposure cannot be counted as independent observations.

Each proposed companion needs a fitted period and RV semiamplitude, a mass estimate with its stellar-mass source, a search-level false-alarm probability, and explicit dispositions for stellar activity, rotation, sampling aliases, fibre configuration, binarity and published counterparts. An RV acceleration and an astrometric anomaly can motivate an orbital fit. Either measurement alone leaves orbital geometry and mass degenerate. Unresolved cases belong in the signal ledger.

## Execution

Archive acquisition and analysis belong to the single manual GitHub Actions workflow `harps-rvbank-survey.yml`. GitHub requires its workflow path to be registered on the default branch before it can receive a manual dispatch. The research branch contains the workflow for review. See [GitHub's dispatch documentation](https://docs.github.com/en/actions/reference/workflows-and-actions/events-that-trigger-workflows#workflow_dispatch).

After registration, run the following from an authenticated checkout with Actions write access.

```sh
gh workflow run harps-rvbank-survey.yml --ref survey/harps-rvbank-companions
```

Use the workflow run's source commit, downloaded-input receipts and dependency lock to identify the search. Retain failed source attempts, excluded targets, rejected signals and unresolved checks when interpreting any completed run.

## Prespecified search and controls

A target enters the blind period search with at least 12 observing nights covering at least 30 days. The search spans 1.5 days to the individual target's baseline. Every eligible target receives a Keplerian fit. The pipeline records up to three sequential periodic components, tests for a residual slope and measures injection recovery over the available baseline. The completeness estimate refers to the periodogram detection stage. [Spectroscopic sources](sources.md) · [RV methods](rv-method.md) · [External checks and astrometry](vet-method.md)

Seven published orbital periods across HD 4308, HD 40307 and HD 69830 serve as positive controls. They enter the blind scan with identical cuts and weights to every other star. The final report lists each recovery, miss and coverage exclusion. A period match requires agreement within both one independent frequency spacing and five percent in frequency. The comparison uses the original reported periods, rounded to the precision given in the control table. [Udry et al.](https://doi.org/10.1051/0004-6361:20054084) · [Mayor et al.](https://arxiv.org/abs/0806.4587) · [Lovis et al.](https://arxiv.org/abs/astro-ph/0703024)

HD 32564 and the HARPS drift of HD 25673 already have entries in Open-Space. Their signals are prior results for this survey. [HD 32564](../hd-32564/README.md) · [HD 25673](../coralie-rv-acceleration/README.md)

## Run products

The run saves separate exposure and target ledgers, all RV fits, activity and alias tests, the injection grid, catalog queries, astrometric fits and matching bibliography records. It tests that the combined shards contain exactly one outcome for every target, preserves recoveries and misses for the controls, and writes a checksum manifest. The read-only workflow uploads these as the `harps-evidence` artifact for inspection before any catalog claim. That artifact includes a proposed root catalog edit with the measured recovery of circular 30-day, 10 m/s injections and the corresponding line 9 count change. The patch is generated only after the result exists. The source responses are saved in the same run as separate snapshot artifacts.

Runtime commands and dependency pins are in the [workflow](../../.github/workflows/harps-rvbank-survey.yml) and [requirements](requirements.txt). For an identified run, retrieve its evidence with the following command, substituting the numerical run ID.

```sh
gh run download RUN_ID -n harps-evidence --dir harps-evidence
```
