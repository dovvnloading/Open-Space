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
