# Working as a research agent

Start from a checkout and identify the question, observing unit and falsification test. Read the first investigation to see how an apparent catalog miss can be explained. Select the appropriate operating mode for the available evidence. [Worked example](docs/FIRST_INVESTIGATION.md) · [Science protocol](docs/SCIENCE.md)

- Use a unique campaign root for each instrument selection, code revision and model setup. The campaign coordinator resumes a frozen target page and records its cursor. Give separate plate visits separate directories. [Campaign commands](docs/CAMPAIGNS.md)
- For individual experiments, enqueue a verified observation manifest and claim jobs with the worker. A job stores a copy of the validated rows, so a later input-file edit cannot change a queued experiment. Export evidence with `python -m open_space report JOB_ID`. Place `--db PATH` before the subcommand when using a different ledger.
- A feature identifier is local to one observation and code revision. Always pair it with the observation and run IDs. Use the coverage counts when interpreting an empty feature list. Preserve both the catalog match and observer-frame alternatives when they coincide.
- Before escalating an unexplained feature, check masks, continuum dependence, the observer-frame wavelength, common source transitions, plausible blends, redshift uncertainty and the history of reports on this target. Record completed and pending checks separately. Independent models evaluating the same pixels share the same observational evidence.
- Run calibration and synthetic tests separately from observational screening. Reproduce a measurement from its primary source or describe exactly which input bytes you used. Keep the provenance for negative and inconclusive results.
- Return reviewable contributions through a pull request. Provide the measured result and its limits, including any failures to reproduce it. Changes to a live campaign budget or model endpoint should be explicit operator choices.

## Optical and infrared surveys

Use the discovery coordinator to assign optical, infrared and full-frame observations. Freeze the input products and comparison population. Use an appropriate counterfactual test or a held-out epoch, verify the location of a flux deficit, and preserve competing source associations and catalog query receipts. Agent rankings are proposals for another measurement. Inspect completed child evidence before changing a scientific disposition. The released target sequence, field coordinates and next unexamined rank are in the operator guide. [Operating guide](docs/DISCOVERY.md) · [Scientific evidence](docs/DISCOVERY_RESULTS.md)

## Reuse the independent-field validation

Start with the published eight-prediction holdout and the separate calcium control if a question depends on the first campaign. Reproduce the entire list, including uncovered wavelengths. Preserve the discovery run when deriving residual spectra. Check whether a feature survives subtracting the survey model and whether its observed wavelength also recurs on the same plate. An in-sample fit can absorb genuine structure. Independent measurements still matter. [Holdout and fit subtraction](docs/HOLDOUT_RESULTS.md) · [Blind residual search](docs/RESIDUAL_SEARCH.md) · [Same-target repeat](docs/REPEAT_OBSERVATION.md)
