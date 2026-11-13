# First completed observing visit

The coordinator processed every primary BOSS main-program galaxy meeting the published quality and redshift cuts on SDSS plate 6497, MJD 56329. The selection returned 12 distinct sky positions across four nonempty pages, followed by an empty page. A separate query for up to 100 fibers on the visit returned the same 12 IDs. Temporary SDSS HTTP 503 errors interrupted archive access during the run. Reinvoking the coordinator resumed the saved pages and completed the visit.

| Measurement | Result |
| --- | ---: |
| Galaxies screened | 12 |
| Narrow features found | 333 |
| Wavelength matches to the bundled source table | 92 |
| Wavelength matches to the three-line airglow screen | 0 |
| Outside both reference tables | 241 |
| Same-sign groups in at least two distinct targets, rest frame | 47 |
| Same-sign groups in at least two distinct targets, observer frame | 36 |

The recurrence calculation used a maximum group span of 150 km/s, and each target contributed at most one detection per group. Different observers looking at the same public coadds can reproduce these measurements. All targets share a plate and observing visit, the source table emphasizes emission lines, and the airglow screen covers only three lines. The comparison between frames has no calibrated null distribution. Each of the 47 groups is a screening result. One explored group is consistent with known Ca II K absorption. [Ranked groups and member feature IDs](first-visit-groups.csv) · [Measurements, source digests and selection query](first-visit.json) · [Calcium investigation](../docs/FIRST_INVESTIGATION.md)

## Reproduce the complete visit

From the checkout used for the recorded software fingerprint, install the archive extra and run the campaign in a new directory. The first visit used deterministic screening and zero provider calls. The model-assisted investigation of the initial three spectra was performed separately by research agents. The worker and model setup are documented for further campaigns. [Campaign runbook](../docs/CAMPAIGNS.md) · [Agent setup](../docs/MODELS.md)

```sh
python -m pip install -e '.[archive]'
python -m open_space.campaign --root runs/replicate-first-visit \
  --plate 6497 --mjd 56329 --limit 3 --pages 5 --max-calls 0
```

Use the same command after an archive outage. Verify that `state.json` records `exhausted: true`. Compare each report's source and CSV digests with the published record before comparing the counts. The completed campaign writes the full membership of each group, including its originating observation and feature IDs, to `recurrence.json`. [Method and limitations](../docs/RECURRENCE.md)


## Follow-up experiments

- The fixed-wavelength test applies eight original predictions and a Ca II K control to three new visits. The JSON records every hit, covered miss and unavailable measurement under each mask. [Readout and reproduction](../docs/HOLDOUT_RESULTS.md) · [Per-observation evidence](holdout.json)
- The fitted-spectrum check subtracts the SDSS model and measures which of the original detections survive at the same threshold. [Measurements](template-residuals.json)
- The exact and randomized redshift reassignments compare the observed degree of grouping with chance alignments under specified redshift-exchange assumptions. [Analysis](../docs/NULL_EXPERIMENT.md) · [Distributions](redshift-null.json)
- The repeat case compares independently acquired spectra of one original target with two alternative continuum and exposure-mask checks. [Analysis](../docs/REPEAT_OBSERVATION.md) · [Evidence](repeat-observation.json)
- A fresh whole-spectrum search starts from observed flux minus the survey fit and ranks surviving groups, retaining the observer-frame and cross-visit controls. [Analysis and residual candidate register](../docs/RESIDUAL_SEARCH.md)
