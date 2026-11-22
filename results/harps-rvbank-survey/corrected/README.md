# HARPS radial-velocity companion survey

**The search fitted 458 periodic signals in 375 of 983 catalog targets.** The targets that met the cadence cuts contributed 41,856 usable spectra. A fitted peak can arise from noise, activity, binarity, sampling or a known companion. Individual novelty claims require review of the source papers and object identities.

The input was J/A+A/683/A125/table4 (verified corrected subset). See [input checksums](inputs.csv), [exposure dispositions](epoch-ledger.csv.gz) and [target eligibility](stars.csv). The [RV fit table](signals.csv) records periods, semiamplitudes, false-alarm estimates and diagnostic checks. The [external vetting table](vetting.csv) records companion-mass constraints and object-level checks, including unresolved evidence. The [acceleration table](accelerations.csv) records long-term trends.

20 fits passed all internal RV tests. The [RV survivors](rv-survivors.csv) include their catalog and literature dispositions. The [associated bibliography](survivor-literature.csv) contains indexed references. Another 11 fits met the statistical thresholds and retained unresolved checks in both the RV and external ledgers. See [inconclusive periods](inconclusive-periods.csv). Of 52 planned external queries, 47 returned usable responses. [Coverage by source](source-coverage.csv) · [Query receipts](catalog-inputs.csv)

The initial period search recovered 3 of 7 prespecified published control periods above the adopted FAP and model-improvement thresholds. The denominator includes every prespecified control, including misses and targets outside coverage. Each planet counts once across its [independently fitted archive series](published-control-series.csv). See [control recovery](published-controls.csv).

The [injection grid](completeness-grid.csv) summarizes 227,616 injected orbits in observed residuals. Its denominator varies with orbital period and available baseline. It measures recovery of the blind periodogram peak before external vetting. The [individual injection cells](injections.csv) retain both recoveries and misses.

## Signal dispositions

- inconclusive: 14
- prior signal: 62
- rejected: 382

## Reproduce the run

Use the source at commit `fc5a5d90b1e34ff3190fd2396acaf61a8ef45dca`, the [survey workflow](../../../.github/workflows/harps-rvbank-survey.yml) and [installed software versions](packages.txt). The workflow records its commands, including shard count and selection criteria. Compare the exact archive and catalog input checksums between runs. [Actions run](https://github.com/dovvnloading/Open-Space/actions/runs/44988368934)
