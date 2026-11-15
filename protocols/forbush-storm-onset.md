# Forbush decreases versus independent geomagnetic-storm onset

Status: **preregistered experiment, no scientific result claimed at this checkpoint**. Created 2026-11-15.

## Unoccupied question

The existing `results/` investigations include individual solar-wind events and solar-system/stellar work, but **no population-level test of neutron-monitor Forbush decreases against independently measured geomagnetic-storm onsets**. This is a new dataset and question, not a follow-up to a previous Open-Space target.

Test whether epochs in the IZMIRAN FEID catalog have a genuine, time-ordered association with the *subsequent* onset of intense geomagnetic activity in an independent hourly Dst series, beyond season- and solar-cycle-matched random epochs. A positive result would motivate direct, real-time neutron-monitor testing. It would **not** by itself establish operational forecasting.

## Frozen guardrails before seeing results

1. **Different data streams:** Use FEID event epochs and the *independent* Kyoto/OMNI hourly Dst time series; do not use `dst_min_nt`, `kp_max`, `ap_max`, minimum Bz, or other geomagnetic labels included in the FEID catalog as outcomes. These catalog fields are not independent validation.
2. **Timeliness:** A completed-event Forbush amplitude is retrospective and must never be described as a real-time predictor. Use only event epoch for the primary lag test. Also investigate whether the FEID epoch denotes shock arrival, detection start, or an analyst-assigned time; this limits advance-warning claims.
3. **Timing and labels:** Define storm onset as the first crossing of hourly Dst below -100 nT, with no crossing in the preceding 72 hours. Examine subsequent windows [0,6), [6,12), [12,24), and [24,48) hours separately. Reject events already in an active storm. Record time conventions in UTC.
4. **Independence:** Group Forbush catalog epochs into episodes at least 72 hours apart; show results with and without this exclusion. Separate discovery years 1957–2000 and held-out years 2001–2016; never tune a threshold after opening the holdout.
5. **Controls:** Draw exposure-matched random epochs preserving year/month and Dst eligibility, and repeat using at least one shifted-time negative control. Resample event episodes rather than correlated hourly measurements. Report actual counts, confidence bounds and sampling coverage.
6. **Falsification:** A retrospective correlation, a null, an ambiguous epoch definition or a target-label leak is **not a new astrophysical discovery**. Publish the result, including non-support, without claiming a previously unknown mechanism.

## Planned reproducible artifacts

- `results/forbush-storm-onset/README.md`: measured result and competing interpretations
- `results/forbush-storm-onset/summary.json`: exact counts, controls, inputs and configuration
- `examples/forbush_storm_onset.py`: independent analysis and data-retrieval recipe
- `tests/test_forbush_storm_onset.py`: artificial timing fixtures and leakage checks

## Sources to verify

- Primary FEID/IZMIRAN catalog: http://spaceweather.izmiran.ru/eng/dbs.html
- Public FEID redistribution (7,097 records, 1957–2016): https://huggingface.co/datasets/juliensimon/forbush-decreases
- Hourly Dst: Kyoto World Data Center, available via NASA OMNI2 https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2_all_years.dat
- NASA OMNI data caveats and fill values: https://omniweb.gsfc.nasa.gov/html/ow_data.html

Do not use this checkpoint as evidence that the hypothesis is supported.
