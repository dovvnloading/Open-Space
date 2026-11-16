# Systematic HK-first PJ74 OBS sample: executed negative result

Selection was frozen in `SYSTEMATIC_SAMPLE_PLAN.md` at commit `da6a32cb`, before the new photon intervals were inspected. Eight intervals in two independent FITS products were retained, including the previously inspected 20-second control. No threshold was tuned on the new results.

## Actual coverage, not nominal time

The eight intervals span **106 nominal seconds**, but only **three contain any valid nightside photons**. Two contain no photons at all; three more contain only off-disc photons. Thus it would be scientifically wrong to report 106 seconds of nightside survey exposure.

| Product segment / HK indices | Nominal seconds | All photons | Valid night photons | Night 170–200nm | Largest 20ms night bin | Triggers |
|---|---:|---:|---:|---:|---:|---:|
| S02 /145–148 |6|0|0|0|0|0|
| S02 /428–435 |14|16525|7696|1067|10|0|
| S02 /729–739 (prior control) |20|9787|3016|404|5|0|
| S03 /3464–3470 |12|0|0|0|0|0|
| S03 /4257–4262 |10|19903|0|0|0|0|
| S03 /4646–4655 |18|11247|0|0|0|0|
| S03 /5066–5076 |20|11901|242|37|3|0|
| S03 /5364–5367 |6|1576|0|0|0|0|

Total: **70,939 extracted photons;10,954 valid nightside photons;1,508 nightside longwave photons;zero candidate triggers**. This adds 7,938 night photons and1,104 night longwave photons beyond the old control.

The three night-photon-supported spans are5.660,10.287 and0.722 seconds. Their occupied20ms bins number284,513 and36 respectively. These are **photon-support diagnostics, not geometrical or sensitivity-corrected exposure**. Night-photon support starts8.339 seconds into the14-second interval and9.712 seconds into the prior20-second interval. Large portions of telemetry-clear wallclock do not observe the night disc.

## Instrument checks and controls

All selected HK endpoints/runs satisfy the frozen overflow/checksum/HVPS/safety criteria. Actual frame tags51717,51718,51722–51725 map to edgemsk1. Empty intervals have no corresponding active frame tag in the interval. Masks use conservative union over actual tags; ancillary STIM centers use16-pixel half-width exclusions. No night photons fall in these STIM boxes. Exact masked/removed counts, quality bytes and source-row hashes are retained in `systematic/systematic-results.json`.

The original fixed detector requires≥20 longwave photons in20ms and≥60% in one2-degree cell. No night bin reaches even the count cut. Twenty conditional-label diagnostics per interval also yield zero triggers (largest randomized night bin14), but **these are not calibrated false-alarm probabilities**. They preserve local100ms/time, detector-Y segment and night status while breaking the exact detector-X/wavelength relation. Eleven off-disc 20 ms bins reached the 20-photon count threshold, with a maximum of 61 longwave photons. The off-disc control records each bin and its count. Planetary-cell localization applies only to on-disc photons. The earlier analysis clipped off-disc coordinates into cell 0 and incorrectly labeled all eleven bins as spatial triggers.

## Interpretation

This is a reproducible **negative feasibility sample**, not an astronomical discovery and not an impact-rate constraint. The important operational finding is that clean HK selection alone wastes most nominal observing time: the search must select actual nightside footprint coverage before scaling. Eight intervals are not eight independent full-disc exposures. No injection/recovery completeness or response-calibrated physical null has been established. We do not infer an absence of impacts from zero triggers.

## Sources and replay

Products under https://pds-atmospheres.nmsu.edu/PDS/data/jnouvs_3001/DATA/ORBIT-74/ :
- `UVS_S02_805506529_2025192_P74OBS_V01.FIT`
- `UVS_S03_805506529_2025192_P74OBS_V01.FIT`

Run the existing `code/screen_pj74_housekeeping.py`, then `code/systematic_sample.py`, then `code/systematic_analyze.py`. The extractor dynamically validates the photon schema/offsets, bounds timestamp searches, verifies range lengths and Content-Range headers, and retains empty intervals without issuing invalid requests. The source binary cache is ignored by git. Environment reconnection changed network access mid-run; public archive retrieval was resumed after approved access restoration. Final extraction and analysis completed successfully.
