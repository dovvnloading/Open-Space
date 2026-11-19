# 4XMM-DR14 QPE control-recovery audit — 2026-11-19

**Verdict: POSITIVE CONTROL CALIBRATION ONLY. No validated new QPE, no independent galaxy-nucleus archive screen, no calibrated null or occurrence-rate limit. Keep PR #205 draft.**

## Recovery denominator and results

The denominator is **six explicitly QPE-active XMM-Newton EPIC-pn observations**, with **29 bright published eruption peaks**: four controls and counts from Webbe & Young (2023) Table 1 ([DOI](https://doi.org/10.1093/rasti/rzad015)), plus eRO-QPE3/4 from Arcodia et al. (2024; [DOI](https://doi.org/10.1051/0004-6361/202348881)). The RX J1301.9+2747 2020 observation additionally has **10** model-fitted bright/weak/partial QPEs in Giustini et al. (2024; [DOI](https://doi.org/10.1051/0004-6361/202450861), Table B.1); the eight bright-segment benchmark is kept separate.

| QPE-active source | XMM obsID | Published bright | Threshold peaks | Measured consecutive peak intervals (hours) | Published benchmark |
| --- | --- | ---: | ---: | --- | --- |
| GSN 069 | 0851180401 | 5 | 5 | 8.333, 9.333, 8.583, 9.333 | ~8–10 h; 5 bright eruptions |
| RX J1301.9+2747 | 0864560101 | 8 | 8 | 1.333, 5.000, 2.917, 5.417, 3.250, 5.583, 3.917 | irregular 1–6 h; **8/10** fitted Giustini peak timestamps matched to ±1200 s; two fitted weak/partial events not recovered by pn |
| eRO-QPE1 | 0861910201 | 2 | 2 | 11.583 | two broad overlapping peaks; second is **truncated at the end of exposure**; XMM within-epoch separation is not the 18.5 h NICER mean recurrence |
| eRO-QPE2 | 0872390101 | 9 | 9 | 2.417, 2.417, 2.500, 2.333, 2.500, 2.333, 2.500, 2.333 | nine published flares at roughly 2.4 h |
| eRO-QPE3 | 0883770101 | 2 | 2 | 20.417 | ~20.4 h |
| eRO-QPE4 | 0883770401 | 3 | 3 | 9.750, 14.750 | **9.82 ±0.03** and **14.70 ±0.02** h; differences expected with 300-second sampling |

**Operational recurrence-source recovery:** **6/6 = 100%** of the six named positive-control XMM epochs show at least two detected peaks and an interval in each published broad range. **Bright-eruption peak-count overlap:** **29/29 = 100%** *capped count diagnostic* (not one-to-one event-level completeness). **Independently verified eruption-level detection efficiency: NOT ESTABLISHED.** In particular the 29/29 number cannot be promoted to a measured sensitivity, 2020 RX has only 8/10 table-fitted timestamps matched, one eRO-QPE1 flare is partly observed, and independent ground-truth timing is not tabulated/matched for every published flare. This is in-sample calibration rather than holdout validation.

## Data provenance and correction to failed PR #205 run

[Original baseline run #44470633598](https://github.com/dovvnloading/Open-Space/actions/runs/44470633598) did not recover controls: an overly restrictive 120 MiB download ceiling, incorrect native 1–13-second bin peak finding, NaN JSON and unreliable 4XMM-to-PPS source indexes. It screened just twelve requested time-series, with nine parsed, and never validated a new event.

**Current retrieval:** ESA XMM Science Archive PPS `level=PPS&extension=FTZ&instname=PN&name=SRCTSR`; each target has a 3 GiB stream ceiling, SHA-256 receipt, bounded 14-minute wall budget, with at most three targets per push-triggered <55-minute workflow. The source product is the per-detection merged `PNX000SRCTSR8###.FTZ` observation-source light curve, independently selected by FITS `SRC_RA/SRC_DEC` within 15 arcsec of the published X-ray position. QPE3 notably uses **P0883770101PNX000SRCTSR8007.FTZ**, **2.28 arcsec** from the published source; its historical DR14 source number 5 did not map to the present PPS source identifier. Position checks and raw tar/product SHA-256 hashes are in the per-batch receipts and outputs, including the separately run eRO-QPE1 astrometric recheck.

**Light curves:** energy-resolved `RATE1+RATE2+RATE3`, the 0.2–2 keV background-subtracted `HDUCLAS2=NET` EPIC-pn data; rate errors combined in quadrature; `BACK1V+BACK2V+BACK3V` local scaled background; exposure-weighted **300 s** bins with >=65% good exposure and no filling GTI gaps. Thresholds in [`../control_thresholds.json`](../control_thresholds.json): >=4 typical-error prominence, >=3.5-sigma 300s excess, >=7 integrated-sigma in fixed 450, 850, 1300 s Gaussian kernels, >=3x baseline, >=60min separation, scaled-background flare veto. The end-of-observation masked padding in replay restores the original good/bad-bin distinction without fabricating light-curve measurements; eRO-QPE1's truncated second flare must not be counted as a complete eruption.

**Actions receipts/results:** [batch 1](controls_batch_01.json) ([run](https://github.com/dovvnloading/Open-Space/actions/runs/44482844158)), [batch 2](controls_batch_02.json) ([run](https://github.com/dovvnloading/Open-Space/actions/runs/44483930261)), [batch 3](controls_batch_03.json) ([run](https://github.com/dovvnloading/Open-Space/actions/runs/44485930501)), [eRO-QPE1 header audit](controls_batch_04.json) ([run](https://github.com/dovvnloading/Open-Space/actions/runs/44488897010)); [no-network, all-six replay](control_replay.json) with [machine-readable gate](control_gate_replay.json). A one-line report-key error in run #44488076115 was repaired, and the corrected no-network replay run #44488549801 passed. These data products, per-source `control_curve_*.csv` and SHA-256 receipts allow peak and recurrence checks. Original pre-control outputs remain historical and must NOT be interpreted as calibrated negative or discovery results.

## Independent-referee acceptance gate (NOT PASSED)

1. **Individual published event matching is incomplete:** table-based one-to-one comparisons were performed for the ten fitted RX J1301 peaks (8/10, with special treatment needed for the MOS-only/partial ones); comparable independent onset/peak lists and uncertainty budgets are not yet aligned for all other epochs, particularly truncated/overlapping eRO-QPE1. Demonstrating in-sample count agreement is not enough.
2. **False-alarm and completeness model missing:** only bright-source and synthetic background-flare checks exist; real blank-field, stellar, AGN and injection/recovery negative controls, proper multiple-testing allowance across 33,488 detections and an out-of-sample validation set are required.
3. **Instrumental vetoes incomplete:** PPS NET and local BACK1-3 reject some particle events but no independently checked `FBKTSR`/high-energy background GTI, aperture contamination, EPIC-pn pile-up, or spectral-temperature/hardness-lag/hysteresis analysis exists. These would be mandatory for any *candidate*.
4. **Parent galaxy-nucleus sample is not defined:** 33,488 soft/long 4XMM-DR14 detections are X-ray detections, **not** a confirmed or optically cross-matched census of galaxy nuclei. Cross-match galaxy centers, host redshifts and Gaia/stellar counterparts, register all surveyed observations and nondetections, freeze the cut and injection-inferred sensitivity, then screen under the validated rule. An archive-wide null or QPE occurrence rate cannot be inferred from the earlier 12 requested source time series.
5. **No new candidate has been scientifically vetted:** no post-gate host screen occurred, hence zero independently checked new candidate records, and no TNS/ATel/ADS/prior-QPE/TDE/AGN cross-epoch novelty checks can support a claim. Existing known QPEs are positive controls, never novel sources.

**AT2019qiz and Ansky** are documented *coverage exclusions*, not missed positive controls: their later QPE detections were not a qualifying QPE-active observation in the frozen DR14 XMM epoch. Neither should artificially dilute or inflate the six-epoch denominator.

**Disposition:** retain draft PR #205; **do not merge**, do not add `[REMARKABLE]` or `docs/REMARKABLE_FINDINGS.md` entry, and do not label this a calibrated negative survey. Gate remains `archive_screen_permitted=false` until above checks genuinely pass.
