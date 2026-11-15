# Forbush-event epochs and **independent** geomagnetic storm onsets

**Research status (2026-11-15): completed archival null/label-leakage audit; *not* a new forecasting discovery.**

This is a previously uninvestigated population in Open-Space. The frozen protocol was pushed before measurement in [`protocols/forbush-storm-onset.md`](../../protocols/forbush-storm-onset.md). The full 1963–2016 development/holdout comparison was replayed with a **post-result review correction** to Dst-hour eligibility and matched event minutes in [GitHub Actions #43777026230](https://github.com/dovvnloading/Open-Space/actions/runs/43777026230). The SHA-256-locked inputs, UTC windows, sample counts, shifted-time negative controls, and 3,000 matched null replicates appear in [`summary.json`](summary.json). The initial frozen protocol is preserved, but the correction happened *after* viewing the initial numbers; the holdout is therefore an era split rather than a completely untouched final confirmatory test.

## Quantified result

| FEID epoch population | Eligible 72-hour episodes | Dst < -100 nT onset in next 48h | Year/month-matched mean | Monte Carlo upper-tail p |
|---|---:|---:|---:|---:|
| Development, 1963–2000 | 1,908 | 78 | 65.06 | 0.0460 |
| **Held-out era, 2001–2016** | **768** | **24** | **18.87** | **0.1280** |
| Combined, 1963–2016 | 2,676 | 102 | 84.16 | 0.0220 |

**The pooled association does not pass the held-out test.** Do not claim a cosmic-ray-only precursor.

The important and reproducible anomaly is the **reversal by FEID onset class**:

| Label (first event in episode) | 1963–2016 events | Observed onset matches | Matched null mean | Observed / expected | Matched-tail p |
|---|---:|---:|---:|---:|---:|
| 1: shock + sudden storm commencement | 442 | 68 | 19.68 | 3.46 | <=0.000333 (upper) |
| 9: no shock detected | 2,150 | 30 | 59.92 | 0.50 | <=0.000333 (lower) |
| **2001–2016 only, class 1** | **87** | **18** | **4.78** | **3.77** | **<=0.000333 (upper)** |
| **2001–2016 only, class 9** | **643** | **3** | **12.48** | **0.24** | **0.001666 (lower)** |

Here the event code is **not a measurement of independent galactic cosmic rays**: the catalog explicitly defines class 1 via shock + SSC and class 9 via absence of detected shock. The FEID team describes onset selection as involving SSC/shock, solar-wind changes and sometimes cosmic rays. Such a population is *selected* using parts of the same geospace disturbance that precedes a Dst threshold crossing; post hoc event magnitude and the catalog's own Dst values would add further target leakage. The striking split is therefore evidence that naive prediction studies based on this catalog would be contaminated, not evidence of a new physical effect.

## Independent measurement and temporal guardrails

- Input A: 7,097 FEID annotated event epochs (1957–2016), [public Parquet mirror](https://huggingface.co/datasets/juliensimon/forbush-decreases), SHA-256 `6d88f455a97c5e2baf9a2d6e63f0f6b549865907fe8a2a0a0eff24e3a5a4914f`.
- Input B: 561,024 OMNI hourly index slots (starting 1963), [NASA-derived Parquet mirror](https://huggingface.co/datasets/juliensimon/omni-solar-wind-parameters), SHA-256 `e97d99bd5ea47786d2ab39dd3f65b3f2db193c48ef4c8ba992a6c5862e13662a`. Only event epochs in 1963–2016 enter the result; don't treat calendar slots beyond the observation date as historical measurements.
- Outcome comes **only** from independent hourly `dst_index_nt`, never from FEID `dst_min_nt`, Kp, or post-event decrease amplitude.
- Storm onset is an hourly Dst strictly below -100 nT after **72 completely observed preceding hours with no Dst below -100 nT**. Missing Dst hours veto onset and eligibility. FEID episodes are clustered when successive events are within 72h; the first timestamp/first onset class is retained.
- **No same-hour lookahead:** for FEID event hour H, eligibility requires fully observed, nonstorm Dst hour H−1 and the preceding quiet interval. The average Dst for hour H is *not read* to screen the event; consequently, an event coinciding with a developing storm *within hour H* cannot be ruled out with hourly data. Future onsets are tested after the exact FEID minute in [0,48) hours; bins [0,6), [6,12), [12,24), [24,48) are stored. Subhour physical ordering remains ambiguous.
- The null samples **eligible event-hour slots in the same year and month**, keeping the original FEID event's minute within the sampled hour (seeded 3,000 repetitions). Eligible slots are based exclusively on the prior completed Dst hour, not the uncompleted current hour. Separate negative controls shift event epochs ±7 or ±14 days **within the same month**, retaining only eligible shifted epochs. Exposure and coarse solar-cycle seasonality are matched, but contemporaneous solar shocks, catalog selection, serial dependence and post hoc timing remain.
- The development comparison 1957–2000 was shortened to **1963–2000** because the available independent hourly Dst series starts in 1963. The held-out years are 2001–2016, as frozen.

## Review-stage shifted-time falsification control

After discovering the within-hour lookahead in the original implementation, the candidate-event epoch test was recomputed with the preceding completed Dst hour and the original FEID minute. Eligible episodes and matched counts remained 2,676/102 overall and 768/24 for 2001–2016, but year/month-matched null means rose modestly (82.28 → **84.16** overall; 18.56 → **18.87** in the later era), reducing the apparent pooled excess.

The registered ±7 and ±14-day negative controls were also added at review stage. In the **2001–2016** era, class-1 FEID episodes yield **18/87 = 20.7%** matches at the cataloged epochs, versus **2/55 = 3.6%** at +7 days, **5/60 = 8.3%** at −7 days, **2/46 = 4.3%** at +14 days, and **0/45** at −14 days. These shifted groups exclude month-crossing and ineligible exposure slots and therefore have **different denominators**. They do not constitute a randomization p-value or establish an independent cosmic-ray precursor, but show that the original timing alignment weakens when moved away from the shock/SSC-annotated epoch. Full controls by class and era are in [`summary.json`](summary.json).

**Research integrity:** This time-bin correction and added controls were implemented after the first results were observed, so the later era remains a transparent robustness split rather than a pristine, prospectively untouched confirmatory test. No new forecasting or astrophysical discovery is asserted.

## A second, independent data-provenance failure

A proposed raw-neutron-monitor check would normally be valuable, but the separate [`juliensimon/neutron-monitor-cosmic-rays`](https://huggingface.co/datasets/juliensimon/neutron-monitor-cosmic-rays) Parquet archive fails a full-file station-identity audit. **All 497,835 rows** use one of **five HTML DOCTYPE fragments** in both `station` and `station_name`. Each bad identifier appears 99,567 times; 20,863 rows lack count rates. The exact five identifiers, SHA-256, and source URL are preserved in [`neutron-monitor-identity-audit.json`](neutron-monitor-identity-audit.json), reproduced by [GitHub Actions #43769771145](https://github.com/dovvnloading/Open-Space/actions/runs/43769771145) and [`examples/audit_neutron_monitor_identity.py`](../../examples/audit_neutron_monitor_identity.py).

This makes station-level attribution impossible **in that public archive snapshot**. The numeric measurements may be legitimate but misattributed; the audit does not establish fabrication. Do not use its station names for cross-station validation. Obtain verified, pressure-corrected hourly counts directly from the [Oulu station](https://cosmicrays.oulu.fi/) or [NMDB](https://www.nmdb.eu/), with provider acknowledgement, before claiming independent cosmic-ray skill.

## Replay

```sh
python -m pip install pyarrow
python -m unittest discover -s tests -p test_forbush_dst_independent.py -v
python -m unittest discover -s tests -p test_audit_neutron_monitor_identity.py -v
python examples/forbush_dst_independent.py --fetch --draws 3000 --output /tmp/forbush-result.json
python examples/audit_neutron_monitor_identity.py --fetch --output /tmp/neutron-identity.json
```

The workflow attaches both machine-readable reports and **fails CI** if a source SHA-256 or any reported measurement differs from the committed packet. A source update needs explicit review and a new reference measurement. No raw archives are committed to Open-Space.

## Prior art and interpretation

- FEID team account of analyst-derived onset times and mixed upstream measurements: [Belov et al., IZMIRAN abstract (2025)](https://seminar.pgia.ru/abstracts_2025.pdf).
- [Forbush-effect onset may use SSC or solar-wind indicators, MNRAS 2023](https://academic.oup.com/mnras/article/521/3/4544/7081343).
- [Interacting-ICME Forbush timing and moderate FD/Dst links, MNRAS 2026](https://academic.oup.com/mnras/article/552/1/stag1664/8785133).
- NASA primary series: [OMNI hourly data](https://omniweb.gsfc.nasa.gov/html/ow_data.html).

**Conclusions:** real, large-sample evidence against interpreting this post-annotated event catalog as a clean forecasting instrument. Neither a first-in-literature precursor nor operational forecasting skill is claimed.
