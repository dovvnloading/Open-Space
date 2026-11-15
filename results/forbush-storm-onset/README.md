# Forbush-event epochs and **independent** geomagnetic storm onsets

**Research status (2026-11-15): completed archival null/label-leakage audit; *not* a new forecasting discovery.**

This is a previously uninvestigated population in Open-Space. The frozen protocol was pushed before measurement in [`protocols/forbush-storm-onset.md`](../../protocols/forbush-storm-onset.md). The full 1963–2016 development/holdout comparison was actually run in [GitHub Actions #43769425860](https://github.com/dovvnloading/Open-Space/actions/runs/43769425860); real data, input SHA-256 hashes, exact UTC windows, sample counts, and 3,000 matched null replicates appear in [`summary.json`](summary.json).

## Quantified result

| FEID epoch population | Eligible 72-hour episodes | Dst < -100 nT onset in next 48h | Year/month-matched mean | One-sided permutation p |
|---|---:|---:|---:|---:|
| Development, 1963–2000 | 1,908 | 78 | 63.42 | 0.0337 |
| **Held-out era, 2001–2016** | **768** | **24** | **18.56** | **0.1203** |
| Combined, 1963–2016 | 2,676 | 102 | 82.28 | 0.0127 |

**The pooled association does not pass the held-out test.** Do not claim a cosmic-ray-only precursor.

The important and reproducible anomaly is the **reversal by FEID onset class**:

| Label (first event in episode) | 1963–2016 events | Observed onset matches | Matched null mean | Observed / expected | Matched-tail p |
|---|---:|---:|---:|---:|---:|
| 1: shock + sudden storm commencement | 442 | 68 | 19.26 | 3.53 | <=0.000333 (upper) |
| 9: no shock detected | 2,150 | 30 | 58.67 | 0.51 | <=0.000333 (lower) |
| **2001–2016 only, class 1** | **87** | **18** | **4.69** | **3.84** | **<=0.000333 (upper)** |
| **2001–2016 only, class 9** | **643** | **3** | **12.14** | **0.25** | **0.001999 (lower)** |

Here the event code is **not a measurement of independent galactic cosmic rays**: the catalog explicitly defines class 1 via shock + SSC and class 9 via absence of detected shock. The FEID team describes onset selection as involving SSC/shock, solar-wind changes and sometimes cosmic rays. Such a population is *selected* using parts of the same geospace disturbance that precedes a Dst threshold crossing; post hoc event magnitude and the catalog's own Dst values would add further target leakage. The striking split is therefore evidence that naive prediction studies based on this catalog would be contaminated, not evidence of a new physical effect.

## Independent measurement and temporal guardrails

- Input A: 7,097 FEID annotated event epochs (1957–2016), [public Parquet mirror](https://huggingface.co/datasets/juliensimon/forbush-decreases), SHA-256 `6d88f455a97c5e2baf9a2d6e63f0f6b549865907fe8a2a0a0eff24e3a5a4914f`.
- Input B: 561,024 OMNI hourly index slots (starting 1963), [NASA-derived Parquet mirror](https://huggingface.co/datasets/juliensimon/omni-solar-wind-parameters), SHA-256 `e97d99bd5ea47786d2ab39dd3f65b3f2db193c48ef4c8ba992a6c5862e13662a`. Only event epochs in 1963–2016 enter the result; don't treat calendar slots beyond the observation date as historical measurements.
- Outcome comes **only** from independent hourly `dst_index_nt`, never from FEID `dst_min_nt`, Kp, or post-event decrease amplitude.
- Storm onset is an hourly Dst strictly below -100 nT after **72 completely observed preceding hours with no Dst below -100 nT**. Missing Dst hours veto onset and eligibility. FEID episodes are clustered when successive events are within 72h; the first timestamp/first onset class is retained.
- Candidate events during a storm, within 72h of a severe Dst reading, or with missing coverage are excluded. Outcomes require an onset *after* the event time and in [0,48) hours; bins [0,6), [6,12), [12,24), [24,48) are stored. Clock resolution near zero remains ambiguous.
- The null samples **eligible UTC hours in the same year and month** as each FEID episode (seeded 3,000 repetitions), not arbitrary hours during active storms. Exposure and solar-cycle seasonality are matched but other confounders, contemporaneous solar shocks, catalog selection, serial dependence and post hoc timing remain.
- The development comparison 1957–2000 was shortened to **1963–2000** because the available independent hourly Dst series starts in 1963. The held-out years are 2001–2016, as frozen.

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

The workflow attaches both machine-readable reports. A source revision may change counts; compare SHA-256 hashes first. No raw archives are committed to Open-Space.

## Prior art and interpretation

- FEID team account of analyst-derived onset times and mixed upstream measurements: [Belov et al., IZMIRAN abstract (2025)](https://seminar.pgia.ru/abstracts_2025.pdf).
- [Forbush-effect onset may use SSC or solar-wind indicators, MNRAS 2023](https://academic.oup.com/mnras/article/521/3/4544/7081343).
- [Interacting-ICME Forbush timing and moderate FD/Dst links, MNRAS 2026](https://academic.oup.com/mnras/article/552/1/stag1664/8785133).
- NASA primary series: [OMNI hourly data](https://omniweb.gsfc.nasa.gov/html/ow_data.html).

**Conclusions:** real, large-sample evidence against interpreting this post-annotated event catalog as a clean forecasting instrument. Neither a first-in-literature precursor nor operational forecasting skill is claimed.
