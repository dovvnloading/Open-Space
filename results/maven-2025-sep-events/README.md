# Extreme 2025 August Mars energetic-particle event: full MAVEN final-season audit

**Result class: independently reproducible major Solar-System event measurement and pre-event Mars-specific forecast concordance; physical CME attribution not yet established.** 
**Observation window:** 2025-08-01 through 2025-12-05 UTC, with gaps. 
**New analysis:** 2026-11-15 (UTC 2026-11-16). 
**Location:** Mars; MAVEN orbiter, not Mars surface. 

This investigation is a **new** scientific lane independent of the existing Mars crater, impact-date, ion-tail, comet, and SPHEREx investigations. It adds new quantitative results without changing prior results or code.

## Main result

The newly indexed/repackaged MAVEN Key Parameters reveal an **extreme, sustained, dual-species high-energy particle enhancement at Mars from 2025-08-22 21:00 through 2025-08-29 18:00 UTC** (earliest and latest selected clock hours, with gaps). In the fully searched data:

| Evidence | Measured number |
|---|---:|
| SEP source rows checked | **1,189,114** |
| Distinct SEP times with SPICE ephemeris match within 9 s | **1,189,114** |
| Observed UTC dates / calendar span | **99 / 127** (2025-08-01 through 2025-12-05) |
| Hours with ≥100 observations | **2,376** |
| Predeclared dual-species candidate hours | **176** |
| Main event's selected hours | **125** |
| Main event's high-significance dual-species samples | **56,032** |
| Of these, samples at >3,000 km spacecraft altitude | **38,385** |
| Distinct MAVEN orbits containing joint event detections | **46** |
| Strongest ion and electron hourly medians, at **2025-08-24 05:00 UTC** | **157,000 and 616,500 (cm² s sr)⁻¹** |
| Survey reference of hourly median ion/electron flux, ≥100-sample hours | **30.55 and 25.2125 (cm² s sr)⁻¹** |
| Peak / survey reference | **~5,140× ions, ~24,450× electrons** |

The peak is an **hourly median**, not an isolated telemetry sample. During the August 24 05:00–05:59 UTC clock hour, **450/450** examined SEP rows passed the simultaneous 3-of-4 directional-view flux and signal-to-statistical-uncertainty criteria; **205** were above 3,000 km. August 24 02:00–10:00 includes numerous consecutive extreme-flux clock hours. This is not a single-pixel or single-instrument glitch. Directional views *share telescopes* and are **not four independent instruments**.

### Independent, pre-event external prediction: a substantial physics lead

NASA's **Moon-to-Mars (M2M) DONKI** forecast record for the **2025-08-21 08:23 UTC far-side halo CME** shows two sets of WSA–ENLIL predictions for **Mars** submitted **2025-08-21 13:15 UTC**, before the event:

| Independently forecast Mars arrival (UTC) | Model characterization | Offset of MAVEN's 2025-08-24 05:00 peak |
|---|---|---|
| 2025-08-24 **00:21** | 2 AU leading edge | +4 h 39 m |
| 2025-08-24 **04:20** | 5.5 AU leading edge | **+40 min** |
| 2025-08-24 **09:38** | 2 AU shock front | −4 h 38 m |
| 2025-08-24 **09:40** | 5.5 AU shock front | −4 h 40 m |

**The independently recorded MAVEN flux peak occurs within the entire independently pre-specified range of these Mars-arrival forecasts** (00:21–09:40 UTC). The Aug 21 forecast was archived before either the peak or this reanalysis, making the date association more informative than a retrospective search alone. The NASA DONKI record also links subsequent energetic-particle events, but this alone does not establish a cause.

**Crucial competing association:** NASA DONKI separately records a **2025-08-20 05:12 UTC CME**, with an independently issued forecast **2025-08-20 11:51 UTC** that put its Mars arrival at **2025-08-23 22:00 UTC** (7 hours before the measured Aug 24 flux maximum). The competing event prevents confidently attributing the MAVEN enhancement uniquely to the Aug 21 CME without shock signatures, particle spectra, and proper heliospheric propagation modeling.

- [NASA M2M DONKI Aug 21 CME #40685 and its original Mars arrivals](https://kauai.ccmc.gsfc.nasa.gov/DONKI/view/CME/40685/1)
- [NASA M2M DONKI Aug 20 CME #40665 and its original Mars arrivals](https://kauai.ccmc.gsfc.nasa.gov/DONKI/view/CME/40665/1)
- [NASA CCMC DONKI access and contemporaneous records](https://ccmc.gsfc.nasa.gov/tools/DONKI/)

The event **does not** establish new-onset solar physics, an exact shock-arrival time, an unreported CME, any radiation dose to rovers, or a causal explanation for MAVEN's end-of-2025 communications loss. It **does** identify and quantify an extreme, multi-orbit Mars radiation-environment episode in newly available 2025 data, and links it to a genuinely prior, Mars-directed solar eruption forecast. The precise source and relative timing of energetic particle acceleration and shock impact remain testable questions.

## Other full-survey candidate events and negative controls

The same *frozen* detection threshold, without per-date tuning, found:

| Period, selected clock hours UTC | Selected hours | Joint SEP rows | High-altitude joint rows | Distinct orbits | Disposition |
|---|---:|---:|---:|---:|---|
| Aug 22 21:00 – Aug 29 18:00 | 125 | 56,032 | 38,385 | **46** | **Strongest episode; pre-event Mars-forecast temporal corroboration** |
| Oct 30 07:00 – Nov 01 03:00 | 34 | 14,394 | 10,214 | **13** | Independent additional event candidate; source unassociated |
| Nov 03 14:00 – Nov 04 11:00 | 15 | 6,459 | 4,741 | **7** | Independent additional event candidate; source unassociated |
| Oct 29 23:00 – Oct 30 00:00 | 2 | 759 | 709 | **1** | **Reject as independent multi-orbit confirmation; track as possible precursor/single-orbit** |

All values are reproducibly copied from [the derived machine-readable evidence](summary.json) and [full day-wise counts](daily.csv). The *third-party mirror* includes no rows for **2025-11-05 through 2025-12-02**; the full 127-day calendar period has only 99 observed days, not an unbroken dataset. Do **not** infer that no radiation events occurred during its 28 absent days.

Independent environmental context for the long August episode, from *different MAVEN sensors* (matched by UTC hour):
- MAG event-window median hourly field magnitude ~**15.16 nT**, versus survey-hour median ~**9.33 nT**.
- SWIA event-window median hourly dynamic pressure ~**0.609 nPa**, versus survey-hour median ~**0.527 nPa**.
- Because MAVEN changes altitude and magnetic-region geometry, these comparisons are **descriptive context**, not tests of an upstream interplanetary shock, and do not independently confirm the CME's identity.

## Frozen automated selection

1. Screen **every available** SEP row in the 2025-08-01 to 2025-12-05 time series. Retain valid, finite, nonnegative flux, removing spurious >10¹⁰ fills. Require **≥3 valid directional views out of four** for each species.
2. Use the companion `sep_ion_fov*_quality` and `sep_electron_fov*_quality` fields as **1σ Poisson statistical flux uncertainties**, **not a binary 0-good data-quality flag** (the third-party dataset card erroneously describes these as quality flags; refer to primary MAVEN SEP/KP documentation). Require flux/uncertainty **≥8** in **≥3/4 views of each species**, with per-sample median ion ≥400 and electron ≥1,000.
3. Require ≥100 observed samples per UTC hour, ≥75 qualifying joint samples, ≥60% of that hour jointly qualifying, and ≥20 qualifying samples at spacecraft altitude ≥3,000 km.
4. Link selected hours if separated by ≤5h, require ≥2 selected hours. Count distinct MAVEN SPICE orbits; **report but demote single-orbit cases**. Apply the same selection to the entire survey; no event-specific threshold changes.
5. Preserve original sensor UTC timestamps; use nearest SPICE geometry within 9 s, not interpolated artificial measurements. Coarsen independent MAG/SWIA context into hourly medians, not falsely synchronous 8-second claims.

This is a high-significance screening rule (especially for repeated, large amplitudes), not a calibrated global false-positive probability. There are multiple selection criteria and search across all dates. Further independent data are needed for causal claims.

## Independent primary NASA PDS4 product verification — **PASSED**

After the initial mirror-based survey and forecast match were locked, an independent source-file replay retrieved the original NASA PDS **MAVEN in-situ KP collection 41.0**, product **`mvn_kp_insitu_20250824_v26_r01`** (created June 10, 2026). This is a **NASA-produced, labeled, fixed-record 45,698,463-byte file**, not a repackaging supplied by the third-party mirror.

- The downloaded primary `.tab` **matched NASA's published MD5 `377f757a180b34c7af46b2937831c181`**, byte-for-byte identity of the primary archive file. Primary PDS SHA-256: `ad8807ac1922dfbb92a806cd49af8ae1a02cab7a18f46030129ec98395a4f723`.
- The original `.xml` PDS4 label was separately downloaded and parsed to derive the byte offsets and widths of the **8 SEP species/directional-view flux fields and 8 statistical uncertainty fields**. No numeric column positions were taken from the mirror.
- The NASA PDS file contains **12,106** Aug 24 records; the rehosted Parquet contains **12,105** records that date (one-row day-level coverage difference, not a flux disagreement).
- During the decisive **05:00–05:59 UTC** peak hour, **450/450 timestamps** matched between primary NASA PDS and the independent rehosting, with **450 positive matched measurements in each of all 16 channels**.
- The **median and 95th-percentile absolute fractional differences were exactly zero for each of all 16 flux/uncertainty channels**. This independently confirms the primary NASA archive amplitude and timestamp provenance of the peak and the correctness of interpreting `*_quality` as numeric uncertainties.

[Primary NASA PDS product detail and original `.tab`/`.xml`](https://search-pdsppi.igpp.ucla.edu/search/view/?f=no&id=pds%3A%2F%2FPPI%2Fmaven.insitu.calibrated%2Fdata%2F2025%2F08%2Fmvn_kp_insitu_20250824_v26_r01&o=1) · [Machine-readable original-archive audit, official hash and per-channel comparisons](primary_pds_audit.json)

The independent download and comparison are reproducible via `examples/maven_pds_kp_primary_audit.py`, `tests/test_maven_pds_kp_primary_audit.py` and [the dedicated primary PDS check workflow](../../.github/workflows/maven-pds-kp-primary-audit.yml). The audit passed its PDS parser tests, checksum, source-schema validation, and numeric cross-check and published its result. **This is cross-validation of NASA's original KP data representation, not a second independent instrument or the separately calibrated 2026-10-30 SEP L2 data**. It validates the *measured peak*, not unique CME association, exact impact time, surface exposure, or all survey days.

## Primary sources, hashes, reproducibility

**Primary spacecraft archive:** NASA MAVEN Science Data Center at University of Colorado LASP, https://lasp.colorado.edu/maven/sdc/public/data/sci/kp/insitu/ . The analysis used an openly available **2026-updated third-party Parquet conversion of NASA LASP MAVEN KP**, https://huggingface.co/datasets/juliensimon/nasa-maven-kp-insitu (updated **2026-10-01**, as of this analysis), with the following exact input-file SHA-256 checksums:

| File | SHA-256 |
|---|---|
| `sep.parquet` | `e3ea044212db893dcdcc050c5b1aa28257b05f483c9952258bc715e9631dec03` |
| `spice.parquet` | `997fa7e3bce857722232531fa18beed4a843eca27aa0ecde910ca01ee5bf1ee0` |
| `mag.parquet` | `45b2e91ead4d5aba3aa68a2a2170f4f224681f3c47117c011aa65b8f20b8ce8c` |
| `swia.parquet` | `186c741f2f04938ce46ebcfe3d2886f77f24b9b794b9872fc3d859268ac05a87` |

Sensor-data source attribution: **NASA/University of Colorado LASP MAVEN mission**; mirror curator: `juliensimon`. File receipts and byte counts are in `summary.json`. There is a newer **2026-10-30** dedicated *PDS calibrated SEP L2 archive release*, but **those Level-2 files were not used to produce these numbers**; they are a pending, genuinely independent product-level follow-up. The mirror publication/refresh date does not mean the actual spacecraft measurement date was 2026.

To repeat independently:

```sh
python -m pip install numpy==2.1.3 pandas==2.2.3 pyarrow==17.0.0
python -m unittest discover -s tests -p test_maven_sep_final_season.py -v
python examples/maven_sep_final_season.py --output /tmp/maven-summary.json --daily /tmp/maven-daily.csv
```

The branch workflow [MAVEN final-season four-look particle-event screen](../../.github/workflows/maven-final-season-sep.yml) successfully ran the **real full-data scan**, tests, and published these compact results plus the [machine-verifiable independent NASA DONKI forecast comparison](forecast_alignment.json) with source hashes. Other repo investigations and earlier measurements are unmodified.

## Critical falsifiable next measurements

1. **Independent calibrated SEP L2:** use the PDS late-2026 release and both sensor telescope directions to check the August 24 absolute flux peaks, species, energy bands, data flags, and telescope orientation against the KP conversion. If they disagree, retract the amplitude claim.
2. **In-situ shock detection:** look near **2025-08-23 22:00 to 2025-08-24 12:00 UTC** for independently timed, sharp, physically connected MAG field, SWIA bulk speed and dynamic pressure discontinuities, with upstream-only geometry from SPICE. SEP onset can precede a shock; absence of a shock signature weakens the simple impact explanation.
3. **Distinguish Aug 20 and Aug 21 CMEs:** trace plasma arrival and field directions using NASA DONKI model files, spacecraft heliolongitude, multiple eruption sources; do not force a unique identity from time concordance.
4. **Independent observations:** compare MAVEN IUVS high-altitude auroral data, Mars Odyssey / MSL RAD energy-resolved radiation signatures where observation windows overlap, or appropriate Mars surface records. Cross-mission detections must be independently calibrated in physical units.
5. **Global null test:** predeclare all forecasted Mars hits and matched non-hit dates in the entire observed window. Quantify coincidence enrichment, controlling solar activity, orbital altitudes, and the missing 28 days, before reporting a significance or causal risk.

**Bottom line:** A real and very large dual-species MAVEN enhancement, measured by this new, reproducible full-archive analysis, coincides with *prior* NASA Mars-directed CME forecasts. The quantitative measurement and forecast association add useful evidence, but the physical identity of the driver and any downstream effects remain unconfirmed.
