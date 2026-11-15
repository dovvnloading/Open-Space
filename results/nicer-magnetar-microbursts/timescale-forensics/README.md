# Time-scale forensics: new holdouts and a discordant earlier burst table

**Scope:** This is an additive follow-up to [PR #23](https://github.com/dovvnloading/Open-Space/pull/23) and its
[catalog-wide NICER clock test](../catalog-clock/README.md). No files in that
analysis are altered. The new measurements here are (1) disjoint *calendar-year
strata* of the existing event-based count table and (2) a previously unused
ten-burst comparison to Shao et al. (2024), Table 4. The second comparison
**does not agree** with the frozen PR #23 local absolute-time reconstruction;
this is a potential competing clock/convention problem, not confirmation of
a second paper having the same 69.184-second mistake.

## 1. Trace the time reference and separate competing mechanisms

Chu et al. (2026), Section 2.1, state that `nicerl2` was followed by
`barycorr` with JPL DE430. Their Table A1 has a column literally called
**TDB Start**; the machine-readable table has 1,132 rows, although the paper
reports 1,130 bursts. Their Figure 3 labels 2018-11-22T11:58:14.425441 and
2022-07-24T23:59:07.273212 as TDB. The independent PR #23 reconstruction
placed the corresponding event peaks 69.184 seconds later.

[NASA's NICER timing guide](https://heasarc.gsfc.nasa.gov/docs/nicer/analysis_threads/time/)
defines the following for *calibrated unbarycentered* event files:

- `TIMESYS='TT'`, `TIMEREF='LOCAL'`, `TASSIGN='SATELLITE'`.
- `MJDREFI=56658`; `MJDREFF=0.000777592592592593` day = **67.184 s**,
  which is TT−UTC at the **2014** NICER epoch, not at the burst date.
- The documented calibrated UTC elapsed-seconds formula is
  `UTC_elapsed = TIME + TIMEZERO - LEAPINIT`, where `LEAPINIT=2` after the
  2015 and 2017 leap seconds. Hence the **2017–2025** TT−UTC offset is
  `67.184 + 2 = 69.184 s = 32.184 + 37 s`.
- For true TDB timestamps `barycorr` changes `TIMESYS` to `TDB` and
  `TIMEREF` to `SOLARSYSTEM`, updating time references, TIME, and GTIs.
  A `barytime=YES` variant adds a BARYTIME column instead of replacing TIME,
  and needs its own metadata handling. See [NASA's barycorr guide](https://heasarc.gsfc.nasa.gov/docs/nicer/analysis_threads/bary_corr_phase_resolved_spectra/).

Two mechanisms fit an apparent +69.184 s shift **equally well over the
observation years**: a TDB instant rendered to a UTC calendar string but
labelled TDB, or an inappropriate hard-coded subtraction of 69.184 seconds at
the final reporting stage. A TT-as-UTC mix-up *after* the appropriate
barycentric correction is another possible code path. These are conjectures
about the **author pipeline**, not established descriptions of its code.
The catalog, paper and original FITS headers identify where a scientific
consumer encounters the ambiguity—**the serialized TDB Start and Figure 3
zero-points**—but not the exact generating source-code statement.

The alternatives predict different residuals: ignoring the 2014 epoch
fraction alone changes the epoch by 67.184 s; omitting LEAPINIT alone changes
local UTC conversion by 2 s; mishandling the NICER TIMEZERO keyword produces
an observation-dependent offset (often about 1 s); missing barycentering
introduces a source/date dependent delay up to hundreds of seconds, not the
observed cross-source 69.184 s. A fixed timezone error could be constant but
would be a whole-hour shift, not 69.184 s. TDB−TT contains a seasonal
term of order 1–2 ms, so **the six-decimal Figure 3 strings could in principle
distinguish a calendar-rendering error from a hard-coded constant**, once
independently reconstructed with the same source coordinates, solar-system
ephemeris and original detection definition. This investigation has **not**
resolved that millisecond distinction.

Run `python header_audit.py /path/to/verified/PR23/archive/cache` to print
the real 2018 and 2022 EVENTS and GTI keywords. The script refuses an
uncalibrated time system; it does not substitute a generic header template
for evidence from these particular FITS files. An audit of the actual
*post-barycorr* intermediate event files and publication serialization code
from the authors is still required to identify the offending program step.

## 2. New year-separated stability check (same underlying photons)

From the frozen [PR #23 paired-window file](../catalog-clock/windows.csv),
partition by year and by SGR 1935+2154 versus other targets. Retain only
rows with good-time coverage for both 0.25-s-margin windows. These rows
weren't reselected by their positive result. The analysis is *post-hoc*
stratification of the previously examined table, **not** an independently
unseen archive or a new pre-registered statistical test.

| Calendar year | Target stratum | Later window wins / eligible | Distinct visits |
|---|---|---:|---:|
| 2017 | Other | 1/1 | 1 |
| 2018 | Other | 1/1 | 1 |
| 2019 | Other | **1/2** | 2 |
| 2020 | Other | 64/64 | 35 |
| 2020 | SGR 1935+2154 | **235/332** | 7 |
| 2021 | Other | 58/58 | 30 |
| 2021 | SGR 1935+2154 | 62/62 | 16 |
| 2022 | SGR 1935+2154 | **401 wins, 1 tie, 2 losses / 404** | 19 |
| 2023 | Other | 6/6 | 4 |
| 2024 | Other | 14/14 | 6 |

The late 2024 rows (two targets, six visits) therefore provide a useful
chronologically separate check, even though the photons were already part
of the first PR. The 2019 loss, the 2020 storm's 97 losses, and incomplete
exposure remain explicit counterexamples to any claim that *every* individual
catalog row is uniquely established by a shifted aperture. Source and visit
correlations prevent treating 944 windows as independent Bernoulli trials.
Replay with `python verify.py`.

## 3. NEW independent-paper conflict — Shao et al. (2024)

Shao et al., *GTC Optical/Near-infrared Upper Limits and NICER X-Ray Analysis
of SGR J1935+2154 for the Outburst in 2022*, ApJ 976:99, Table 4, published
ten UTC burst starts from NICER ObsID **5020560106**. They precede Chu et al.
(2026). We associated ten events by the **distinctive relative spacings** of
their published starts; the pair list and full printed UTCs are in
[shao_2024_table4_pairs.csv](shao_2024_table4_pairs.csv).
[Original Shao Table 4](https://www.researchgate.net/publication/385816869_GTC_OpticalNear-infrared_Upper_Limits_and_NICER_X-Ray_Analysis_of_SGR_J19352154_for_the_Outburst_in_2022)

For these ten pairs, the Chu *printed TDB* calendar strings are **9.562 to
9.968 seconds later** than the Shao *printed UTC* calendar strings, which
is **not** a direct comparison of two like time scales. To test absolute time
without suppressing either paper's stated convention, convert the frozen
PR #23 predicted *local* MET to UTC with NASA's calibrated TT reference and
the +69.184-s shift (and LEAPINIT=2). The Shao UTC times then turn out to
be **61.228–61.495 seconds later**, with median **61.275 s**.

This is a **frame-dependent comparison, not a diagnosed second clock
error**. Shao et al. explicitly state in their Section 2.2 that they
barycentered the NICER data using the **DE200** ephemeris, while their Table 4
labels its dates **UTC**; Chu et al. used DE430 and label theirs **TDB**.
Neither table indicates unambiguously whether the final calendar string was
generated from the spacecraft-local or solar-system-barycentric timestamp.
Subtracting Shao's UTC-labelled time from PR #23's predicted *local UTC*
therefore compares quantities that may not refer to the same spatial clock
location. The >60-s difference cannot by itself establish an erroneous UTC
conversion in Shao.

We explicitly test the spatial-reference ambiguity in
[frame_controls.py](frame_controls.py): for matched event **Chu #443**,
the corrected Chu barycentric TDB is 17:41:46.293, whereas PR #23's
predicted NICER-local TT is 17:40:35.043428. The implied **SSB-TDB minus
local-TT** difference is **+71.249572 s** (approximately the Rømer delay).
Shao's UTC-labelled 17:40:27.157 is **+61.297572 s** after the predicted
NICER-local UTC 17:39:25.859428. These two offsets differ by **9.952 s**,
numerically the difference between the raw published Chu and Shao calendar
strings. This identity is algebraic: it is **not** proof that either paper
has the correct absolute timing. The TDB–TT periodic correction, precise
ephemerides, detector screening and photon-onset definition still require
a same-event, same-frame reconstruction. DE200 versus DE430 should not be
assumed to explain a 10-second residual without computing it.

The ten distinct interburst spacings strongly establish matching *events*.
They do **not** establish a second publication error. **No third independent
clock was identified here for these ten particular events.** Until the
original Shao barycentered event TIME/BARYTIME header, its UTC conversion
and the original NICER FITS rows are checked, do not promote the ~61.3-s
frame-mixed value into a correction for Shao or Chu.

## 4. Affected and unaffected claims; negative controls

| Material | Status / caution |
|---|---|
| Chu et al. Table A1 machine-readable *TDB Start* and Figure 3 absolute epochs | Directly implicated by the event/Fermi tests in PR #23. Test/relabel individual starts before a timing match. |
| Its two microsecond-burst attributions | Their *absolute timestamps* shift; a global 69.184-s shift **cannot** change the microsecond duration, the fast/slow detector-chain evidence, the 2018 particle-like coincidences, or the 2022 off-axis geometry. These remain **candidate instrumental events** rather than disproven magnetar emissions. |
| NICER/GBM coincidence searches | At least seven recovered matches from seven different visits outside the April 2020 storm appear at +69.184 s in PR #23. Not all are published missed coincidences, and the crowded storm is not diagnostic. |
| Borghese et al. (2022) SGR 1935 burst times | Two independent same-observation TDB epochs corroborate the new clock, already measured in PR #23. No evidence that Borghese inherited this catalog's error. |
| Shao et al. (2024) separate NICER UTC table | **New frame-dependent cross-paper test** above. No proof of the same offset or of a new error; different barycentric versus local clock conventions must be resolved. |
| GECAM, FRB, INTEGRAL, XMM, NuSTAR, Swift, Aql X-1 and unrelated NICER catalogs | Must be tested **individually** using matching event rows, coordinates, GTIs and the stated native clock. We have not diagnosed a shared error there. The GECAM October 2022 burst table does **not** supply a contemporaneous event for the ten Oct 12 Shao–Chu pairs. |
| FRB 200428 (2020-04-28, around 14:34 UTC); FRB 221014 (2022-10-14, around 19:21 UTC) | No Chu `TDB Start` rows on those dates within five minutes of the stated radio UTC hour. These aren't evidence of a missed radio counterpart in this specific burst list; proper barycentric conversion and coverage checks are still required. |

**Falsifiers and decisive next checks:** (i) a third, correctly referenced
independent instrument clocks these bursts at the unshifted Chu TDB epochs
after spacecraft delays; (ii) a replay of original calibrated headers,
`TIMEZERO` and `LEAPINIT` plus JPL barycorr places the event TDB times at
the published, rather than shifted, labels; (iii) the corrected offset
depends strongly on date/source (after positional errors are included);
(iv) an author-supplied intermediate TDB FITS event list has a `TDB Start`
matching the printed time and the independent photon and GBM comparisons
can be reconciled without the 69.184-s correction; (v) a future leap-second
boundary, if data were available, fails the one-second jump predicted by an
actual TT−UTC misuse. A **fixed 69.184-s literal subtraction** would *not*
jump, making the alternatives observationally degenerate in NICER's
2017–2024 analyzed span without the original serialization code.

### Reproduction

```sh
python results/nicer-magnetar-microbursts/timescale-forensics/verify.py
python results/nicer-magnetar-microbursts/timescale-forensics/frame_controls.py
# Additional header-level audit; needs astropy and the original archive files:
python results/nicer-magnetar-microbursts/timescale-forensics/header_audit.py /path/to/cache
```

Primary references: [Chu et al. 2026, ApJS 284:25](https://doi.org/10.3847/1538-4365/ae5655);
[Shao et al. 2024, ApJ 976:99](https://doi.org/10.3847/1538-4357/ad822f)
([journal listing](https://www.iaa.es/en/publications/gtc-opticalnear-infrared-upper-limits-and-nicer-x-ray-analysis-sgr-j19352154-outburst));
[NASA NICER timing guide](https://heasarc.gsfc.nasa.gov/docs/nicer/analysis_threads/time/).
