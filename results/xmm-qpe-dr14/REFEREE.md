# Independent referee report — XMM QPE DR14, PR #205 (2026-11-19)

**Recommendation: KEEP DRAFT. No merge, no new-QPE claim, no population-null claim.** The predecessor pilot was scientifically invalid as a search; the repaired branch is now a **working calibration demonstrator** for six positive controls, not an archive survey.

## What actually passed

1. All six QPE-active XMM EPIC-pn controls are downloaded from ESA PPS without hitting the former 120 MiB cap. Each actual source has been identified by `SRC_RA/SRC_DEC` from the FITS header, with offsets **0.67–3.27 arcsec**. eRO-QPE3's true PPS identifier is 7 rather than historical DR14 source 5.
2. Soft 0.2–2 keV band `RATE1+RATE2+RATE3` with `HDUCLAS2=NET` is background-subtracted; background estimator `BACK1V+BACK2V+BACK3V` and fractional exposure are tracked. 300s GTI-aware rebins and fixed source-only multi-width peak templates recover the control peak trains.
3. Measured source-level recurrence-range recovery **6/6**; published bright eruption peak-count overlap **29/29** (count-only upper bound; in-sample). RX J1301.9+2747 **8/10** specific Giustini 2024 fitted peak times matched within 1200s. eRO-QPE3 **20.417h** vs ~20.4h; eRO-QPE4 **9.75/14.75h** vs **9.82/14.70h**. Full table and per-target evidence in [the control audit](outputs/CONTROL_RECOVERY.md), [the replay](outputs/control_replay.json), and [threshold specification](control_thresholds.json).
4. Fixed 3 GiB/14-minute per target, maximum three per under-55-minute Actions job; raw tar/product SHA-256 receipts. Synthetic source and correlated-background flaring fixtures pass. Earlier failed replay report output key was corrected; last complete replay was successful.
5. eRO-QPE1's second source-associated flare is only partially observed; its 11.583h within-epoch separation is **not** a measurement of the published ~18.5h NICER mean cycle.

## Rejection criteria still failed

- **Not independent validation:** the finder is tuned on the same six positive controls; the 29/29 count match is **not** 29/29 timestamp-verified sensitivity, and two RX fitted weak/partial eruptions remain unmatched (the first is MOS-only). Some broad recurrence windows are permissive.
- **Background and pile-up incomplete:** the PPS NET series and embedded BACK1–3 are useful, but no separately verified filtered high-energy particle/FBKTSR background GTI and no `epatplot` pile-up/mos confirmation or independent spectral fit/hardness evolution test. Synthetic flare rejection alone is not enough.
- **No negative controls or calibrated trials:** no blank background, variable Galactic stars, TDE/AGN negatives, injection-recovery curve, multiple-testing correction/FDR, or held-out threshold benchmark.
- **No validated *nuclear* parent sample:** 33,488 DR14 soft/long detections are X-ray detections. Galaxy host/nuclear offsets, redshifts, counterpart colors/Gaia and star rejection, per-obs completeness and field-coverage denominators have not been assembled; no meaningful archive screen has occurred after control validation. The old twelve attempted downloads (nine parsed) cannot substitute for this.
- **No candidate vetting:** zero new postgate candidates examined; no optical/IR host assessment, pile-up or background-flare check, literature/TNS/ATel novelty evidence for a newly selected source. It would be improper to claim zero new QPEs *in the archive* from an unexecuted survey.
- **Control scope exclusions:** AT2019qiz/Ansky QPE activity is outside the qualifying DR14-era XMM control epochs and cannot be scored as misses (nor as hits).

## Decision rule

Remain draft until (a) every control's published eruption times and expected interval are matched with literature-cited uncertainty, with full exposure/window accounting; (b) independently verified particle background, pile-up, energy evolution and negative controls establish low false-positive rate and sensitivity; (c) optical galaxy-center host sample and frozen thresholds are defined; (d) full archive screen is batched under 55-minute jobs with receipts; and (e) every alleged candidate is vetted with cross-matches, TNS/ATel/literature and artifact checks. If no candidate passes those tests, report a bounded survey null **with sensitivity and exposure denominator**; if a genuinely novel QPE passes, only then consider `[REMARKABLE]`. None of (a)–(e) currently holds in full.

**PR #205 stays draft. Merge SHA: NONE.** Original 33,488-entry catalogue prefilter and failed pilot remain auditable history, not a completed survey.
