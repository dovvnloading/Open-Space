# Archival XMM-Newton QPE search: 4XMM-DR14 feasibility pilot

**Result (reviewed, 19 November 2026): INCONCLUSIVE. No new nuclear QPE identified; no calibrated null established.** GitHub Actions [run 44470633598](https://github.com/dovvnloading/Open-Space/actions/runs/44470633598) succeeded at 07:20 UTC and committed real data receipts and screening outputs. That operational success **is not scientific validation**. A prior attempt, [run 44468260187](https://github.com/dovvnloading/Open-Space/actions/runs/44468260187), timed out.

## Fixed observations and provenance

- Full **4XMM-DR14 v1.0** detection FITS from the HEASARC fixed-release URL: HTTP 200, 1,005,466,736 downloaded bytes; published MD5 `47c6bd113cb8510db418665472c81b48` verified; downloaded SHA-256 `6d8703ad0ac3b8f737e1c293651905c2c7a4a1007c430add67d4cdf95b916597`.
- Catalogue: **1,035,832** detection rows, **372,313** rows flagged TSERIES. Registered soft/long prefilter left **33,488** detections, frozen in [frozen_detections.csv](outputs/frozen_detections.csv) with SHA-256 `c2793211c74f9a85992253fba76743a94bcdf2cffb2137f29a3037907d4e56a5` (digest recorded in the run manifest; not independently redownloaded in this review).
- Screen cuts: TSERIES, SUM_FLAG ≤ 1, EP_8_CTS ≥ 250, at least one camera ontime ≥ 21,600 s, EP_HR3 < −0.3, EP_HR4 < −0.4. This is an **X-ray spectral/exposure** selection, not a galaxy-nucleus selection.
- Selected first **eight ranked detections** plus four catalogue-matched known-source controls, at most 12 PPS requests. **9/12** requests yielded parsable EPIC-pn source light curves. The three unavailable products were truncated at the job's **120 MiB per-request** limit. All 12 download attempts, including three repeated downloads of observation 0011420101 for separate sources, are recorded with timestamps, URLs, byte counts, hashes or explicit errors in [receipts.json](outputs/receipts.json).
- Source scripts and [frozen_config.json](frozen_config.json) are preserved. The complete [frozen_selection.json](outputs/frozen_selection.json) contains some bare `NaN` values and is **not strict JSON**; this is an identified reproducibility defect to fix before a stronger release. The observed catalogue and PPS downloads are not permanently vendored; reruns depend on their public archive endpoints.

## Temporal screen and all flagged nuclear candidates

The eight ranked sources were 4XMM J144207.4+352622, J114323.6+195817, J123600.5−395215, J161742.0−833750, J051541.4+010441, J064059.8+095149, J064048.1+095339, and J064022.6+094946. The pilot parsed EPIC-pn broad-band light curves for **all eight**; none triggered the scripted repeated-peak detector. There were **zero unvetted repeat flags**, not zero QPEs in the parent population. [Light-curve checks](outputs/lightcurve_checks.json) · [Empty flags](outputs/unvetted_flags.json).

**There are no candidate nuclei with validated repeated soft flares to report.** No source has a measured QPE recurrence period, peak-to-quiescent amplitude, or thermal-temperature evolution from this run. EP_HR3/HR4 catalogue cuts are not temperature measurements; a parsed broad-band flare flag is not evidence of a nuclear QPE. None of the eight was cross-matched to a Legacy Surveys/SDSS/2MASS galaxy nucleus or passed Gaia stellar rejection.

## Known-QPE controls: catalogue presence ≠ flare recovery

| Published control | DR14 position-matched detections (within 12″) | Pilot PPS test | Actual QPE-flare recovery |
| --- | ---: | --- | --- |
| GSN 069 | 33 | Obs 0851180401: 120 MiB PPS cap exceeded | **No** |
| RX J1301.9+2747 | 5 | Obs 0864560101: 120 MiB cap exceeded | **No** |
| eRO-QPE1 | 2 | Obs 0861910201 parsed, 25.264 h light-curve span; detector found zero peaks | **No — known-source false negative / sensitivity unverified** |
| eRO-QPE2 | 3 | Obs 0872390101: 120 MiB cap exceeded | **No** |
| eRO-QPE3 | 2 | Position matches only; not in capped PPS requests | **Not tested** |
| AT2019qiz | 0 | No 12″ position match at supplied coordinates; no PPS retrieval | **Not tested** |
| eRO-QPE4 | Not queried | Position not in frozen control table | **Not tested** |
| Ansky (SDSS1335+0728) | Not queried | Position not in frozen control table | **Not tested** |

**Positive-control flare recovery is 0 validated, not a measured recovery fraction.** Four named controls were requested, but only eRO-QPE1 yielded a parsed light curve; even that source's known QPE activity was not recovered. No inference of QPE absence can be made from the other download failures, limited epochs, or from a detector that fails its positive control.

## Prior literature / novelty audit

A targeted cross-check verifies that named controls were **already known in the literature**, not newly discovered in this run: GSN 069, RX J1301.9+2747 and eRO-QPE1–2 are discussed in [Giustini et al. (2024)](https://doi.org/10.1051/0004-6361/202450861); eRO-QPE3 and eRO-QPE4 are published in [Arcodia et al. (2024)](https://doi.org/10.1051/0004-6361/202348881); [AT2019qiz's QPEs](https://www.nature.com/articles/s41586-024-08023-6) were published in 2024; [Ansky's QPEs](https://doi.org/10.1051/0004-6361/202555258) were published by 2025. These are **external prior reports**, not recovered flares in our output. Further sources, e.g. eRO-QPE5 and J2344, are discussed in [later literature](https://doi.org/10.1051/0004-6361/202558531) and were not part of this frozen control set. **No fresh candidate passed the time-series screen, hence no candidate-specific ATel/TNS/AGN/TDE prior-art clearance was performed.** It must not be represented as a completed literature exclusion of new QPEs.

## Referee verdict and blockers

**REJECT as a QPE discovery or calibrated null; retain as an inconclusive pilot only.** [Independent-referee audit](REFEREE.md) documents: no validated known-QPE flare recovery; 3/4 requested controls missing or truncated; untested eRO-QPE3/4, AT2019qiz and Ansky; no background flare/GTI/pile-up vetting or soft-band spectral fits; no host identification, Swift/Chandra/eRASS1 cross-epochs, injections, or false-positive rate. **Recovery fraction and false-positive rate are undefined**, not zero. Strict JSON serialization of the frozen selection requires repair.

A future substantive search must use source-specific PPS retrieval without the 120 MiB whole-archive cap, energy-resolved source/background light curves and scientifically motivated time binning; show actual published-control flares with measured recurrence and amplitude; then freeze a broad nuclear-host sample, inject flares, test stellar/background controls, and vet any source against prior announcements and other X-ray epochs.

**Classification:** unclassified feasibility pilot; not a measured null and **not [REMARKABLE]**. The findings policy in `docs/REMARKABLE_FINDINGS.md` does not permit an entry.

## Files / replay

- [Run summary](outputs/summary.json) · [machine receipts](outputs/receipts.json) · [full frozen list](outputs/frozen_detections.csv) · [selected sample and controls](outputs/frozen_selection.json) · [EPIC-pn results](outputs/lightcurve_checks.json) · [run report](outputs/RUN_REPORT.md).
- `search.py` runs as **one network-enabled GitHub Actions job**, triggered only by pushes to the research branch affecting science code or frozen configuration under this result folder. Generated outputs and documentation-only changes do not re-run it. Scientific dependencies are pinned in `.github/workflows/xmm-qpe-dr14.yml`.
