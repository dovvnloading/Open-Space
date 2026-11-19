# 4XMM-DR14 QPE search — SDSS host audit, failed gated pn-screen pilot

**Current result (19 November 2026):** Known QPE control recurrence remains **6/6 active XMM epochs** with **29/29 published bright-event peak-count overlap (in-sample; not sensitivity)**. A full **33,488-row** soft/long X-ray prefilter was queried against SDSS DR18: **2,510** clean galaxy-centroid associations within 2 arcsec. A published-MD5-verified, pointlike long-pn catalogue gate froze **1,108 provisional detections in 742 XMM observations**. **Three** pre-registered XSA PPS observations failed at transport; **zero new-source pn light curves were screened**. No independently confirmed new QPE, vetted candidate, calibrated false-positive rate or archive-wide null. [PR #205](https://github.com/dovvnloading/Open-Space/pull/205) remains **DRAFT, not merged**.

## Measured positive-control recovery

| Control | QPE-active XMM obsID | Published bright eruptions | Recovered soft peaks | Recovered peak separations (hours) |
| --- | --- | ---: | ---: | --- |
| GSN 069 | 0851180401 | 5 | 5 | 8.333, 9.333, 8.583, 9.333 |
| RX J1301.9+2747 | 0864560101 | 8 | 8 | 1.333, 5.000, 2.917, 5.417, 3.250, 5.583, 3.917 |
| eRO-QPE1 | 0861910201 | 2 | 2 | 11.583 (second flare is partially observed) |
| eRO-QPE2 | 0872390101 | 9 | 9 | 2.333–2.500 |
| eRO-QPE3 | 0883770101 | 2 | 2 | 20.417 (published ~20.4) |
| eRO-QPE4 | 0883770401 | 3 | 3 | 9.750, 14.750 (published 9.82, 14.70) |
| **Total** | **6 epochs** | **29** | **29** | **6/6 controls show recurring peaks** |

The six-epoch source-level **operational recovery is 6/6 (100%)**. Bright-event count overlap is **29/29 (100%)** as a capped, count-only upper bound, **not a measured sensitivity or individually proven 29-event match**. For RX J1301.9+2747, the independently tabulated ten 2020 fitted/partial peak timestamps yield **8/10** matches within ±1200 seconds; the bright-event benchmark is eight. For eRO-QPE1 the measured 11.583-hour within-epoch separation must not be mistaken for the later NICER mean ~18.5-hour recurrence. Literature and quantitative caveats are in [control recovery audit](outputs/CONTROL_RECOVERY.md), [machine-readable replay results](outputs/control_replay.json) and [control gate](outputs/control_gate_replay.json).

## Method, frozen in-sample for control replay

Source-specific XMM PPS merged EPIC-pn `PNX000SRCTSR8###.FTZ` from the XSA, with FITS `SRC_RA/SRC_DEC` crosschecks (0.67–3.27 arcsec across all six). Source soft bands `RATE1+RATE2+RATE3` form the **0.2–2 keV background-subtracted NET** count rate; `BACK1V+BACK2V+BACK3V` form the local scaled background. Exposure-weighted 300-second bins preserve GTI gaps. A fixed 3.5-sigma single-bin + 7-sigma multiscale (450/850/1300 s) peak criterion, minimum 3x baseline and scaled-background flare veto are documented in [control thresholds](control_thresholds.json), implemented in [control_validation.py](control_validation.py); [replay_controls.py](replay_controls.py) reruns saved binned curves without XSA network requests. The 120 MiB XSA budget that caused the pilot control failures was raised to a bounded 3 GiB per target with at most three downloads per 55-minute Actions job, SHA-256 receipts per product and stream; eRO-QPE3's PPS source index was repaired through FITS positional selection rather than assuming the older DR14 source number.

The branch continues to use the existing push-triggered [`.github/workflows/xmm-qpe-dr14.yml`](../../.github/workflows/xmm-qpe-dr14.yml) on `research/xmm-qpe-dr14`; source changes are separate batches, and output-only result commits do not recurse. The last job performs a no-network six-control replay. The historical [original pilot outputs](outputs/summary.json) demonstrate only feasibility; **do not interpret their failed control checks as a negative result.**

## Gated archive host-screen attempt and referee verdict

### Auditable cohort freeze

- **33,488** eligible soft/long 4XMM detections (original X-ray-only prefilter; not galaxy nuclei) were each queried in 419/419 successful SDSS DR18 SQL batches.
- **2,510 detections / 2,519 pairs** met the SDSS PhotoPrimary morphological clean primary galaxy nucleus 2-arcsec positional requirement; 9 detections had ambiguous galaxy matches.
- A checksum-verified 4XMM-DR14 detection FITS, fetched in 60 HTTP 206 range chunks (MD5 47c6bd113cb8510db418665472c81b48), froze **1,108** point-like long-EPIC-pn, PN-photon-count-qualified detections in **742 unique observations** and **897 unique X-ray SRCIDs**. The filter requires a unique SDSS galaxy match, EP_EXTENT<=6 arcsec, EP_EXTENT_ML<4, PN_ONTIME>=21600 s, PN_8_CTS>=250, and clean critical positions in the 12-character PN_FLAG, in addition to the existing soft-hardness/quality cuts. An early zero arising from misparsing PN_FLAG as numeric was **rejected and corrected**. [Receipt](outputs/archive_pn_point_freeze.json).
- **3 pre-registered PPS observations** (0872392901, 0306630101, 0200730401) were tried in one bounded runner job; all requests failed as ESA XSA closed the connection. **0/742 observations and 0/1,108 new-source pn detections** were actually screened. Zero observed new flags is not an archive null. [PPS failure receipt](outputs/archive_pilot_batch_01.json).

### Controls, injection and background

- The original frozen soft-band NET PN finder and its threshold JSON are **unchanged**. Six known QPE epochs pass the published recurrence-range control; bright-event 29/29 *count-only*, not all individual published event timestamps.
- Two-Gaussian quiet-series injections on those same six known-QPE controls recovered both injected event times in **240/324** trials. This is rate-level, not event-list injection or a calibrated survey sensitivity function.
- Shuffling the actual eruption-containing controls yielded **210/210** recurrent flags; this is **not an independent astrophysical negative-control sample or survey false-positive rate**. The scaled in-source BACK1-3 series gave **0/6** recurrence flags, but no independent off-source region or particle-flare GTI was verified. [Diagnostic details](outputs/archive_statistical_checks.json).
- The full Legacy Surveys sky, Gaia proper-motion/star veto, 4LAC/SIMBAD blazar clearance, pn event-level pile-up and particle-background screening, non-QPE AGN/TDE/star negative controls, eRASS1/Swift/Chandra and TNS/ATel/literature candidate vetting are **not completed**.

**Independent methodological referee decision: KEEP DRAFT, DO NOT MERGE.** No validated new QPE and no [REMARKABLE] marker or docs/REMARKABLE_FINDINGS.md entry. The remaining 742 provisional XMM observations have no validated new-source pn curves, FDR and sensitivity are undefined, and an archive null or discovery would be unsupported. Full numerators, hashes, run IDs, and blockers are in the [archive audit](outputs/ARCHIVE_SCREEN_AUDIT.md) and [referee report](REFEREE.md).

**Control-scope exclusions:** AT2019qiz and Ansky do not have qualifying QPE-active DR14-era XMM control epochs. They are not scored as failures.

## Reproduction/evidence

- [Published-count benchmark](https://doi.org/10.1093/rasti/rzad015) (Webbe & Young 2023, Table 1); [eRO-QPE3/4 times](https://doi.org/10.1051/0004-6361/202348881) (Arcodia et al. 2024); [RX fitted peak timing](https://doi.org/10.1051/0004-6361/202450861) (Giustini et al. 2024).
- [Control batches 1–4](outputs/CONTROL_RECOVERY.md), [PPS receipts](outputs/control_receipts_batch_03.json), [replay](outputs/control_replay.json), and [strict frozen selection](outputs/frozen_selection.json).
- [Final host/pn gate](outputs/archive_pn_point_freeze.json), [SDSS crossmatch receipt](outputs/archive_sdss_receipt_00000_33488.json), [surrogate injection/shuffle controls](outputs/archive_statistical_checks.json), [PPS pilot failure](outputs/archive_pilot_batch_01.json), [machine admission decision](outputs/archive_screen_gate.json), and [full runner audit](outputs/ARCHIVE_SCREEN_AUDIT.md).
- [PR #205](https://github.com/dovvnloading/Open-Space/pull/205) for run history and explicit draft status.
