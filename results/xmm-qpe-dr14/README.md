# 4XMM-DR14 QPE search — SDSS host audit, failed gated pn-screen pilot

**Current result (19 November 2026): bounded archive screen, no new QPE.** The frozen 0.2–2 keV EPIC-pn finder recovers known-QPE recurrence in **6/6** control epochs (29/29 published bright peak-count overlap; in-sample, not a measured completeness function). Of **33,488** soft/long 4XMM-DR14 detections, **2,510** matched clean SDSS galaxy centroids within 2″, a point-source/long-pn gate kept **1,108**, SDSS spectra removed **119** QSO/star/blazar rows and Gaia DR3 astrometry removed **22**, leaving **967 detections in 670 XMM observations**. From the NASA HEASARC mirror, **910 source light curves in 637 observations** were screened; **33 observations (57 detections) have no verified pn curve and are not negatives**. The screen raised **6 flags**: four are the published QPE RX J1301.9+2747 recovered blind in four epochs and one is a known Seyfert 2 during a particle flare. The sixth (WISEA J134144.52−005832.5, z = 0.148, obs 0841482101) is **unconfirmed, low-significance variability**: two of its three nominal peaks coincide with elevated particle background and the third is about 3σ in an independent soft-band check ([flag dispositions](vetting/README.md)). **No new QPE is identified.** The sample is limited to SDSS-footprint hosts; the Legacy Surveys complement, event-level pile-up and off-source particle controls, and a calibrated sensitivity or false-positive rate were not completed, so this is not an archive-wide occurrence limit.

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
- The full Legacy Surveys sky, complete Gaia epoch-aware stellar vetting, 4LAC/SIMBAD blazar clearance, pn event-level pile-up and particle-background screening, non-QPE AGN/TDE/star negative controls, eRASS1/Swift/Chandra and TNS/ATel/literature candidate vetting are **not completed**.

**Independent methodological referee decision: KEEP DRAFT, DO NOT MERGE.** No validated new QPE and no [REMARKABLE] marker or docs/REMARKABLE_FINDINGS.md entry. The final 670 provisional SDSS/Gaia-observed XMM observations have no validated new-source pn curves, FDR and sensitivity are undefined, and an archive null or discovery would be unsupported. Full numerators, hashes, run IDs, and blockers are in the [archive audit](outputs/ARCHIVE_SCREEN_AUDIT.md) and [referee report](REFEREE.md).

**Control-scope exclusions:** AT2019qiz and Ansky do not have qualifying QPE-active DR14-era XMM control epochs. They are not scored as failures.


### Additional SDSS contaminant-only audit (not a complete host vet)

A later bounded, successful SDSS DR18 SpecObjAll and near-star crossmatch of the 1,108 detection rows found 205 unique hosts with spectra (419 spectroscopic rows), and **119 source detections** with spectroscopy classified as QSO/STAR or recognized BL Lac/blazar subclass; it found **0 nearby clean SDSS stars within 1 arcsec**. This excludes the 119 rows conservatively and leaves **989 provisional detections, 795 unique X-ray SRCIDs and 688 unique XMM observations** in the [post-exclusion manifest](outputs/archive_sdss_after_exclusions.csv). The run made all SQL chunks successfully and retried PPS access for the three pilot observations plus a known QPE control; all four ESO/XSA connection attempts failed again. **No Gaia proper-motion/parallax, complete SIMBAD/4LAC blazar rejection or new pn light curves were obtained.** The 1,108-to-989 count is an SDSS-only reduction, *not* a validated galaxy-nucleus survey or screening result. [Independent exclusion receipt](outputs/archive_sdss_stellar_blazar_audit.json).

### Gaia DR3 astrometric veto — completed on SDSS survivors, still provisional

A full **25/25 successful Gaia DR3 TAP upload batch audit** checked **989/989** remaining post-SDSS sources within 1 arcsec of their X-ray positions, with **408** Gaia positional counterparts. **22 distinct X-ray detections** met a conservative Gaia contamination flag (positive parallax >=5σ or at least one proper-motion component >=5σ). Excluding those flags yields **967 detections, 784 SRCIDs, 670 XMM observation IDs** in a frozen, still-provisional [SDSS+Gaia host CSV](outputs/archive_provisional_nuclei_sdss_gaia.csv) with a [machine-readable checksum/gate](outputs/archive_final_host_cohort.json). These flags are not conclusive classifications of every associated host: astrometry of extended galaxies can be spurious, epoch-dependent offsets are unresolved, and full 4LAC/SIMBAD blazar/AGN and Legacy Surveys coverage remain absent.

The 3 PPS retrieval failures preceded Gaia and did not return any new-source light curves; hence **0 of 670 final provisional observation IDs** were screened with the frozen peak finder. The [runner audit](outputs/ARCHIVE_SCREEN_AUDIT.md) keeps the background-region, pile-up, event-list injection and FPR caveats explicit. **Referee gate remains closed.**

## Reproduction/evidence

- [Published-count benchmark](https://doi.org/10.1093/rasti/rzad015) (Webbe & Young 2023, Table 1); [eRO-QPE3/4 times](https://doi.org/10.1051/0004-6361/202348881) (Arcodia et al. 2024); [RX fitted peak timing](https://doi.org/10.1051/0004-6361/202450861) (Giustini et al. 2024).
- [Control batches 1–4](outputs/CONTROL_RECOVERY.md), [PPS receipts](outputs/control_receipts_batch_03.json), [replay](outputs/control_replay.json), and [strict frozen selection](outputs/frozen_selection.json).
- [Final host/pn gate](outputs/archive_pn_point_freeze.json), [SDSS crossmatch receipt](outputs/archive_sdss_receipt_00000_33488.json), [surrogate injection/shuffle controls](outputs/archive_statistical_checks.json), [PPS pilot failure](outputs/archive_pilot_batch_01.json), [machine admission decision](outputs/archive_screen_gate.json), and [full runner audit](outputs/ARCHIVE_SCREEN_AUDIT.md).
- [PR #205](https://github.com/dovvnloading/Open-Space/pull/205) for run history and explicit draft status.
