# 4XMM-DR14 search for X-ray quasi-periodic eruptions — control-only calibration

**Current result (19 November 2026):** Known-QPE soft-X-ray recurrence **recovered in six of six selected QPE-active XMM observations**, with **29 of 29 published bright-eruption peak counts overlapped**. This is **not** independently verified flare-level completeness, **not** a new QPE detection, **not** an archive-wide screen of nuclear galaxies and **not** a calibrated null. Draft [PR #205](https://github.com/dovvnloading/Open-Space/pull/205) must **remain draft; do not merge.**

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

## Referee decision: not ready for archive survey or merge

The original frozen 4XMM-DR14 prefilter has **33,488 soft/long X-ray detections**; these have **not** been crossmatched to an optically selected, completeness-controlled census of **galaxy nuclei**. No post-control candidate screen has therefore run, and **zero new candidate QPEs have been independently checked**. It would be inaccurate to report no new QPEs *in the catalogue*: the archive has not been searched with validated nuclear-host selection.

**Required before screening or an eventual merge:** independently align the full set of published eruption timestamps (especially eRO-QPE1's incomplete second peak), validate particle flare GTIs and separate source/background light curves, pile-up and soft spectral/hardness evolution, negative stellar and non-QPE AGN controls, injection/recovery and a false-alarm distribution; freeze host crossmatching and survey completeness; then screen archived nuclear detections using the frozen thresholds and vet candidates against Gaia, host offsets, XMM background variability, TNS, ATel, published QPEs and transient/AGN literature. See [independent referee report](REFEREE.md) and [control recovery audit](outputs/CONTROL_RECOVERY.md).

**Coverage exclusions:** AT2019qiz and Ansky have no qualifying QPE-active XMM epoch in this frozen DR14-era control setup; they are **not** counted as false negatives. **No** `[REMARKABLE]` label, `docs/REMARKABLE_FINDINGS.md` entry or merge commit is warranted.

## Reproduction/evidence

- [Published-count benchmark](https://doi.org/10.1093/rasti/rzad015) (Webbe & Young 2023, Table 1); [eRO-QPE3/4 times](https://doi.org/10.1051/0004-6361/202348881) (Arcodia et al. 2024); [RX fitted peak timing](https://doi.org/10.1051/0004-6361/202450861) (Giustini et al. 2024).
- [Control batches 1–4](outputs/CONTROL_RECOVERY.md), [PPS receipts](outputs/control_receipts_batch_03.json), [replay](outputs/control_replay.json), and [strict frozen selection](outputs/frozen_selection.json).
- [PR #205](https://github.com/dovvnloading/Open-Space/pull/205) for the run history and the explicit draft review status.
