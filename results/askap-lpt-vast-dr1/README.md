# VAST Extragalactic DR1: 149 catalogue variability leads, no established long-period transient

**Result (19 November 2026, UTC): a real public ASKAP search recovered 149 compact, strongly variable 888-MHz source-catalogue entries, not 149 discovered transients.** The frozen screen encompassed **548,113 sources**, **2,945 images** and **6,449,203 Parquet measurement rows**. Sixty highest-scoring leads were followed through **693 published epoch measurements**; **14** have at least two >=8-sigma Selavy detections and one <3-sigma forced measurement. These 14 patterns merit inspection, but **no new minute-to-hour periodic radio source is validated, and this is not a calibrated null for long-period transients.** The data product has months-long time gaps and mostly 10–11 epoch flux measurements, with a few much closer observation pairs. Two detections thirteen minutes apart are not a measured period.

## Selection and evidence denominator

The search is a **new dataset**: [2026 VAST Extragalactic DR1](https://doi.org/10.25919/nh9d-t846), not a reanalysis of this repository's GLEAM-X component-spacing, stellar Stokes V, or VAST TDE studies. Code and thresholds were committed on this branch before any public source table was screened. The runner retrieved all three release tables via the authors' DAP API, preserved their SHA-256 digests, then froze the complete **2,945-image image-ID, date, field, SBID and quality-metadata list** in `output/frozen_images.jsonl` **before ranking**. The [exact manifest](output/frozen_images.jsonl), [initial threshold protocol](output/frozen_protocol.json) and [network receipts](output/network_receipts.json) are part of this result. The published image and source table counts reproduce; the measurement Parquet metadata reports **6,449,203**, exactly **one more row** than the 6,449,202 stated in the paper. We preserve the discrepancy as a metadata observation, with no interpretation as an extra astrophysical source.

An entry passes the initial source-level screen if its averaged compactness (integrated/peak) is between 0.75 and 1.40, there are at least five measurements and two Selavy detections, no component-sibling flags, a reported maximum peak SNR >=8 and peak >=5 mJy/beam, plus a positive maximum/minimum peak-flux ratio >=5 or `eta_peak >=100` together with `v_peak >=0.8`. Negative forced flux is **not** treated as an infinite flux ratio. **149/548,113** meet this heuristic; the [full 149-row selection](output/all_screened_sources.csv) is saved and ranked with no position chosen after seeing the data. This selection is *not* a sky-completeness-corrected transient search: source-catalogue construction excludes close neighbours, faint/extended sources, and failures of the upstream source association.

The [top 60](output/ranked_sources.csv) were exhaustively revisited in their [693 released measurements](output/ranked_epoch_measurements.csv), with flags, forced/Selavy origin, image IDs and epoch times retained. Independent reimplementation of the frozen gates and image-ID joins appears in the [referee report](output/referee.json). It separately identifies **14 catalogue on/off patterns** (at least two above-8-sigma Selavy catalogue measurements plus a forced measurement below 3-sigma), and records the handful of observation pairs under an hour apart. A bright low/bright-high pair can be ordinary AGN variability, refractive scintillation, image artefact or a stellar flare; none is a confirmed pulse or physically vetted transient. Per-epoch local image-noise and native FITS cutout inspections remain missing. The dynamic ranges are not claimed to be astrophysical until pixel checks.

## Controls, external checks and failure modes

**Known LPT controls:** the frozen positions of GLEAM-X J162759.5−523504.3, GPM J1839−10, ASKAP J1832−0911, ASKAP J1935+2148, ILT J1101+5521 and rounded-name CHIME J0630+25 fall outside the *approximate square image envelopes* of all 2,945 images. This is an **out-of-footprint control**, not a failed instrumental recovery and not a limit on their pulse amplitudes. [Control/geometry receipts](output/known_lpt_footprints.json). The published release specifically emphasizes extragalactic high-latitude fields, so no eligible positive control is available.

**Completed:** actual CSIRO DOI metadata and byte checksums for all three inputs, complete field/epoch list, frozen source-level screen and saved shortlist, native release-measurement rechecks, small Gaia DR3 5-arcsec cone-sample for the top 12, and a conservative independent CSV referee. The archive's legacy `files.json` endpoint returned HTTP 404 but its supported `/files` endpoint and direct downloads worked; both response paths are in the receipt. ATNF is attempted against its older URL and, if needed, the public 7 October 2026 PSRCAT release; its per-position disposition (including a failure if inaccessible) is in `output/atnf_receipts.json`. A Gaia non-match does **not** establish no optical/IR counterpart; no counterpart-free object is claimed.

**Still pending / not established:** independent pixel confirmation of the 14 patterns and calibration/leakage/morphology; RACS-low/mid/high Stokes I/V and linear polarisation; Legacy Surveys and unWISE forced optical/IR counterpart checks; a complete ATNF/pulsar-survey and literature review; eRASS1/Swift X-rays; and especially a periodicity search in public snapshot or 10-second visibility-based light curves. These DR1 epoch tables alone cannot resolve a 10-minute–hours pulse train. **No upper limit on new LPT incidence and no [REMARKABLE] result** is justified. Read the code and [referee report](output/referee.json) before interpreting individual IDs.

For a background contamination context, the DR1 discovery team classified **84 of 117 published strongly variable sources as AGN or radio stars** (44 AGN, 40 stars; 27 known pulsars, 6 other). The fraction **84/117 = 71.8%** describes that published sample, *not* a calibrated false-positive rate for these 149 sources. No random-sky negative-control normalization, image artefact rate or recovery efficiency has been measured, so a survey-wide false-positive estimate is **not available**.

## Input provenance and reproduction

| Published primary input | Verified source SHA-256 | Size |
| --- | --- | ---: |
| `VAST_Extragal_DR1_image_df.csv` | `70ce692b3c8d067d53e248b6265420a26880147bd7bc15327930880bd8395195` | 936,444 B |
| `VAST_Extragal_DR1_source_df.parquet` | `fff50cd07f2a092e50b822a27e2c4de1e229876deb394508e2345ccd8ca737f4` | 76,405,315 B |
| `VAST_Extragal_DR1_measurement_df.parquet` | `a3af873784b4b3d5a727d6d50890a48a5fd1bbc41ebe035106ce83d47cfb9ca6` | 442,056,111 B |

The exact *frozen image JSONL* SHA-256 is `e59daee09e646af6c89e722d752b8d8b1ff2c7b50c181f973d4e6bb3f5cac6ba`. The [original paper](https://doi.org/10.1017/pasa.2026.10171) and [authors' public example code](https://github.com/IrisdeRuiter/VAST_Extragalactic_DR1/tree/c1d014a51660ee31baa93cd3accf199b4260347f) describe the measurement definitions and downloads. This is public CC BY 4.0 research data. No bulk 520 MB source tables are vendored.

Run the single [push-triggered GitHub job](../../.github/workflows/askap-lpt-vast-dr1.yml), or replay with Python 3.12, NumPy, Pandas, PyArrow, Requests and Astropy:

```sh
python -m unittest discover -s results/askap-lpt-vast-dr1 -p 'test_*.py' -v
python results/askap-lpt-vast-dr1/search.py --cache /tmp/vast-dr1 --out results/askap-lpt-vast-dr1/output
python results/askap-lpt-vast-dr1/referee.py results/askap-lpt-vast-dr1/output
```

The workflow commits the compact field list, 149 source rows, 60 lead measurements and receipts in one batch; the generated `output/` files are **excluded from the trigger**. The independent review checks joins, thresholds, ordering, denominators and the documented one-row difference, *not* any claim of pulse periodicity. This is a **catalogue triage / incomplete physical validation**, not a discovery or a well-calibrated LPT population null.
