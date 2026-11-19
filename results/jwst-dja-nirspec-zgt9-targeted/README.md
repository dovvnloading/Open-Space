# DJA v4.4 targeted JWST/NIRSpec z>9 follow-up — 999 real spectra, no validated discovery

**Claim (2026-11-19 UTC): GitHub Actions successfully downloaded and fitted 999 targeted public JWST/NIRSpec spectra, but established **NO securely new z>12 galaxy**, **NO calibrated `P(low-z|D)`**, and **NO survey-wide null.** The 1D likelihoods are exploratory, 2D images are coadds, independent NIRCam and source-level checks remain undone, and a referee found strong false-positive high-z model ranks. This PR must **remain DRAFT and UNMERGED**. Do not mark `[REMARKABLE]` or edit `docs/REMARKABLE_FINDINGS.md`.

## What was actually downloaded and computed

**Observed [GitHub Actions run #44454883069](https://github.com/dovvnloading/Open-Space/actions/runs/44454883069) (2026-11-19 03:56–04:02 UTC, successful):**

| Stage | Observed evidence |
| --- | --- |
| Publisher release/version | [DJA v4.4](https://zenodo.org/records/15472354), 80,367 *spectrum records*, not independent galaxies |
| Summary table | **133,004,545 bytes** Zenodo GET 200; publisher MD5 matched `721c58c446f4a4eb1570aef02c4ab47b`, observed SHA-256 `728afcdd5ec4d164b9d816ea903e648b1d4937ae52c8a271056c9119cabf1ccf` |
| Selection | **999/80,367** catalog rows, any DJA `z_best` or `z_phot` >9 or six published z>12 baseline controls; all six controls found. **79,368 other spectra not fitted** |
| Targeted downloads | **999/999** individual DJA S3 `.spec.fits` retrieved and FITS SPEC1D-validated; **0 failures**, total 2,038,757,760 downloaded bytes; per-URL SHA-256 receipts |
| Model fits | **999/999** uncalibrated high-z Lyman-break/wing vs dusty Balmer/line/stellar-artifact 1D proxy fits; **0 numerical failures** |
| z>12 review figures | **64** SCI/WHT coadded-2D/1D figures, not 64 unique galaxies; every displayed file has SCI, but **not** separate native visit/DQ validation |
| Independent scientific validation | **0** NIRCam crosschecks, **0** full per-source priority/alternative-model posterior validations, **0** securely new discoveries |

Detailed committed outputs: [run and publisher receipt](outputs/run_receipt.json) · [selection receipt](outputs/selection_receipt.json) · [complete compact 80,367-row denominator](outputs/compact_catalog_index.csv) · [selected list](outputs/selected_spectra.csv) · [999 download/hash receipts](outputs/download_receipts.jsonl) · [999 unvetted fits](outputs/fit_results.csv) · [64 review rows](outputs/zgt12_literature_review_queue.json) · [coadded 2D figures](outputs/figures/) · [runner connectivity](outputs/runner_host_probe.json). Raw FITS (~2.04 GB) and 133 MB original CSV are not checked into Git; the full audit artifact is available in [the Actions run](https://github.com/dovvnloading/Open-Space/actions/runs/44454883069).

**Environment contrast:** the prior *local* container had six DNS failures; the **GitHub Actions runner resolved Zenodo/S3 and obtained all 999 selected products**. MAST API HEAD returned HTTP 500 while the MAST archive homepage returned 200; this job retrieved DJA, not native MAST detector exposures.

## Independent-referee outcome and priority

The [2026-11-19 referee report](REFEREE.md) and [64/64 per-row literature/identity ledger](literature_referee_2026-11-19.json) group the review list into **27 sky positions** within 0.25 arcsec: **10 already published z>12 systems** (37 spectra), **5 published lower-z objects** (5 spectra), **2 catalog lower-z cases** (2 spectra), **one wrongly matched neighboring positive control** (10 spectra), and **9 unresolved source groups** (10 spectra). A finite web literature pass is not an exhaustive ADS/arXiv priority search; none of the unresolved objects is established as new.

The severe **low-redshift interloper controls** invalidate claims based on the uncalibrated proxy scores alone: RUBIES-UDS-144195 (published z≈3.36, strong Balmer break) scores nearly 1 for z>12; JADES-GN-11508 (published z=9.93) also scores 1; JADES-GS-20021387 (published z=11.06) scores 0.992. An automated 0.5-arcsec positive-control match also confused JADES MSA **183349** (DJA z=3.47) with **GS-z14-0 MSA 183348**, ~0.37 arcsec away. The **next** target selector tightens the positional window to 0.2 arcsec and adds a synthetic regression test; the **original run receipts remain untouched**.

Published comparison sources: [MNRAS 2026 Table 1](https://academic.oup.com/mnras/article/548/3/stag701/8653929) · [Tang et al. 2026 z=9–14 sample](https://doi.org/10.3847/1538-4357/ae4edb) · [RUBIES low-z Balmer false positive](https://doi.org/10.1051/0004-6361/202555816). Previously published z>12 aliases in these data include MoM-z14, GS-z14-0/1, GS-z13-0, JADES-GS-z13-1-LA, UNCOVER-13077, GHZ2/GLASS-z12, GS-z12-0, UNCOVER-38766 and CAPERS-EGS-65480. These are **not new findings**.

## Method and reproduction

This is a **targeted** continuation of [merged PR #195](https://github.com/dovvnloading/Open-Space/pull/195), whose frozen v4.4 scope was 80,367 combined spectrum records but fit zero pixels. The new code selects the publisher-MD5-verified *catalog first*, then fetches individual FITS via DJA S3 HTTPS; it **does not** need the nine combined Zenodo spectral FITS shards. Selection cannot detect objects that the DJA catalog assigns both `z_best<=9` and `z_phot<=9` unless they match published controls; it is **not** a complete independent fit of all 80,367 spectra. The 1D model reuses PR #195's illustrative templates, not calibrated physically marginalized interloper evidence.

To reproduce on a runner with external network, see [WORKFLOW.md](WORKFLOW.md) and `run_dispatch.py`. On this branch the existing [JWST workflow](../../.github/workflows/jwst-dja-nirspec-preflight.yml) supports manual dispatch **and future automatic pushes that touch this result folder**, excluding `outputs/**` to avoid self-trigger loops. Future attempts write `outputs/runs/<run_id>/` without overwriting the first successful run. The initial manual run remains frozen in `outputs/`. The only validated numerical method tests are synthetic fixtures; no native per-exposure detector/slit checks or NIRCam filter-convolved SEDs have been performed.

**Required before readiness/merge:** independently validate native 2D integrations and DQ/contamination, measure NIRCam fluxes/limits and convolve both high-/low-z models, calibrate stellar/Balmer/line-confusion posterior odds against published low-z controls and injected false breaks, complete alias/coordinate literature checks for all nine unresolved groups (ADS, arXiv, DJA/JADES/CEERS/UNCOVER/MoM/GLASS). An actual `z>14` galaxy requires coherent independent spectroscopic and photometric evidence, not a model rank.


## 2026-11-19 lead follow-up — independent nine-group referee

**Updated evidence, not a discovery:** GitHub Actions [#44459389472](https://github.com/dovvnloading/Open-Space/actions/runs/44459389472) and [#44460415546](https://github.com/dovvnloading/Open-Space/actions/runs/44460415546) downloaded/reviewed the 10 spectra for the nine initially unresolved sky groups, checked their coadd SLITS metadata, matched JADES DR4 catalogue identities, and ran bounded 1D surrogate interloper injections. The new [full referee](REFEREE_FOLLOWUP_2026-11-19.md) and [machine-readable nine-group ledger](FOLLOWUP_NINE_2026-11-19.json) supersede the previous nine-way *unresolved* tally as follows:

| Sky group | DJA program/source | Evidence-bounded disposition | Discriminating observation |
| --- | --- | --- | --- |
| 1 | `3073_50252` | Unresolved; no independently validated redshift | No independent physical redshift evidence |
| 4 | `6368_211132` | Unresolved; no independently validated redshift | No independent physical redshift evidence |
| 6 | `6368_151513` | Unresolved; no independently validated redshift | No independent physical redshift evidence |
| 11 | `1181_19411` | JADES DR4 target, quality E; no measured redshift | DR4 `goods-n-mediumjwst_19411`, flag E |
| 16 | `4233_22104` | Unresolved; no independently validated redshift | No independent physical redshift evidence |
| 18 | `4233_7420` | Unresolved; no independently validated redshift | No independent physical redshift evidence |
| 22 | `3215_20095279` | JADES DR4 target, quality E; no measured redshift | DR4 `goods-s-ultradeep_20095279`, flag E |
| 24 | `1286_176143` | JADES DR4 z=3.60656, quality D (tentative, not secure) | DR4 z=3.60656, flag D; 3/20 dusty interlopers pass false-high gate |
| 26 | `1180_b37` | Empty-background shutter B37; not a galaxy | 3/3 SLITS records refer to `background_37`; source RA=Dec=0 |

**Summary:** **1 sky-background artefact** (`1180_b37`; all raw-source SLITS rows marked `background_37`, RA=Dec=0), **1 tentative lower-z** (`1286_176143`, JADES DR4 z=3.60656 flag D), **2 JADES DR4 catalogued survey targets with redshift flag E** (`1181_19411`, `3215_20095279`), and **5 still unresolved**. Published survey target does **not** imply published measured high-z. All native-integration DQ checks, independent NIRCam flux/upper-limit measurements, complete filter-convolved alternatives, calibrated physical redshift posteriors and exhaustive ADS/arXiv/survey alias priority checks **remain missing**. The 20-per-family 1D injection diagnostics cannot calibrate physical redshift odds; 3/20 dusty interlopers already fool the near-1-scored `1286_176143` fit. **Zero validated new z>12 galaxies. PR #197 remains DRAFT, UNMERGED, with no merger SHA. No `[REMARKABLE]` label.**
