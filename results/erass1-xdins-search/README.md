# eRASS1 DR1 XDINS search — 1,330 soft sources, no independently validated new XDINS

**Result (19 November 2026 UTC; bounded catalogue screen, not a calibrated population null):** The [official eROSITA-DE western-hemisphere eRASS1 DR1 main v1.2 catalogue](https://erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/Catalogues_dr1/) contains **930,203 sources**; the hard v1.0 contains **5,466**. A pre-registered-for-this-run soft, compact, reliable-source screen retained **1,330**. The independent [Salvato et al. eRASS1 counterpart releases](https://erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/Catalogues_dr1/) assign a high-probability (NWAY p_any ≥ 0.8) association in Gaia DR3, Legacy Surveys DR10 or CatWISE2020 to **1,056**; **274** lack an association *at that specific probability threshold*, **not** confirmed multiwavelength nondetections. Among the top **60** by X-ray counts × softness, **25** have an association and **35** do not. **Zero newly identified XDINS candidates passed the required optical-limit, spectrum, variability, archival and novelty checks, because those checks are incomplete. No new XDINS is claimed, nor is a calibrated null.**

The source list contains **known** soft neutron stars: all **five Magnificent Seven (M7) XDINS in the released western footprint** are recovered in the frozen source selection; the other two M7 lie outside the public DR1 hemisphere. The previously reported eRASSU J065715.3+260428 and the [already confirmed 2026 XDINS](https://doi.org/10.1051/0004-6361/202659673) eRASSU J131716.9−402647 are also recovered. These are positive controls, **not discoveries**. The [control separations and individual decisions](output/controls.json) are exact outputs; epoch/proper-motion differences remain.

## Frozen source selection and actual denominators

The [frozen protocol](output/protocol.json) and [runner](run.py) were committed before any main catalogue input was processed on the network-enabled GitHub Actions runner. The screen requires a main-band flux ≥ 10⁻¹³ erg cm⁻² s⁻¹, broad-band detection likelihood ≥12 and ≥20 net counts, position error ≤10 arcsec, point-like extent (EXT=0; EXT_LIKE≤6), ≥5 band-P1 net counts and (P1+P2)/(P1+P2+P3) ≥0.90. It excludes catalogue optical-loading/bright-point/overdensity flags and sources with strong hard-catalogue associations. A >5× brighter eRASS neighbour within 60 arcsec is separately flagged. These are screening thresholds, not an exposure/absorption-corrected thermal spectral model; eRASS flux assumes a power-law conversion and is **not** a blackbody fit.

| Stage | Count |
|---|---:|
| Whole DR1 main catalogue | 930,203 |
| Detection likelihood, broad-band counts and flux | 75,645 |
| Compact/positional cuts | 61,624 |
| Softness and P1-count cuts | 1,659 |
| Catalogue contamination flags removed | 1,330 |
| Strong hard associations excluded | 1,330 |
| Brighter-neighbour screen retained | 1,319 |
| Ranked, multiwavelength-follow-up subset | 60 |
| High-probability Gaia/Legacy/CatWISE association, among 1,330 | 1,056 |
| No high-probability association, among 1,330 | 274 |
| High-probability association, among 60 | 25 |
| No high-probability association, among 60 | 35 |
| Validated newly discovered XDINS | **0 (validation incomplete)** |

Counts and full 1,330-row source-level evidence: [candidate_screen.json](output/candidate_screen.json) · [screened_main.csv](output/screened_main.csv) · [official_counterparts.csv](output/official_counterparts.csv) · [official_top60.csv](output/official_top60.csv). **Do not interpret 274 or 35 as new neutron-star candidates.**

## Public data hashes, provenance and receipts

The [original main/hard download receipts](output/receipts.json) contain timestamps, URLs, redirects, ETags, Last-Modified, byte lengths, catalogue row counts and both compressed and extracted SHA-256 checksums. Raw archives are intentionally not committed (main FITS exceeds 1 GB).

| Frozen upstream archive | Downloaded bytes | SHA-256 of archive |
|---|---:|---|
| eRASS1 Main **v1.2** FITS tar.gz | 643,072,855 | `a4af6173c93bd4367d8a37465bdf5061b6092f02b8f522fc68074afff0e59131` |
| eRASS1 Hard **v1.0** FITS tar.gz | 2,057,573 | `ef1d2dbcff2c193f1434d58cb6dc97e6370a1e8af24752a379b265313cc48837` |
| Salvato eRASS1 **LS10** COLFITS tgz | 490,390,929 | `368eb4994328725a63222d330f8315d144cf4778fe68876ef038c6c4087aa29c` |
| Salvato eRASS1 **Gaia DR3** COLFITS tgz | 398,248,982 | `d1d85833783c4f24b9769cbd0d3e6909315bd941e20934561b59fe868f920754` |
| Salvato eRASS1 **CatWISE2020** COLFITS tgz | 350,361,463 | `f24f2d543ffa1bb4ceeb78d9e49c2bffcf68235f8150b4e15c2164210c314af3` |

The [official counterpart receipts](output/official_receipts.json) also include the **extracted FITS hashes**, row counts (LS10 966,194; Gaia 959,261; CatWISE 1,001,527) and missing-field list. The analysis uses an **exact eROSITA UID join**, not a probabilistic positional cone-match masquerading as an identification. The published NWAY p_any score is a probability, not a guarantee. We call ≥0.8 a *high-probability association* and retain weaker/alternative matches in the source table.

## Evidence completed versus pending

**Completed and deposited:** official full 930,203-row main / 5,466-row hard downloads and row-count checks; SHA-256 of archives and extracted FITS; fixed X-ray gate counts for *all main rows*; 1,330-row and top-60 lists; brighter-neighbour flags; exact-UID independent join to all three June-2026 official optical/IR counterpart catalogues; LS photometric-redshift and AGN-colour/coronal tags where a valid association row exists; source-level NWAY likelihood; five-of-five in-footprint M7 and two published eROSITA INS positive controls; 11 independent output consistency checks ([referee.json](output/referee.json)) and a separate [adversarial scientific review](REFEREE.md). The original 10 CDS XMatch requests **all timed out at TCP connect**. Their initially empty results represent **failed queries, not non-detections**. [Exact request/error receipts](output/receipts.json) remain unmodified; the direct official releases were obtained in a subsequent batch.

**Unresolved:** actual *no-counterpart* classifications and coverage/depth limits for Gaia, Legacy DR10, Pan-STARRS, DES, unWISE and 2MASS at source positions; independent PS1/DES/2MASS/unWISE catalogue queries (CDS service unreachable from runner); proper-motion-aware matching; robust (f_X/f_{opt}>10^4) constraints; radio/IR colour, catalog star/AGN and cluster rejection where the NWAY association is weak; systematic X-ray light-curve **non-variability**, eRASS event-level background checks, source images/spurious-confusion testing, blackbody (kT=40–120) eV plus absorption fits, XMM-Newton/Chandra/Swift archive detections and **pulsation timing**; current ATNF radio catalogue; the full 33-object [Kurpas et al. 2024](https://doi.org/10.1051/0004-6361/202450248) novelty list and 2025–26 follow-up literature; injection/recovery/contamination completeness. Kurpas' full VizieR retrieval failed and is logged; the two explicitly named published eROSITA sources are accounted for but **not** a comprehensive novelty audit.

**Why the optical-flux requirement remains open:** only **18/60** top-ranked sources have a reported LS10 5σ g depth on the matched counterpart row; no source has a calibrated >10⁴ X-ray/optical lower bound from the current data. The notebook's (g=V) Maccacaro formula is expressly an illustrative *magnitude required* (or uncertain proxy) and cannot supply a physical limit in place of actual depth, bandpass and Galactic extinction. A missing or low-probability match can occur without survey coverage.

## Replay and scientific status

The **single** [.github/workflows/erass1-xdins.yml](../../.github/workflows/erass1-xdins.yml) job runs on push to `research/erass1-xdins-20261119` only when analysis scripts/dependencies under this result folder change. It does not run for all pull requests or on generated `output/**` commits. It installs pinned [requirements](requirements.txt), obtains upstream data and computes outputs in one job, commits artifacts in a batch and executes [referee.py](referee.py). First-run inputs were cached as committed small, SHA-receipted outputs so the separate **official_followup.py** job run reuses the frozen X-ray selection without an unnecessary 1 GB re-download. To independently replay the full source screen, remove `output/screened_main.csv` and run `python run.py`; to refresh the official tables, `python official_followup.py`; then `python referee.py`.

The published search by **[Kurpas et al. (2024)](https://doi.org/10.1051/0004-6361/202450248)** already reported **33** candidates and estimated only 30–50% completeness. A validated *new* object would need exclusion of those candidates and subsequent follow-up, not just absence in catalogues. The 2026-confirmed eRASSU J131716.9−402647 is an instructive control, **not new today**. This investigation is a bounded, reproducible catalogue/association audit only. Under [docs/REMARKABLE_FINDINGS.md](../../docs/REMARKABLE_FINDINGS.md), **no [REMARKABLE] tag or new log entry is justified**.
