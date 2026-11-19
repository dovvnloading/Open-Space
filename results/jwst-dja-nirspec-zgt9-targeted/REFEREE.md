# Independent-referee audit — DJA v4.4 targeted z>9 screen

**Referee decision (2026-11-19 UTC): CHANGES REQUIRED; KEEP PR #197 DRAFT / UNMERGED.** This is an adversarial second analysis of the independently obtained GitHub Actions runner output, not formal independent peer review. We could audit committed tabular receipts and literature, but the local environment cannot open the binary FITS or rendered PNG figures to visually inspect their pixels. Do **not** interpret this as a 2D-exposure or photometric validation.

## Observed run and provenance

- **Run:** [GitHub Actions 44454883069](https://github.com/dovvnloading/Open-Space/actions/runs/44454883069) finished **successfully**. GitHub-hosted Actions has DNS and HTTPS for the relevant archive hosts, unlike the prior local container. [Runner probe](outputs/runner_host_probe.json): Zenodo/S3 HTTPS 200, DJA documentation 200, MAST archive 200; MAST API HEAD 500 (an unsupported HEAD response does not establish that API product retrieval failed).
- **Pinned table:** 133,004,545 bytes, Zenodo HTTP 200, publisher MD5 `721c58c446f4a4eb1570aef02c4ab47b` matched; observed SHA-256 `728afcdd5ec4d164b9d816ea903e648b1d4937ae52c8a271056c9119cabf1ccf`.
- **Denominator/selection:** 80,367 v4.4 catalog spectra; only `z_best`/`z_phot` are present among the selected high-z hint fields. 999 spectra selected. **79,368 catalog rows were NOT independently 1D-fit**; a high-z source incorrectly assigned low/no catalog redshift would be omitted unless a baseline control. No survey-wide completeness/null is established.
- **Individual spectra:** 999 S3 downloads, **0 failed**, totaling 2,038,757,760 observed bytes, with per-URL SHA-256 and FITS SPEC1D checks in [download receipts](outputs/download_receipts.jsonl). The runner made 999 unvetted 1D fits and 64 SCI/WHT 2D PNG figures (all 64 SCI present). Binary FITS stayed in ephemeral runner storage and are **not** retained as an independent replayable artifact; Git SHAs and output hashes alone do not prove model validity.
- **Source files:** [run receipt](outputs/run_receipt.json), [selection receipt](outputs/selection_receipt.json), [999-row fit table](outputs/fit_results.csv), [80,367-row compact denominator](outputs/compact_catalog_index.csv), [64-row queue](outputs/zgt12_literature_review_queue.json), and [figure directory](outputs/figures/). Original large uncompressed source files remain outside Git. Audit artifact `dja-v44-pr197-observed-audit` from run #44454883069 is a separate downloadable Actions asset.

## Independent priority/identity check

The 64 z>12 *review rows* are **not 64 distinct candidate galaxies**. Grouping the published coordinates within 0.25 arcsec produces **27 sky groups**; the complete per-row 64/64 audit, identity notes, published references and pending queries are in [`literature_referee_2026-11-19.json`](literature_referee_2026-11-19.json).

| Distinct sky groups | Rows | Independent bounded literature/identity disposition |
| --- | ---: | --- |
| 10 published z>12 galaxies | 37 | **Already known:** MoM-z14, GS-z14-0, GS-z14-1, JADES-GS-z13-1-LA, GS-z13-0, UNCOVER-13077, **GHZ2/GLASS-z12, GS-z12-0, UNCOVER-38766, CAPERS-EGS-65480**; the last four were not among the PR #195 six frozen controls |
| 5 published z<=12 objects | 5 | **Already known lower-z:** GLASS-10003 (z=9.3127), JADES-GS 10058975 (z=9.438), RUBIES-UDS-144195 (z~3.36), JADES-GN-11508 (z=9.93), JADES-GS-20021387 (z=11.06) |
| 2 catalog lower-z cases | 2 | DJA z_best 1.985/3.682; no independent z>12 confirmation |
| 1 false positive-control identity cluster | 10 | JADES MSA **183349**, DJA `z_best=3.4737`, was incorrectly labelled published **GS-z14-0 (183348)** on 0.5-arcsec positional match. These are distinct IDs about **0.37 arcsec** apart. The next selector uses 0.2 arcsec and adds a regression test |
| 9 unresolved sky groups | 10 | **Not known to be new.** Bounded public-web searches did not establish a unique published identity; these are an unvetted priority queue, requiring native 2D/NIRCam and exhaustive ADS/arXiv literature work |

**Falsifying counterexamples for the apparent high-z model confidence:** the published **RUBIES-UDS-144195** at z≈3.36 receives `p_z_gt12_proxy≈1`, though the published paper specifically notes that its Balmer break mimics a photometric z~14 Lyman break; **JADES-GN-11508** at z=9.93 receives proxy 1 and z≈17.65, **JADES-GS-20021387** at z=11.06 receives proxy 0.992 and z≈14.81, and **GLASS-10003** at z=9.3127 receives proxy 0.467 at a median z≈13.05. Model values slightly above 1 in the CSV are floating-point artifacts, not physical probabilities. These show that the simplistic stellar, line and Balmer model **is not calibrated for priority, false-positive control or credible high-z odds**. An isolated grating/proxy peak is not a robust galaxy redshift.

### Published identity anchors

- [MNRAS 548(3), stag701 (2026), Table 1](https://academic.oup.com/mnras/article/548/3/stag701/8653929): ten exact source/sky matches including published z>12 controls, GHZ2, GS-z12-0, UNCOVER-38766 and CAPERS-EGS-65480.
- [Tang et al., *The JWST Spectroscopic Properties of Galaxies at z=9–14* (2026)](https://doi.org/10.3847/1538-4357/ae4edb), Table B1: JADES-GN-11508 and JADES-GS-20021387 are known at z<12.
- [Hviding et al., A&A 702, A57 (2025)](https://doi.org/10.1051/0004-6361/202555816): RUBIES-UDS-144195 has a broad Balmer break at z≈3.36, directly refuting the proxy high-z classification.
- [GLASS-10003/Gz9p3 spectroscopy](https://escholarship.org/content/qt80h7j3rn/qt80h7j3rn_noSplash_e869eaa1a99b4248f7972397c5ff9da6.pdf): z=9.3127, not z>12.
- [Cameron et al. 2023 JADES NIRSpec target table](https://uhra.herts.ac.uk/id/eprint/10807/2/aa46107-23.pdf): JADES-GS+53.11243-27.77461 = 10058975 at z=9.438.

## Mandatory scientific failures

1. **Not an exhaustive z>9 1D scan:** 79,368 rows were never fitted. The selected 999 were chosen by catalog/photometric hints, so no completeness/survey-wide null or full-source denominator.
2. **No calibrated competing-model odds:** the PR #195 five-family diagonal-noise Lyman-step/damping-wing proxy is not a physically complete dusty Balmer, emission-line and stellar atmosphere marginal evidence. The published low-z controls above demonstrate serious high-z false positives, so no scientific `P(low-z|D)` or robust redshift interval follows.
3. **No independent native 2D validation:** the 64 displayed 2D coadds lack separate detector integrations, DQ/saturation/contamination and shutter/neighbour checks. Also this referee could inspect PNG/FITS existence and hashes through GitHub metadata but **not visually audit binary pixels** in the text-only connected environment.
4. **No NIRCam photometry/SED test:** zero independently checked NIRCam fluxes, limits, aperture/identity matches or filter-convolved comparisons.
5. **Literature priority incomplete:** the bounded Table-1/Table-B1/web check identifies known objects and disqualifies obvious false positives, but every unresolved target still requires independent ADS, arXiv, DJA, JADES, CEERS, UNCOVER, MoM and GLASS checks by ID and coordinates before any "new" label. Absence of a search result is NOT evidence of novelty.
6. **Coverage version:** DJA v4.4 is a 2025 archive snapshot; objects from the 2026 literature may be absent. The 64-row queue is the pipeline shortlist, not every z>12 galaxy in the Universe.

## Release verdict and push policy

The new push job will be **branch/path filtered** to `results/jwst-dja-nirspec-zgt9-targeted/**`, excluding generated `outputs/**` to prevent loops. Each future run writes `outputs/runs/<run_id>/` so receipts do not overwrite run #44454883069. This configuration is staged in the same PR and not a separate claim of an additional run.

**Approve only an internally useful exploratory workflow+receipts, not a scientific discovery. Under the user's evidence-gated merge instructions, do NOT mark ready or merge #197.** No new secure z>12 object, no survey-wide null, no `[REMARKABLE]` or `docs/REMARKABLE_FINDINGS.md` entry is justified.


## 2026-11-19 lead follow-up (supersedes the earlier nine unresolved tally)

The nine positions were re-examined in two further successful runner passes. Read [the independent full updated decision](REFEREE_FOLLOWUP_2026-11-19.md) and the [per-group machine-readable dispositions](FOLLOWUP_NINE_2026-11-19.json). **One `1180_b37` is empty-sky background, one `1286_176143` has a tentative quality-D JADES DR4 z=3.60656, two other JADES positions have quality-E no-redshift matches, and five remain unresolved.** None is validated as new z>12. Native detector DQ, complete shutter contamination, NIRCam flux/upper limits, physical redshift posteriors and exhaustive literature priority still fail the scientific gate. **CHANGES REQUIRED; KEEP DRAFT / UNMERGED; no merge SHA.**
