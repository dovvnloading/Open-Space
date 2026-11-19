# Skeptical referee audit — targeted DJA z>9 continuation

**Date:** 2026-11-19 UTC. **Decision: CHANGES REQUIRED / KEEP DRAFT AND DO NOT MERGE.** This is an explicitly adversarial second-pass *author-side* audit, **not** an independent human peer review. It examines only committed source code, host receipts, and synthetic local tests; real astrophysical source data are absent.

## Findings

| Gate | Disposition and test evidence |
| --- | --- |
| Host access, actual archive retrieval | **FAIL/BLOCKED.** Six endpoint-level HEAD probes: zero HTTP responses; Zenodo, DJA path-style S3, virtual-host S3, DJA docs, MAST API and MAST archive all failed DNS resolution or timed out. See `host_probe_2026-11-19.json`. Tests ran locally, not on GitHub Actions. |
| Dataset scope | **DOCUMENTARY ONLY.** PR #195 cites 80,367 DJA v4.4 source-*spectra*. In this continuation no publisher MD5 was recalculated because no summary bytes arrived. No denominator was independently established. |
| Selection completeness | **NOT DEMONSTRATED.** `select_targets.py` rejects wrong publisher MD5 and unexpected row count, writes an all-row audit and selects any available named redshift hint greater than 9 or a published-control ID/position. Its rule and output ledger passed local synthetic tests. However the actual redshift-column schema and any high-z posterior tail are unverified; entries without hints could harbor a >9 solution. This is **not** a complete posterior-based search. |
| Positive controls | **NOT RECOVERED.** Six literature examples in PR #195's baseline are known published objects; matching by filename suffix or 0.5-arcsec coordinate was tested only against made-up rows. A missing control must be reported, not counted as a non-detection. v4.4 is a 2025 snapshot and cannot be assumed to include every 2026 publication. |
| Targeted retrieval | **NOT EXECUTED.** `download_targets.py` only *plans* hashes-checked selections and individual S3 `.spec.fits` URLs. Local synthetic test passed traversal rejection and duplicate-URL de-duplication. No remote FITS was retrieved or parsed; raw bytes, instruments and 2D HDUs are unverified. |
| 1D fit / low-z interlopers | **ZERO.** The earlier PR #195 method is an uncalibrated proxy and consumes combined FITS shards, not the per-source files in this continuation. A calibrated physical Lyman break/damping wing, Balmer/dust, stellar atmosphere and line-confusion comparison has **not** run. No `P(low-z|D)`, redshift interval or sensitivity is defensible. |
| 2D independent exposure review | **ZERO.** An individual DJA `.spec.fits` may include a SCI 2D extension, but that is not necessarily an independently reduced visit. No native detector DQ, MSA shutters, separate integration checks, contamination model, or slit-trace overlay. |
| NIRCam / literature priority | **ZERO object-specific validations.** No flux-limit catalog, NIRCam passband-convolved SED, ADS/arXiv object search, crossmatch to JADES/CEERS/UNCOVER/MoM/GLASS, or proof of originality. The prior literature baseline is *not* an independent per-candidate search. |
| Remarkable and release decision | **PASS the no-overclaim gate.** No `[REMARKABLE]`, no modification to `docs/REMARKABLE_FINDINGS.md`, no discovery or null claim. Leave draft and unmerged. |

## Local executable checks

- `python select_targets.py --self-test`: **PASS** (synthetic row flags; positive control by exact ID and spherical separation).
- `python download_targets.py --self-test`: **PASS** (per-source URL allowlist and traversal rejection).
- Small synthetic 5-row gzip-catalog end-to-end selector: **PASS**, 3 selected spectrum rows and 6 audit lines including header, *with intentionally overridden test-only checksum/count*, which do **not** demonstrate the published 80,367-row inputs.
- Small synthetic two-row downloader `--plan-only` test: **PASS**, one unique URL; no network operation and no FITS bytes.

These are tests of program logic, not 1D fitting, 2D detection or astrophysical correctness.

## Mandatory before ready/merge

Obtain and validate the **actual v4.4 summary bytes** (published MD5 and locally computed SHA-256). Inspect real columns and the full 80,367-row denominator, then publish the auditable z>9 set and demonstrate recovery of all known controls present in release. Independently fit each downloaded source with calibrated high-/low-z/stellar alternatives, inspect genuinely separate native 2D exposures, check NIRCam, and perform per-object up-to-date literature priority searches. Reject any sample-wide null until completeness and selection effects are measured. If DNS access remains blocked, **do not merge**.

**Referee verdict:** the repository has a reproducible blockage receipt and a cautious future target-selection/retrieval method, **not scientific results**.

## Follow-up workflow code review (before live run)

One gated `workflow_dispatch` job is now staged on the PR branch by extending the existing default-branch dispatchable preflight workflow. The new code will download the publisher-MD5-verified summary, select available z>9 hints/controls, stream selected individual DJA FITS, evaluate the previously documented **uncalibrated** five-family 1D templates and attempt to visualize `SCI`/`WHT` 2D coadds. This is a **source-code review only**: the new self-test and live data path have not been executed here and no runner FITS receipt has yet been inspected. Automatic `UNREVIEWED_NOT_NEW` labels deliberately stop short of a literature priority verdict; there is no independent NIRCam photometry, native exposure/DQ inspection, or selection-function calibration. If a runner GET fails or its artifact has missing objects, reject null/discovery interpretation. **Verdict remains KEEP DRAFT / DO NOT MERGE** until actual outputs are available and reviewed.
