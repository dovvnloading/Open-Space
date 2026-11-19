# Targeted DJA v4.4 z>9 NIRSpec screen — blocked data-access continuation

**Claim (2026-11-19 UTC): No astronomical spectrum was downloaded, fitted, or validated in this continuation. There is NO new z>12 discovery, no calibrated interloper probability, and no survey-wide null result.** This is a blocked continuation of [merged PR #195](https://github.com/dovvnloading/Open-Space/pull/195), not a replacement of its frozen 80,367-spectrum scope.

## GitHub Actions data attempt (new, not yet executed)

A **single manually dispatched** GitHub-hosted runner is staged in the existing JWST preflight workflow, branch-gated to PR #197. See [WORKFLOW.md](WORKFLOW.md) for the exact Run workflow steps, pinned Zenodo summary and individual DJA S3 FITS strategy, downloaded SHA-256 receipts, uncalibrated five-family 1D fit, attempted 2D cutouts and candidate literature-review queue. The job will commit observed output receipts in **one** push. **This source-code addition is not evidence of any successful download, astrophysical fit, literature priority check or validated galaxy.** Until `outputs/run_receipt.json` exists with real bytes, the only observed network result remains the local DNS failure described below. The job itself does not establish independent NIRCam photometry, separate native 2D exposures or novelty; its `z>12` queue requires a separate skeptical human-level review.

## Why the targeted analysis was locally blocked

At 2026-11-19 03:23:24 UTC, the local compute container failed DNS/HTTPS HEAD tests of the exact [Zenodo v4.4](https://zenodo.org/records/15472354) summary CSV, the [DJA S3 extraction host](https://dawn-cph.github.io/dja/general/api_summary/) (path-style and virtual-host style), the DJA documentation host, and two MAST hosts. See [observed machine-readable probe](host_probe_2026-11-19.json) and [reproducible probe script](host_probe.py). All six had HTTP code `000`; curl exited 6 (could not resolve host) or 28 (DNS timeout). The old PR already recorded an independent failed retrieval. Web-indexed descriptions of the archive are not data bytes and **cannot** substitute for FITS and hash receipts.

## Fixed targeted selection and validation protocol

1. Download **only** the publisher-MD5-verified DJA v4.4 *summary* `dja_msaexp_emission_lines_v4.4.csv.gz` (133 MB; **not** the nine 8.4 GB combined spectrum shards). Freeze exactly 80,367 source-spectrum rows, sorted in original catalog order and with observed SHA-256.
2. Preserve **every row** with any finite high-redshift hint above 9 (DJA `z_best`, `z_phot`, available `zfit`/redshift interval bounds, or existing independent quick-fit upper support if actually available). Do not apply a visual-grade/SNR cut at selection. Add all published z>12 positive controls from [PR #195's literature baseline](../jwst-dja-nirspec-highz-screen/literature_baseline.json), matching against exact row identifiers *or* astrometric coordinates, with unmatched controls reported rather than invented as matches. Retain the full selected-row and exclusion audit.
3. Obtain **individual** DJA `https://s3.amazonaws.com/msaexp-nirspec/extractions/{root}/{file}` spectra for only those selected rows, in checked path form. Use publisher archive identifiers and SHA-256 of each actual file, with no wildcard download. Do not treat a catalog row as an independent galaxy.
4. Independently fit actual 1D and 2D data against Lyman-break plus varying IGM damping wings, redshift grid z>9, dusty Balmer/4000 Å galaxies, strong-line confusion, brown-dwarf spectra and artefact/noise models. Vary priors, model assumptions, extraction apertures and data-quality masks. The simplified `screen.py` in PR #195 produces **uncalibrated proxy scores**; these are not physical posterior odds.
5. Inspect native SCI/ERR/DQ 2D pixels and disjoint MAST/DJA exposures for break consistency, saturation/bad pixels, traces and neighboring MSA shutters; verify NIRCam blue-band limits/red-band detections using measured filter throughput and the same SED templates. Search ADS/arXiv and JADES/CEERS/UNCOVER/MoM/GLASS + DJA, using coordinates and aliases, for **every** surviving z>12 solution. Re-reported controls are **known**, not new discoveries.

## Observed vs required

| Gate | Observed status |
| --- | --- |
| DJA summary downloaded and publisher-MD5 checked | **BLOCKED (DNS)** |
| Full 80,367-row summary screened / all DJA z>9 targets selected | **NOT DONE** |
| Targeted individual FITS downloaded and SHA-256 verified | **0** |
| Physical 1D fits / interloper comparisons | **0** |
| Independent 2D pixel inspection, NIRCam checks and object-specific literature priority | **0** |
| Claim/remarkable label | **NONE** |

The [published positive controls](../jwst-dja-nirspec-highz-screen/literature_baseline.json) include MoM-z14, JADES-GS-z14-0 and other *already known* z>12 objects; their existence **does not** establish an independent recovery here. Treat this branch as **draft/unmerged** while host retrieval is blocked. Never add `[REMARKABLE]` or alter `docs/REMARKABLE_FINDINGS.md` without independent discovery-level evidence.

## Local commands (not Actions)

The following tests ran **in the compute container** using only synthetic rows (no remote spectra): `select_targets.py --self-test` **PASS**; `download_targets.py --self-test` **PASS**; 5-row selector fixture **PASS** (3 selected); 2-row plan fixture **PASS** (1 unique URL). These do not test real source coverage or the physical fitting algorithm.

When DNS is restored, use the **existing pinned PR #195 downloader** for **catalog-only** retrieval; this deliberately skips every 8.4-GB combined FITS shard. Do not run the old full-shard `screen.py` on per-object `.spec.fits`.

```sh
python results/jwst-dja-nirspec-zgt9-targeted/host_probe.py > /data/dja-v44/host_probe.json
python results/jwst-dja-nirspec-zgt9-targeted/select_targets.py --self-test
python results/jwst-dja-nirspec-zgt9-targeted/download_targets.py --self-test
python results/jwst-dja-nirspec-highz-screen/retrieve.py download --catalog-only --cache /data/dja-v44
python results/jwst-dja-nirspec-zgt9-targeted/select_targets.py \
  --catalog /data/dja-v44/dja_msaexp_emission_lines_v4.4.csv.gz \
  --out /data/dja-v44/selection
python results/jwst-dja-nirspec-zgt9-targeted/download_targets.py \
  --selected /data/dja-v44/selection/selected_spectra.csv \
  --selection-receipt /data/dja-v44/selection/selection_receipt.json \
  --out /data/dja-v44/targeted
```

Expected generated *uncommitted* artifacts, **only after a real catalog is verified**: `catalog_audit.csv` (the complete denominator), `selected_spectra.csv` (every eligible catalog row), `selection_receipt.json` (hashes/counts/known-control matches), `targeted_plan.json`, and `download_receipts.jsonl` plus actual `fits/` bytes if network retrieval succeeds. Do not confuse the planned URLs or locally generated SHA-256s with publisher-certified source-level redshifts. Actual 1D/2D fits, separate detector exposures, NIRCam filter photometry and per-survivor literature priority remain separate **unimplemented validation gates**; they cannot be substituted by the synthetic checks.

The [skeptical referee report](REFEREE.md) formally **rejects scientific acceptance and merging** while the data gates are absent.
