# Targeted DJA v4.4 z>9 NIRSpec screen — blocked data-access continuation

**Claim (2026-11-19 UTC): No astronomical spectrum was downloaded, fitted, or validated in this continuation. There is NO new z>12 discovery, no calibrated interloper probability, and no survey-wide null result.** This is a blocked continuation of [merged PR #195](../../../../pull/195), not a replacement of its frozen 80,367-spectrum scope.

## Why the targeted analysis is blocked

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
