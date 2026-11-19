# DJA v4.4 z>9 NIRSpec — manual and automatic Actions runner

The existing GitHub workflow `.github/workflows/jwst-dja-nirspec-preflight.yml` is configured on the PR branch for **one** networked job `targeted-dja-v44-zgt9`. It uses GitHub-hosted Actions runners, not the offline local Python sandbox.

## Triggering

- **Manual:** Actions → **JWST DJA NIRSpec high-z preflight** → **Run workflow** → branch `research/jwst-dja-v44-zgt9-targeted`. This worked for [run #44454883069](https://github.com/dovvnloading/Open-Space/actions/runs/44454883069).
- **Automatic future runs:** a `push` to `research/jwst-dja-v44-zgt9-targeted` changing files under `results/jwst-dja-nirspec-zgt9-targeted/**` or the workflow file triggers the **same** job. **Output-only changes under `results/jwst-dja-nirspec-zgt9-targeted/outputs/**` are excluded by the negative path pattern**, preventing the Actions bot's automatic commit from rerunning itself. The legacy #195 preflight job is skipped on pushes to the targeted branch.
- `workflow_dispatch` does not support native `paths` filters; the workflow applies branch/job guards. There is one research job, no matrix or separate downloader/fit jobs. The Actions run deposits receipts with one Git commit/push using a narrowly scoped `contents: write` token.
- Future runs write to **`outputs/runs/<run_id>/`** rather than overwriting the original [run #44454883069](outputs/run_receipt.json) output set at `outputs/`. Its observed job finished successfully and committed 999 FITS download receipts, 999 unvetted 1D proxy fits and 64 plotted coadded 2D figures.

## Fixed data and computation

The code downloads only the **133,004,545-byte** publisher-MD5-verified DJA v4.4 **summary CSV.gz**, checks all 80,367 rows, and selects `z_best>9`, `z_phot>9`, and six literature baseline controls. It downloads only those individual `spec.fits` products by verified DJA S3 URL, using a documented alternate CANUCS bucket where applicable, with HTTPS GET, byte counts, SHA-256 and FITS HDU validation. The initial successful run selected and downloaded 999 spectra totaling **2,038,757,760 bytes**. The nine combined Zenodo spectral shards are **not** downloaded.

The 1D model uses the prior PR #195 Lyman continuum break/damping-wing versus dusty Balmer/4000 Å, line-confusion, approximate stellar/brown-dwarf and artefact templates; it is an **uncalibrated proxy**, not validated posterior odds. Its 2D figures show individual coadded SCI (and WHT, if present) pixels, not independent detector exposures. The source names, 1D fit values, observed download digests, 2D images and a z>12 review queue are committed; large raw FITS and the 133 MB original catalog are left outside Git.

## Data gate and scientific release gate

The workflow emits `run_receipt.json` even on GET/MD5/fitting failures, and marks a run failed on blocked retrieval or partial downloads. A green run only means the **observed scripted retrieval and unvetted 1D computation** completed. It cannot establish new z>12 galaxies, a fully sampled survey null, NIRCam red/blue SED support, clean independent native exposures or literature novelty.

The [referee report](REFEREE.md) and [27-position identity ledger](literature_referee_2026-11-19.json) found published lower-z false positives including RUBIES-UDS-144195 and JADES-GN-11508, as well as a too-wide 0.5″ positive-control crossmatch. The next selector uses **0.2″** with an explicit nearby-source regression test. No `[REMARKABLE]` entry. **PR #197 remains draft/unmerged until full independent validation.**
