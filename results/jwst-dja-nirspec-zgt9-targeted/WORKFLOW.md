# One dedicated manual GitHub Actions runner for PR #197

The default branch already contains the workflow
`.github/workflows/jwst-dja-nirspec-preflight.yml` with
`workflow_dispatch`. PR #197 extends this EXISTING file with a single
gated science job. It can therefore run on the PR branch **without**
putting the research code on main first.

To launch: GitHub repository → **Actions** → **JWST DJA NIRSpec
high-z preflight** → **Run workflow** → select
`research/jwst-dja-v44-zgt9-targeted` → **Run workflow**.
No extra input is needed. This is deliberately **manual only**: it
does not start on every push. `workflow_dispatch` has no native path
filter; its one research job is guarded by the exact branch ref and
touches only `results/jwst-dja-nirspec-zgt9-targeted/`.
The legacy PR #195 preflight is not run on manual dispatch.

The job uses the GitHub-hosted network, not the networkless Python
sandbox. It fetches the **exact v4.4 summary CSV.gz** (publisher MD5
`721c58c446f4a4eb1570aef02c4ab47b`), verifies all 80,367
catalog rows, and selects any spectrum with a DJA redshift hint >9
plus published z>12 positive controls. It then fetches **only those
individual .spec.fits** from the DJA S3 HTTPS host, with the
documented alternate CANUCS bucket where appropriate; **no
8.4-GB Zenodo spectral shard is downloaded**. Each GET, SHA-256,
failure and FITS check is recorded.

The 1D model reuses PR #195's *exploratory* high-z
Lyman-break/wing versus dusty Balmer, line-confusion, stellar/dwarf
and artefact template comparisons, now including z=9–10; the
numerical scores are **NOT calibrated astrophysical probabilities**.
For every catalog z>12 solution, matched published control and
proxy-supported z>12 review item, it attempts a figure displaying
the actual FITS SCI (and WHT, when available) 2D pixels plus 1D
spectrum. A co-added extraction is NOT a disjoint independent
exposure and the 2D x-axis may be approximate. If SCI is absent,
the record explicitly says so.

A compact 80,367-row denominator, target selection list, observed
download and redshift-fit receipts, all available 2D figures, and an
object-by-object z>12 priority-review queue are committed **once**
to `outputs/` on the PR branch by the workflow's scoped token.
The full audit is additionally an Actions artifact; raw FITS and
the 133 MB summary are **not** committed. All non-control objects
remain `UNREVIEWED_NOT_NEW`. The job does NOT invent an ADS/arXiv
priority verdict, independent native-exposure/DQ check or NIRCam
passband-convolved photometric validation.

**Blocking behavior:** A failed publisher MD5, DNS/HTTPS GET,
partial individual FITS retrieval or numerical error is recorded in
`outputs/run_receipt.json` and the job fails. A complete unvetted
1D run may still have a green workflow, but **PR #197 remains a draft**
until a separate referee has examined source-specific 2D/native
exposures, NIRCam flux limits, interloper priors and the literature
for *every* surviving z>12 case. No `[REMARKABLE]` or discovery
claim follows from this automation.
