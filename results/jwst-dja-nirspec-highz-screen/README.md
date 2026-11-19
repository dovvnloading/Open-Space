# JWST DJA + MAST NIRSpec extreme-redshift search (2026-11-18)

**Claim (provisional): No new high-redshift galaxy is established by this investigation at opening.** The intended test is whether public, independently reduced JWST/NIRSpec prism or grating spectra support a galaxy at `z > 14`, or an unreported galaxy at `z > 12`. **A selected source, catalog redshift, photometric dropout, or single continuum discontinuity is not a discovery.** The default conclusion is *unresolved* unless all spectroscopic, image-level, photometric and literature checks below are complete.

## Frozen scientific protocol
- **Primary archive:** public [DAWN JWST Archive](https://dawn-cph.github.io/dja/) NIRSpec 1D/2D prism and grating spectra and associated msaexp redshift fits; supplement only with **public** [MAST](https://mast.stsci.edu/) NIRSpec exposures with explicit retrieval receipts. Record the release/version, access UTC, full selection query, per-file URL/identifiers, count, status, byte size and SHA-256; do not silently substitute releases. The list of all screened spectra must precede adaptive candidate selection.
- **Fit:** evaluate every eligible spectrum on a fixed high-z grid with a Lyman-α (1215.67 Å) break, IGM damping-wing variants and continuum slope/noise nuisance terms; compare integrated marginal likelihood against dusty Balmer/4000-Å break galaxies at `z ~ 3–5`, strong-line masquerades, brown-dwarf/stellar templates and artefact models. Report the *normalized posterior* `p(z|D)` and `P(low-z|D)` under explicit, sensitivity-tested priors. A best-fit z or a photometric z is not a posterior probability.
- **Validate each short-listed target:** inspect native 2D break in distinct integrations and detector/slit quality flags, bad-pixel and contamination maps, neighboring/overlapping MSA shutters and extraction choices; independently check NIRCam fluxes, limits and passbands against the predicted break and all alternate SEDs. Require coherent 1D/2D evidence, not only a pipeline fit.
- **Priority:** search object aliases, coordinates and redshift in ADS, arXiv, DJA high-z lists and JADES, CEERS, UNCOVER, MoM and GLASS literature; record query/date and known-vs-new judgment for each surviving target. A literature non-match is not proof of originality.
- **Referee gate:** no `[REMARKABLE]` tag or entry in `docs/REMARKABLE_FINDINGS.md` unless the independent checks have established discovery-level evidence. Preserve failures and limitations.

## Deliverables and audit state
| Gate | Opening status |
| --- | --- |
| Archive release/version and complete frozen spectrum manifest | PENDING |
| Per-input hashes and HTTP/query receipts | PENDING |
| Per-spectrum independent high-z and interloper posterior fits | PENDING |
| Native 2D / independent exposure and slit-contamination inspection | PENDING |
| NIRCam consistency and literature crossmatch | PENDING |
| Independent skeptical review | PENDING |

Scripts, frozen inventories, raw receipts, per-spectrum results, and a dated referee report will be deposited here as work is executed. If the released products cannot be obtained or validated, the only warranted finding is **incomplete coverage**, not a null detection or a discovery.

**Separation of work:** this campaign uses only JWST NIRSpec DJA/MAST archival spectroscopy; it does not touch the COSMOS2025 stellar-mass/baryon-bound question in PR #194.
