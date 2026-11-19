# JWST DJA + MAST NIRSpec extreme-redshift search (2026-11-18)

**Claim (provisional): No new high-redshift galaxy is established by this investigation at opening.** The intended test is whether public, independently reduced JWST/NIRSpec prism or grating spectra support a galaxy at `z > 14`, or an unreported galaxy at `z > 12`. **A selected source, catalog redshift, photometric dropout, or single continuum discontinuity is not a discovery.** The default conclusion is *unresolved* unless all spectroscopic, image-level, photometric and literature checks below are complete.

## Frozen scientific protocol
- **Primary archive:** public [DAWN JWST Archive](https://dawn-cph.github.io/dja/) NIRSpec 1D/2D prism and grating spectra and associated msaexp redshift fits; supplement only with **public** [MAST](https://mast.stsci.edu/) NIRSpec exposures with explicit retrieval receipts. Record the release/version, access UTC, full selection query, per-file URL/identifiers, count, status, byte size and SHA-256; do not silently substitute releases. The list of all screened spectra must precede adaptive candidate selection.
- **Fit:** evaluate every eligible spectrum on a fixed high-z grid with a Lyman-α (1215.67 Å) break, IGM damping-wing variants and continuum slope/noise nuisance terms; compare integrated marginal likelihood against dusty Balmer/4000-Å break galaxies at `z ~ 3–5`, strong-line masquerades, brown-dwarf/stellar templates and artefact models. Report the *normalized posterior* `p(z|D)` and `P(low-z|D)` under explicit, sensitivity-tested priors. A best-fit z or a photometric z is not a posterior probability.
- **Validate each short-listed target:** inspect native 2D break in distinct integrations and detector/slit quality flags, bad-pixel and contamination maps, neighboring/overlapping MSA shutters and extraction choices; independently check NIRCam fluxes, limits and passbands against the predicted break and all alternate SEDs. Require coherent 1D/2D evidence, not only a pipeline fit.
- **Priority:** search object aliases, coordinates and redshift in ADS, arXiv, DJA high-z lists and JADES, CEERS, UNCOVER, MoM and GLASS literature; record query/date and known-vs-new judgment for each surviving target. A literature non-match is not proof of originality.
- **Referee gate:** no `[REMARKABLE]` tag or entry in `docs/REMARKABLE_FINDINGS.md` unless the independent checks have established discovery-level evidence. Preserve failures and limitations.

## Frozen public-source scope and what was actually checked

The most recent DJA compilation identified in the public release pages at this freeze is **v4.4** ([Zenodo 10.5281/zenodo.15472354](https://zenodo.org/records/15472354); documented DJA update 5 September 2025). Its **80,367 unique spectrum records** comprise **34,949 prism records** plus **45,418 grating records** across eight grating/filter groupings. **These are spectra, not distinct galaxies or independent exposures.** The published 8.4 GB deposit provides one summary catalog, a column description and nine combined 1D spectral FITS tables. The exact eleven publisher filenames, URLs and reported MD5 hashes are committed in [`frozen_sources.json`](frozen_sources.json).

A **separate** DJA newly-public MAST dashboard snapshot listed 16 NIRSpec program/target/mode rows (the displayed grouped counts sum to 239). They are *observing metadata only*: some are stars, Solar System or transients, and the dashboard rows are neither public-product identifiers nor unique galaxies. The source list also records four high-z-labelled programs, GO-6954, GO-7014, GO-9180 and GO-9645, but **none has been extracted, crossmatched, downloaded or spectroscopically validated here**. The `mast-index` script is a future executable public-product inventory, not evidence that products have already been retrieved.

**Observed execution status as of deposit:** The ChatGPT execution environment could not resolve the public archive hosts for the multi-GB retrieval; zero DJA FITS files, zero source-level manifest rows and zero MAST product bytes have been downloaded and SHA-256 verified. No 1D fits of actual galaxies were performed. The dataset is scoped from publisher metadata, not measured from local spectra. The root README lists this **provisional, uncounted preflight** so that it cannot be mistaken for a completed negative survey.

### Reproduction and required evidence (outside this restricted environment)

Run from the repository root with sufficient network, storage (source files total approximately 8.4 GB, plus results), CPU and memory. No raw FITS data or source catalog bytes are vendored into Git:

```sh
python results/jwst-dja-nirspec-highz-screen/verify_preflight.py
python -m pip install numpy scipy astropy astroquery
python results/jwst-dja-nirspec-highz-screen/screen.py --self-test
python results/jwst-dja-nirspec-highz-screen/retrieve.py download --all --cache /data/dja-v4.4
python results/jwst-dja-nirspec-highz-screen/retrieve.py freeze --cache /data/dja-v4.4
python results/jwst-dja-nirspec-highz-screen/retrieve.py mast-index --cache /data/dja-v4.4
python results/jwst-dja-nirspec-highz-screen/screen.py --cache /data/dja-v4.4 --out /data/dja-v4.4/fit
```

The retrieval scripts must reject any publisher-MD5 mismatch; successful execution would create **new observed SHA-256 receipts** in the local cache (not yet committed). The source-row freeze must check all 80,367 entries and preserve exact FITS-aligned order. The 1D code does not use `z_best` or DJA grade as an input to the fit. Its fitted `p_gt14_proxy` and `p_lowz_gal_proxy` are **illustrative model-comparison scores**, *not* calibrated physical probabilities: the broad source/noise/stellar/line priors and aperture correlations require extra work. The `--self-test` uses synthetic arithmetic invariants, not recovery/false-alarm controls or a real spectrum. CI only runs preflight and numerical invariants; it does **not** run a real-sky scan.

For publication-grade claims, the native 2D pixels, MSA shutters/overlaps, bad-pixel masks and disjoint integrations must be brought into evidence, and all candidate model fits convolved with source-associated NIRCam fluxes and limits. The provenance, individual source results, candidate list, photometry and per-object literature decisions have **not** been produced. See the referee's [`REFEREE.md`](REFEREE.md). Published positive controls and their finite, explicitly incomplete [literature baseline](literature_baseline.json) prevent misreporting already published galaxies as discoveries; no object-specific priority search was performed.

## Deliverables and audit state
| Gate | Opening status |
| --- | --- |
| Publisher release/version and nine-shard counts | COMPLETE FROM PUBLISHER METADATA; NOT INDEPENDENTLY VERIFIED |\n| Complete source-level spectrum manifest | PENDING: raw catalog never retrieved |
| Publisher MD5/URLs and dashboard metadata | FROZEN in JSON; NOT local download receipts |\n| Per-input observed SHA-256 and HTTP/MAST API query receipts | PENDING: no bytes retrieved |
| Per-spectrum independent high-z and interloper posterior fits | PENDING |
| Native 2D / independent exposure and slit-contamination inspection | PENDING |
| NIRCam consistency and literature crossmatch | PENDING |
| Independent skeptical review | COMPLETE as critical documentation; review REJECTS scientific claims |

The source-*scope* inventory, scripts, published-checksum ledger, literature baseline and dated referee report are deposited. A [fail-closed dossier gate](candidate_gate.py) rejects claims without SHA-256 provenance for separate 2D exposures, pixel/slit checks, calibrated alternate-redshift evidence, NIRCam filter/flux records and an ADS/arXiv/DJA/JADES/CEERS/UNCOVER/MoM/GLASS priority receipt; it **does not** substitute for performing those checks. Native source-level inventories, observed raw hashes, 1D fits, 2D receipts and object-specific photometry/literature records are **still absent**. The only warranted outcome is **incomplete coverage and an unexecuted astrophysical screen**, not a null detection or a discovery.

**Separation of work:** this campaign uses only JWST NIRSpec DJA/MAST archival spectroscopy; it does not touch the COSMOS2025 stellar-mass/baryon-bound question in PR #194.
