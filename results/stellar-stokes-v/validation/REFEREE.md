# Independent referee check — Stokes V campaign (2026-11-17)

**Verdict: suitable to merge as a bounded catalogue-recovery and partial-image null result, not a discovery.** The result README states what was actually measured and explicitly refuses both planetary-emission attribution and unjustified nondetection limits. No `[REMARKABLE]` designation or log change is supported under `docs/REMARKABLE_FINDINGS.md`.

## Reconciliation against committed numbers

| Check | Independent result and evidence | Verdict |
| --- | --- | --- |
| Frozen target denominators | `samples/receipt.json` has 4,862 planet hosts, 4,480 with verified Gaia motion, 48,809 nearby Gaia, 5,908 CNS5 and 372 astrometric L/T/Y dwarfs. Gzip host/dwarf datasets have the listed uncompressed SHA-256s. 382 hosts lack validated PM/epoch. | PASS |
| V-LoTSS accepted sources | `vlotss_receipt.json`: 68 accepted sources; `vlotss_matches.csv`: 35 on-target pairs [1,20,14,0], **21 unique catalogue names**, zero 120″ offset pairs. `lofar_known.csv`: 20/21 Yiu Table A.2 name overlaps; **i Boo** remains published in the accepted source catalogue despite no name overlap. No new unique source. | PASS |
| RACS-low2 released sources | `casda_receipts.json`: 221 Stokes V rows, 55 schema rows; `racs_matches.csv`: 91 matched target pairs [1,51,39,0] corresponding to **46 distinct source IDs**; every ID appears in `racs_known.csv` with a published `source_association`/`star` type. `racs_offset_matches.csv`: header only, zero control pairs across 32 shifted centres per source/sample. Radio catalogue associations do not resolve individual Gaia components in multiple matches. | PASS, component caveat |
| Exhaustive pair disposition | `match_dispositions.csv` contains 126 rows = 35+91, **all published at radio-source-record level**, not 126 independent stellar emitters; `host_v_measurements.csv` contains GJ 1151 and COCONUTS-2 A. No new-and-unvalidated or verified chance-alignment source identified; alternative Gaia components explicitly unresolved. `disposition_receipt.json` includes per-file SHA-256s. | PASS with source-level qualification |
| Published host V measurements | GJ 1151: integrated signed V +0.45935628 ±0.05893252 mJy, local V RMS 0.058866967 mJy/beam (5 RMS = 0.294335). COCONUTS-2 A / L 34-26: peak V +1.5404146 mJy/beam, local RMS 0.211 (5 RMS = 1.055). Native units and local 5-RMS thresholds reproduced; **neither detected host has a 5σ upper limit.** | PASS |
| ASKAP host footprint geometry | `footprint_receipt.json`: RACS 7,616 + VAST 8,887 = **16,503** selected products, 136,617 matched host–field visit rows, 4,862 host-status rows. Decompressed `host_coverage.csv.gz`: RACS geometric 4,478 / GOOD 4,413; VAST geometric 1,320 / GOOD 1,293; all GOOD-VAST hosts also GOOD-RACS; 4,477 footprint-only, 384 outside both, one published RACS source detection. **No forced host pixel limits.** | PASS; footprint ≠ usable pixel sensitivity |
| VLASS image V | `inventory_receipt.json` 12 inventory V products; `extraction.json` groups and native cutout pixel/header hashes. `measurements.csv`: 13 sample labels, **11 independent** motion-qualified beam positions and zero host positions; `offsets.csv` 104 labelled offsets, **88 independent**; `summary.json`: 0/11 and 0/88 at ≥5σ, maximum |SNR| 2.790 and 2.178, local 5σ 0.891596–1.045182 mJy/beam. CFBDS J005910.90−011401.3 is one non-detection, not 2 because CNS5 277 has an alternative inherited coordinate epoch. | PASS |
| TESS coincidence and planetary phase | `tess_activity.csv`: six host/radio cases, **zero TESS sector overlap**; gaps at least 103 days, 116 positive excursions elsewhere. `periods.csv.gz` is an inventory only. No causal evidence or same-epoch burst/orbit test. | PASS, not validation of an emitter mechanism |
| False-pair accounting | V-LoTSS 0/8 rotated source-centre controls per sample => empirical 0, Poisson 95% upper **0.374** false pairs / sample search. RACS 0/32 => empirical 0, upper **0.094**. VLASS 0/88 independent positions ≥5σ. Source-centred offsets correlate; these are **not calibrated global FDRs** or proofs of physical counterpart identity. | PASS, conditional |
| DR3 access | `dr3_pointings_receipt.json` raw table 2,551 rows; three relevant small Stokes V archive requests returned HTTP **202**, with response hashes. No DR3 pixel-level nondetection, coverage or limit claim. | PASS, pending download |

## Provenance and executable audit

The sources preserve request URLs, timestamps when returned, input SHA-256s and expected row counts in `samples/receipt.json`, `lofar/data/*receipt.json`, `askap/casda_receipts.json`, `askap/footprint_receipt.json`, `validation/known_receipt.json`, `validation/tess_receipt.json` and `vlass/inventory_receipt.json`. Key frozen input hashes include **V-LoTSS FITS** `031d636728022e89b3ff682d07f7c2f150362f279b3c3f794acc665faeb60ccd`, **RACS-low2 V CSV** `74528235bf8dddd8c60c74ff37ab387fa4ceddc46095f9463209aa0e35d75c83`, **host sample uncompressed CSV** `be3ac577ac17bc25d1db46581e3e55becc3e694f8f143d570319bfde496009e8`, and **DR3 pointing CSV** `f6187289fc39b5e1a04c034c409961ed1c187fc61b918f02e81da3bb440f8c31`. `disposition_receipt.json` SHA-256s were independently checked against both new CSVs, and `vlass/extraction.json` records every original-field header hash, row byte range, row SHA-256, cutout pixel hash and matched I/V product. All analysis scripts already committed on the branch were read against their tables; the source search and image pipeline are appropriately described as distinct tests.

For an **offline independent check in a full clone**, run:

```sh
python results/stellar-stokes-v/validation/referee_check.py
python results/stellar-stokes-v/validation/dispositions.py
# Requires numpy, pandas and astropy for the native VLASS pixel/hash replay:
python results/stellar-stokes-v/vlass/measure.py
```

The verifier checks uncompressed SHA-256s including the large field/visit tables, source/sample denominators, offset counts, host flux thresholds, and archived cutout file presence. The disposition script regenerates its own CSVs and hashes. The VLASS measurement script verifies native pixel hashes before reproducing fluxes and controls. **Scope of this in-chat review:** the four smaller compressed sample/host/pointing tables and all text CSVs/JSONs were inspected/recounted, while two >1 MB compressed CASDA field/visit tables could only be checked through the committed row-count/SHA receipts (the GitHub connector did not return their bytes). A full end-to-end offline replay of external archives and those two compressed files was **not run here**; that is not represented as completed independent validation.

## Residual limitations (not merge blockers for this scoped result)

- No host-level upper limits in the 4,477 RACS/VAST footprint-only entries; `GOOD` geometry is not V pixel photometry, a per-host RMS, or leakage-calibrated coverage. The 384 outside selected polygons are no-coverage *for those products*, not necessarily invisible to every radio survey.
- Sample overlap and RACS multiple-Gaia hits preclude component attribution. Publication reconciliation is to the accepted source catalogues/association labels, not an independent second astrophysical confirmation of each mechanism.
- The DR3 Stokes products remain archive-staged; full survey V images, repeated epochs, Stokes I/V beam/leakage validation and simultaneous optical/radio observations remain follow-ups.
- Original survey response caches and complete full-size image fields are not committed; frozen source catalogues, query receipts, matched data, and compressed VLASS cutouts are. The README does not claim complete public-sky sensitivity.

**Referee decision:** approve a **limited null** catalogue entry and merge PR #156 with no remarkable log entry. A *future* signal must independently survive field/beam, leakage, proper-motion, offset controls and repeat/independent V epoch checks before any discovery label.
