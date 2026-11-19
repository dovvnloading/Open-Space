# Nine-position DJA v4.4 lead follow-up — independent-referee decision

**Decision (2026-11-19 UTC): CHANGES REQUIRED; PR #197 stays DRAFT and UNMERGED.** No galaxy newly established at z>12 or z>14; no physical redshift posteriors. The nine objects remain a **targeted candidate-control queue**, not a discovery list or a complete 80,367-spectrum scan.

## Nine dispositions

| Sky group | DJA program/source | Evidence-bounded disposition | Discriminating observation |
| --- | --- | --- | --- |
| 1 | `3073_50252` | Unresolved; no independently validated redshift | No independent physical redshift evidence |
| 4 | `6368_211132` | Unresolved; no independently validated redshift | No independent physical redshift evidence |
| 6 | `6368_151513` | Unresolved; no independently validated redshift | No independent physical redshift evidence |
| 11 | `1181_19411` | JADES DR4 target, quality E; no measured redshift | DR4 `goods-n-mediumjwst_19411`, flag E |
| 16 | `4233_22104` | Unresolved; no independently validated redshift | No independent physical redshift evidence |
| 18 | `4233_7420` | Unresolved; no independently validated redshift | No independent physical redshift evidence |
| 22 | `3215_20095279` | JADES DR4 target, quality E; no measured redshift | DR4 `goods-s-ultradeep_20095279`, flag E |
| 24 | `1286_176143` | JADES DR4 z=3.60656, quality D (tentative, not secure) | DR4 z=3.60656, flag D; 3/20 dusty interlopers pass false-high gate |
| 26 | `1180_b37` | Empty-background shutter B37; not a galaxy | 3/3 SLITS records refer to `background_37`; source RA=Dec=0 |

**Important classification language:** "Published survey target" means an independently catalogued *observation*, **not a published high-redshift galaxy**; DR4 quality **E** means no redshift. "Lower-z" for group 24 is only a DR4 quality **D** *tentative* value; neither secure low-z nor high-z posterior can be assigned. For five other groups, no honest choice among published/lower-z/artefact/new can yet be justified by the evidence. Their status is deliberately **unresolved**. No unvalidated object is called a new z>12 candidate.

## Newly observed runner provenance

- [Runner #44459389472](https://github.com/dovvnloading/Open-Space/actions/runs/44459389472) **succeeded**: downloaded **10/10** specific DJA FITS for the nine groups, inventoried native-looking filenames and FITS extensions, MAST observations, and retrieved **94,124,160-byte JADES DR4 catalogue**, observed SHA-256 `4b4e7fd4306208ee801406c7ed73fdab895a2e67184af7f639a91e5f9c8ab7d4`. The three exact position matches are JADES PID/ID 1181/19411 (flag E), 3215/20095279 (flag E), and 1286/176143 (tentative z=3.60656, flag D); unique JADES tier IDs were retained to resolve duplicate NIRSpec IDs. [JADES DR4 quality-flag documentation](https://jades.herts.ac.uk/DR4/Readme_DR4_catalogues.md).
- [Runner #44460415546](https://github.com/dovvnloading/Open-Space/actions/runs/44460415546) **succeeded**: examined **10/10** additional SLITS tables (3–18 rows each). The purported JADES `1180_b37` is a *background* extraction: every slit-row filename contains `background_37` and its target `source_ra=0, source_dec=0`; its nominal DJA summary sky position is not a valid galaxy identity. This is a **confirmed artefact/control**, not a new discovery.
- Both passes' generated evidence is committed under `outputs/runs/44459389472/` and `outputs/runs/44460415546/`. JSON source-by-source files: [final machine ledger](FOLLOWUP_NINE_2026-11-19.json), [FITS/HDU/DR4/MAST inventory](outputs/runs/44459389472/followup_probe.json), [1D injection simulations](outputs/runs/44459389472/followup_injections.json), [NIRCam attempted catalogues](outputs/runs/44460415546/followup_nircam.json), [per-slit metadata](outputs/runs/44460415546/followup_slits.json).

## Independent adversarial assessment

1. **Detector/slit contamination gate NOT PASSED.** All ten DJA files have coadded SCI/WHT and SLITS provenance but **none has independent detector-integration DQ**. One cannot inspect cosmic-ray, bad-pixel, persistence, cross-shutter contamination or detector-epoch consistency without native MAST exposures and full MSA neighbour geometry. SLITS rows showing one source ID per extraction do **not** rule out adjacent shutters contaminating a spectrum. MAST cone searches establish product availability, not native DQ verification.
2. **NIRCam/SED gate NOT PASSED.** Three specified public catalogue file URLs (GOODS-N, GOODS-S, Abell2744) returned HTTP 404/short bodies, and the other fields were outside this bounded catalogue probe. Thus **zero measured fluxes and zero upper limits** exist at the nine slit locations from this test. No valid filter-convolved Lyman-break, dusty Balmer, brown-dwarf or line-confusion SED comparison was obtained. Failed downloads or no match are never non-detections.
3. **Interloper/posterior gate NOT PASSED.** An explicit independent-parametrization injection/recovery diagnostic ran **20 per family per spectrum** across dusty-Balmer, line confusion, cool-dwarf surrogate, detector spikes, and high-z controls through **the exact same PR #195 1D proxy pipeline**, with original target errors/masks. For **group 24**, **3/20 dusty Balmer injections (15%)** crossed the apparent high-z ranking and two-sided threshold; that source's original z>12 proxy score is **0.997**, despite the JADES **z=3.60656 (quality D)** catalogue estimate. Several spectra recovered **0/20** of injected high-z surrogates at target-matched S/N, so their apparent nulls are sensitivity limitations, not measured evidence for low-z. These 20-object synthetic families are not a representative astrophysical interloper mixture and cannot calibrate `P(z>12|D)` or a credible redshift interval. Prior known false-positive controls (published RUBIES-UDS-144195 and JADES-GN-11508) remain fatal.
4. **Priority gate NOT PASSED.** The JADES DR4 full catalogue provides three source/tier/position matches and establishes the B37 sky-control identity through the FITS. Bounded exact-ID/coordinate public-web searches and prior literature comparison do **not** constitute a comprehensive ADS, arXiv, CAPERS, RUBIES, UNCOVER, CEERS, GLASS, MoM and JADES alias/priority audit across all nine. An absent exact-ID search result does not imply novelty. Source positions and aliases are preserved in the JSON ledger for a future full check.
5. **No survey-wide null.** Only the original 999 targeted spectra were fitted from 80,367 DJA v4.4 catalogue rows; 79,368 were not independently fitted.

**Referee disposition:** **1 empty-sky artefact**, **1 tentative catalogue lower-z**, **2 previously catalogued JADES survey targets with no measured z**, **5 unresolved**. **0 validated new z>12 galaxies; 0 calibrated source redshift posteriors.** Do not add `[REMARKABLE]` or touch `docs/REMARKABLE_FINDINGS.md`. PR #197 **remains DRAFT/UNMERGED; merge SHA: none**.

The correct next work, if requested, requires independently calibrated native 2D/DQ + overlapping-shutter checks, actual NIRCam forced photometry/limits with published filter throughput, exhaustive interloper population injections and redshift model evidence, and an alias-complete ADS/arXiv/individual survey priority review. The repository does not claim those checks are complete.
