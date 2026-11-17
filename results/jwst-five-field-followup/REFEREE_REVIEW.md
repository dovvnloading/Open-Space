# Independent referee check — five public JWST fields (PR #146)

**Disposition: bounded null; two sources remain unresolved, and no new companion is validated.** This review audits the *frozen* measurements independently of the narrative, and intentionally does not turn a weak single-visit source into a candidate planet.

## Recomputed quantities

| Claim | Referee check against the deposited evidence | Result |
|---|---|---|
| LAWD 66 1.207″ motion | Difference between the first and last F164N host centroids in `lawd66/epochs.csv`, **adding** the per-epoch east/north field-tie translations, yields 1.206612″ over 1,140.916 days. The *unregistered* centroid displacement is 1.221772″ and must **not** be substituted for the field-relative motion. | Matches rounded claim |
| LAWD 66 stationary sources | `lawd66/peaks.csv` has 93 triggers: 92 have a fixed-sky win and none a host-motion win; these 92 cluster into nine groups in `lawd66/sources.csv`. The tenth is one F164N peak at 3.273864″ (E +0.343698″, N +3.255773″ relative to the star) on 2025-04-21. | Matches |
| LAWD 66 native dithers | `lawd66/dither_photometry.csv` gives four positive apertures (0.746, 0.420, 1.297, 0.795 µJy) with zero flagged aperture pixels; only one **visit** contains this peak. The other epochs do not establish co-motion or a constant 0.802 µJy object. | **Unclassified** |
| HD 53143 archive association | The two rows of `nircam-hd53143-hd38206/outputs/NICI_track.csv` give residuals of 0.149391″ and 0.149290″ for the historical fixed-background prediction versus 4.672510″ and 4.683494″ for co-motion. `known_source_color.csv` gives F200W−F444W = −1.1658 and −1.1682 AB. | Background; consistent with published Gemini/NICI source |
| HD 53143 MIRI | Each November 2025 F2550W roll and its combined product has an empty `peaks.csv` across the stated 3–9″ search. Both rolls are one epoch; the public release was November 2026. | Bounded nondetection |
| HD 38206 | Five 2025-01 and three 2025-10 F444W detections in the two `peaks.csv` tables. Two cross-epoch coincidences are at precisely (98,215) and (166,235) in **detector** pixels rather than on the sky; the pipeline detector-quality mask flags those positions. All eight co-moving predictions have second-epoch 3×3 scores below 2.36. Sixteen predicted historical loci were sampled in 20 Gemini H-band frames each. The 35.4 mas fixed-versus-moving baseline does not by itself discriminate the hypotheses. | No repeat companion |
| HD 21997 | `HD21997_peaks.csv` has zero F444W detections. Auxiliary F200W and MIRI detections are reference/detector dependent, not reproducible same-sky sources. | Bounded nondetection |
| HD 202628 | Two F444W peaks (`HD202628_peaks.csv`) fall from 5.44/8.16 to 0.31/−0.50 when the reference changes. The F1140C peak at about 3.283″ has a 6.12 local-MAD pipeline response in one roll that becomes −0.77 with a different native PSF reference. Two same-day rolls supply no proper-motion baseline. | **Unclassified** |
| HD 202628 F1140C sensitivity | Independently intersecting eligible angles by trial ID in the two-roll `F1140C_recovery_trials.csv` gives **9/28** and **0/28** joint recoveries at 10 µJy; **28/28** at 15 µJy with HD 201865 as reference and **28/28** at 20 µJy with HD 22333. The observed 2–3 µJy proxy flux is below both high-recovery limits. | Matches; cannot reject a weaker source |
| Contrast curves in all five fields | Grouping the original per-angle injection ledgers by radius and tested flux, reapplying each method's 90% or 95% rule (including the F444W brighter-level consistency requirement), reproduces **all 57 rows** of `contrast_summary.csv` without a mismatch. These are 14,256 F444W, 3,138 F2550W, 7,945 LAWD 66 single-image and 1,470 paired-motion injections: **26,809 trial rows**. | All five have empirical real-image recovery curves |

The reproducible assertions are in `../../tests/test_jwst_146_referee.py`; run from the repository root:

```sh
python -m unittest discover -s tests -p 'test_jwst_146_referee.py' -v
python -m unittest discover -s tests -p 'test_jwst_five_field_followup.py' -v
```

## Receipts, reproducibility and scope

The referee test rehashes the deposited LAWD 66 MAST query/receipt bytes against `lawd66/inputs.json`, the HD 53143 MIRI CAOM response against `hd53143-miri/receipts/query.json`, and the four HST archive responses against `archives/query_receipts.json`. The HST track ledger independently contains **504, 32, 12 and 168** tested combinations for HD 53143, HD 38206, HD 21997 and HD 202628, respectively (**716** total). The NIRCam, MIRI and detector products have frozen SHA-256 input manifests and MAST identifiers within their field folders.

**Important boundary of the referee rerun:** the repository deposits scripts, small machine-readable measurement and injection ledgers, hashes, and archive responses, **not** all the calibrated raw FITS image bytes. LAWD 66 alone uses 32 products totaling about 5.36 GB. The independent offline run recomputes the statistics, flags the known detector-location coincidences, checks the response hashes and validates each deposited recovery curve. It does not constitute an independent download and pixel-level regeneration of every MAST FITS reduction; that requires fetching the hashed public inputs and running the field-specific scripts. The GitHub claims verifier's generic offline replay also reports unresolved shell-variable script paths and some network/cache drift in the **older #142 project**, even when the workflow finishes successfully. None of those warnings should be misrepresented as a fresh independent raw-FITS replay. These are reproducibility limits, not evidence for a new object.

## Required observations to settle the two unclassified signals

1. **LAWD 66:** A genuinely *later* NIRCam **F164N** dither sequence, for example **2028-04-21**, at the predicted host-relative co-moving offset **E +0.344″, N +3.256″** (3.274″ separation) from the Gaia-propagated white dwarf. Obtain independent image-level detections and 90% recovery for an injected ≈0.8 µJy source at that position, and check the competing **fixed-sky** prediction at RA **257.863643°**, Dec **−14.799260°** registered in 2025. Use updated proper motion/parallax and WCS for the actual observing date; offsets assume negligible orbit curvature. An independent centroid at one of these discriminating positions, or a sufficiently deep nondetection, would settle the present single-visit ambiguity.

2. **HD 202628:** A later two-roll MIRI **F1140C** sequence, for example **2028-06-02**, testing the host-relative co-moving offset **E +1.178″, N +3.065″** (3.283″ separation) from the propagated HD 202628 position **and** the 2023 fixed-sky position. Require reference-star substitution, native-exposure persistence, astrometric calibration, and real-image **joint** injection recovery at **≤2 µJy** at the tested position. The existing 15–20 µJy limit cannot settle a 2–3 µJy source; a later shallow nondetection would not resolve it.

Dates are **proposed tests, not confirmed JWST schedules**. These positions specify *star-relative* co-motion, not a guaranteed precise absolute WCS ephemeris. Neither signal meets the discovery inclusion bar; no `[REMARKABLE]` title or change to `docs/REMARKABLE_FINDINGS.md`.
