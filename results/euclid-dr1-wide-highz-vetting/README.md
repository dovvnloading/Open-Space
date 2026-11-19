# Euclid DR1-Foundation wide-survey VIS-dropout search

**Final result — full public Wide Survey MER HATS catalogue scanned; no new z > 7 quasar validated.** This extends the [seven-cone DR1-F pilot](../euclid-dr1-highz-quasars/README.md), merged as PR #206. It is a *catalogue-screening measurement with an unconfirmed photometric follow-up list*, not a survey of 1,900 deg² of fully quality-masked pixels, a calibrated quasar yield, or a discovery. This does not qualify for the remarkable-findings log.

## What was scanned, and how much sky?

- Dataset: public **12 November 2026 Euclid DR1-Foundation** IRSA Wide Survey MER, [DOI 10.26131/IRSA673](https://doi.org/10.26131/IRSA673). Deep Survey data (IRSA675) are not in the denominator.
- **768/768** HATS leaf partitions, **32 sequential batches of 24**, and **888,152,165/888,152,165** advertised catalogue rows scanned with projected full Parquet row groups. The frozen partition index SHA-256 is **f7940bcaed19b6f67295c5592cc81783d798c0f858062d7293e5e3dfbfa57c2a**. Each batch records source S3 keys, ETags, row-group counts, index/source-column information, strict-row records and output hashes; independent replay passed.
- **Nominal DR1-F Wide Survey imaging footprint: approximately 1,900 deg². Measured unmasked, uniformly sensitive imaging area: unavailable.** The 768 HATS HEALPix *cells* span 3,107.065 deg² geometrically, including partial cells and unobserved sky; this is **not** the imaged area. An HATS row is not proof of clean imaging at any control position, and raw source rows are not independent clean sources.
- The scan used the existing branch-scoped, path-filtered, push-triggered GitHub Actions job. Outputs were checkpointed in order; generated files do not retrigger the job. S3 multipart ETags are recorded as identity receipts, not falsely called whole-object SHA-256 hashes.

## Strict selection, measured counts and control recovery

The **unchanged pilot cut** requires a broad VIS-faint/NIR prefilter, J ≥ 6σ, H ≥ 4σ, |VIS| < 2σ, Y−J ≥ 0.40 (or the specified two-sigma Y lower limit if Y flux ≤ 0), J−H ≤ 0.45, VIS_DET=0, clean J/H flags, ellipticity < 0.30, and no frozen lower J−H limit. Do not retroactively treat subsequent shape/colour diagnostics as part of this frozen selection.

| Stage | Count | Interpretation |
| --- | ---: | --- |
| Published Wide MER rows read | 888,152,165 | All 768 HATS partitions |
| Broad prefilter | 56,662,007 | Permissive, not a quasar sample |
| Strict dropout selection | **462,500** | Photometric/flag selections, not quasars |
| Strict rows with J−H < −0.5 | 62,358 | Pathologically blue, retained for audit |
| Strict rows with −0.5 ≤ J−H ≤ 0.45 | 400,142 | Post hoc colour diagnostic |
| Strict rows with any ground optical band >3σ | 74,209 | Conservative *screen*, individual associations not verified |
| Strict rows in borrowed NISP shape band | 1,375 | Not a validated quasar PSF diagnostic |
| Colour + shape diagnostic + no >3σ ground optical flag, excluding matched known quasars | **478** | Ranked exploratory follow-up list; no confirmed candidates |

**Frozen known controls:** 31 Yang et al. (2026) spectroscopic quasars ([DOI](https://doi.org/10.1051/0004-6361/202658883)), plus 89 previously known z > 6.5 quasars from [Bañados' dated 18 March 2026 census](https://github.com/banados/reionization-era-quasars); 19 spectroscopic dwarf contaminants are separately frozen. Of the 120 known quasars, **24 match within 3 arcsec in the scanned catalogue and 3 pass the strict selection: 3/24 = 12.5%** conditional recovery. Yang-only: 3/22 = 13.6%; catalogue-matched controls at z > 7: **3/8 = 37.5%**. The 96 unmatched controls are **not** counted as failed recoveries: image-pixel coverage and a complete November 2026 census are unverified.

**Major falsification test:** only **1/24** matched known quasars lies in the borrowed **−3.25 ≤ MUMAX_MINUS_MAG ≤ −2.65, ellipticity < 0.2** NISP shape band; **0/3** strict-recovered quasars does. Therefore the band, taken from a published ultracool-dwarf diagnostic, is not a validated quasar morphology filter and the **478-object ranking is severely control-biased**. The catalogue's VIS-based POINT_LIKE_PROB is not valid for VIS dropouts and was never used as such.

## Ranked photometry and quasar-versus-dwarf model

This is a **priority list, not a discovery list**. The fixed rank score weights J signal-to-noise, Y−J, J−H and the *exploratory* NISP shape flag. The leading five by **photometric score**, shown as (RA°, Dec°), are:

| Rank | MER object_id | RA, Dec (deg) | Y−J; J−H (AB) | Score | Two-class model interpretation |
| ---: | --- | --- | --- | ---: | --- |
| 1 | -344582599576564803 | 34.458260, −57.656480 | 3.682; +0.220 | 8.800 | Outside empirical colour support |
| 2 | -640623252224721376 | 64.062325, −22.472138 | 5.174; +0.012 | 8.709 | Outside empirical colour support |
| 3 | -687809343399542971 | 68.780934, −39.954297 | 3.625; +0.331 | 8.579 | Outside empirical colour support |
| 4 | -828861472530046374 | 82.886147, −53.004637 | 3.824; −0.174 | 8.413 | Outside empirical colour support |
| 5 | -698973702349061069 | 69.897370, −34.906107 | 3.281; +0.261 | 8.352 | Near-zero rounded two-class likelihood; not an astrophysical probability |

**Do not assign calibrated quasar probabilities to those extreme-colour leaders.** Separate from photometric rank, the leading five **within supported labelled-colour space** are below. Values are explicitly *conditional on the source being either a quasar or the selected dwarf type*, a hypothesised quasar prior of **0.001**, fixed-bandwidth two-colour KDE, and an unrepresentative matched control set. Under exactly those assumptions, P(dwarf | Q or dwarf) is 1−P(Q | Q or dwarf); **neither is the actual field probability**.

| Rank within model | MER object_id | RA, Dec (deg) | Conditional P(Q) | Conditional P(dwarf) |
| ---: | --- | --- | ---: | ---: |
| 1 | -797174192418604492 | 79.717419, −41.860449 | 0.09488% | 99.90512% |
| 2 | -649582714470637195 | 64.958271, −47.063720 | 0.08956% | 99.91044% |
| 3 | -602426184503595079 | 60.242618, −50.359508 | 0.07311% | 99.92689% |
| 4 | -808747778455215049 | 80.874778, −45.521505 | 0.07311% | 99.92689% |
| 5 | 2801290696627858342 | 280.129070, +62.785834 | 0.07027% | 99.92973% |

The model uses **24 matched labelled quasars and 14 labelled dwarfs**, Gaussian KDE bandwidth **0.32 mag**, with prior-sensitivity table also at 1% Q prior. Its leave-one-out Brier score is **0.629861**, far too poor for purity or population classification; galaxies, artefacts, blends and other stars are omitted. There are **238** nominal two-class conditional calculations, **89** colour out-of-support cases withheld, and the remainder with no defensible numerical model value. These decimals are **not calibrated posterior probabilities**. See [model diagnostics](output/probability_model.json), [probability-supported shortlist](output/probability_supported_candidates.json), [full ranked register](output/ranked_vetted.jsonl.gz), and [top 100 with external receipts](output/top_candidates_with_checks.json).

## Attempted external vetting and what remains open

- **NIR morphology / dwarf locus:** ellipticity and MUMAX_MINUS_MAG were recorded; control tests show the adopted shape band **rejects essentially every known quasar**, so it cannot be used to claim completeness or candidate purity. Empirical two-colour dwarf/quasar likelihoods are uncalibrated as explained above.
- **unWISE W1/W2 / CatWISE:** 24 priority objects were queried. CatWISE 2020 returned **one position-matched catalogue source** with W1−W2(Vega) ≈ **0.606 mag** at MER id **2793077860674821947**, and no >5σ astrometric component. Other CatWISE requests yielded 17 catalogue-empty and 6 service errors; unWISE 2019 had 1 returned-row response with **no 3-arcsec match**, 16 empty and 7 errors. A catalogue-empty result is **not** a forced-flux upper limit; unWISE W1/W2 deblending and object-specific photometric limits remain unresolved.
- **Proper motion:** 24 Gaia DR3 requests returned no source rows at the tested positions. A missing Gaia row does **not** establish zero motion; CatWISE provides only the single above motion diagnostic, without a significant component. No reliable all-object proper-motion exclusion.
- **Image artefacts:** actual Euclid SIA metadata were found for the **top eight**; attempted bounded FITS science cutouts returned HTTP 400, leaving **zero image pixels verified**. VIS ghosts, persistence, image-edge proximity, blends and multi-exposure NISP/VIS defects are **not cleared**. [Pixel-level query/error receipts](output/image_pixel_receipts.json).
- **Literature and public spectra:** SIMBAD queries failed (24/24); DESI, SDSS and ESO metadata requests failed (24/24 each), and MAST returned metadata for 7/24 with 17 errors. This is **not** a clean spectrum/non-association finding. A 2026 all-literature object-by-object check and actual public-spectrum inspection remain open.
- **Error accounting:** all 24 externally queried rows have at least one external-service error. All query results, response hashes where available, returned coordinates, errors, and conditional motion/IR diagnostics remain in [external receipts](output/external_receipts.json) rather than being converted into upper limits or nondetections.

## Expected z > 7 abundance, independent validation, disposition

The [illustrative z > 7 luminosity-function integral](qlf.py) uses a Matsuoka et al. (2023, [ApJL 949 L42](https://doi.org/10.3847/2041-8213/acd69f)) double-power-law M1450 fit, flat ΛCDM (H₀=67.7 km s⁻¹ Mpc⁻¹, Ωm=0.31), k=−0.78 density evolution, **7 ≤ z ≤ 9**, and J_AB ≤ 23 approximated without a detailed SED or bandpass model. It yields **18.710 expected quasars across the nominal 1,900 deg²** before masking and strict-cut completeness; see [numerical assumptions and output](output/qlf_expectation.json). A separate empirical Yang et al. (2026) extrapolation (12 z≥7 in roughly 3,000 deg² scaled to 1,900 deg²) is **7.6**, with different selection and redshift/flux support. Neither number should be compared to the 462,500 dropouts or 478 diagnostic-ranked objects as though they were quasars. There is no effective-area map, injection-recovery model, selection function or meaningful survey rate/upper limit.

The network-free [independent referee replay](referee.py) verified the frozen index and every manifest/candidate hash, row-group partition accounting, all 462,500 stored strict rows against independently implemented cuts, 62,358 blue objects, the 478-object cohort, and 3/24 control recovery. [Referee machine report](output/referee.json) status: **INDEPENDENT_RECEIPT_REPLAY_PASS_NO_DISCOVERY_CLAIM**. The recorded optical-SNR threshold is rechecked at the pipeline's fixed three-decimal reporting precision. [Referee decision](REFEREE.md): acceptable as a reproducible **inconclusive catalogue-screening measurement**, **not** as a validated new z > 7 quasar sample.

**Scientific conclusion:** zero new spectroscopically confirmed z > 7 quasars; no validated discovery, calibrated candidate purity, effective area, clean artifact assessment, or population null constraint. Nothing belongs in [the remarkable-findings log](../../docs/REMARKABLE_FINDINGS.md).

### Reproduce
[RUNBOOK.md](RUNBOOK.md) · [scan.py](scan.py) · [aggregate.py](aggregate.py) · [vet.py](vet.py) · [image_review.py](image_review.py) · [referee.py](referee.py) · [tests.py](tests.py) · [qlf.py](qlf.py) · [32 batch receipts](output/batches/) · [summary](output/summary.json) · [controls](output/controls.json) · [final vetting status](output/vetting_status.json).

This work has made use of the Data Release (DR1) data from the Euclid mission of the European Space Agency (ESA).
