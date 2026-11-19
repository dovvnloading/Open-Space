# Euclid DR1-Foundation Wide: high-recall control recovery and source-level photometric follow-up

**Claim — measured but inconclusive, no new quasar:** the original Euclid DR1-F strict VIS-dropout screen missed most known quasars. The frozen, control-tuned relaxed *eligibility* cut recovers **21/24 (87.5%)** catalogue-matched spectroscopic z>6.5 quasars instead of **3/24 (12.5%)**, including **7/8 at z>7**; it retains **4/14** matched labelled dwarf contaminants and yields **1,959,717 raw eligible MER catalogue rows**. This is a more sensitive catalogue *triage*, **not** a purity-controlled quasar sample, an 87.5% complete 1,900-deg² survey, or a discovery.

## Why each control failed

The [complete per-object failure table](output/control_failures.tsv) reports all 24 matched quasars' old/new cut outcomes, VIS_DET, signed 3-FWHM-aperture VIS/J/H S/N, Y−J and J−H, flags, ellipticity and NISP MUMAX_MINUS_MAG; [raw source rows and ETags](output/control_audit.json) make every reason independently reproducible. The audit re-read **35 exact, ETag-verified HATS source partitions**, recovered all **24 quasars and 14 dwarfs**, and independently reproduced the old 3/24. Original cut failures *overlap*: Y−J≥0.40 fails **20/24**, |VIS S/N|<2 fails **10/24**, VIS_DET=0 fails **2/24**, J/H clean-flag conditions fail **2/24**, and ellipticity<0.30 fails **3/24**. A borrowed narrow NISP shape band retains only **1/24**, so it is not a validated quasar morphology discriminator. Of the 120 frozen dated literature quasars, **96 lack a catalogue match within 3 arcsec** and are excluded from the recovery denominator because the released clean-image footprint is not verified.

The [frozen v1 rule](freeze.json) requires J flux>0.65, J≥6σ, H≥4σ, |VIS 3-FWHM aperture S/N|<7 (without hard VIS_DET veto), Y−J≥0 (or the 2σ Y lower limit when Y≤0), −0.20≤J−H≤0.48, ellipticity<0.40, −1.45≤NISP MUMAX_MINUS_MAG≤−0.75, and spurious probability <0.7 or missing. J/H flags, ground optical >3σ, and VIS_DET are preserved for ranking/review rather than silently turning known quasars into rejects. It misses three matched controls: EUCL J052209.82−512709.2 (J−H≈−1.41), EUCL J091639.93+683652.9 (VIS S/N≈16.3, abnormal morphology), and J024401.02−500853.74 (NISP MUMAX≈−3.06).

**Conditional catalogue-matched recall: 21/24 = 0.875; Wilson 95% interval 0.690–0.957. Labelled dwarf retention: 4/14 = 0.286; Wilson 95% interval 0.117–0.546.** The cut was *tuned on these very same labels*, so this is resubstitution, **not held-out or injection-recovery completeness**. The 14 known dwarfs are not a representative survey contaminant population; the 28.6% retention is **not field purity**. See [machine calibration](output/calibration.json) and [referee replay](output/referee.json).

## Full-source re-selection

All **768/768 original Euclid DR1-Foundation Wide MER HATS leaves and 888,152,165/888,152,165 rows** were re-read, using the original S3 source keys and ETags (not merely the old 462,500 strict rows), in **32 committed sequential batches of 24 partitions**. Every batch's top register and source count is hashed, with one path-filtered, push-triggered networked Actions job and output commits excluded from its trigger. There were **1,959,717 eligibility rows**, including **708,402** without a >3σ ground optical flag and **1,796,764** without a nonzero J/H flag. The exact global photometric top 200 is preserved by the 200 highest scores per partition; the broader retained rank register has **147,522** rows (and excludes two literature-matched control IDs for candidate follow-up). Those counts are **not quasar counts**.

[Full survey summary](output/summary.json) · [32 source batches](output/batches/) · [full ranked top register](output/ranked_candidates.jsonl.gz) · [top 100](output/top_candidates.json).

## Ranked photometric priorities, not validated quasars

The first eight **by the frozen heuristic photometric score** are below. Y−J/J−H are AB colours; conditional Q-vs-dwarf likelihoods use a hypothetical Q prior of 0.001. A rounded 0 means *vanishingly small in that arbitrary two-class KDE*, not zero true quasar probability; “outside” means a numerical value is withheld.

| Rank | MER object ID | RA°, Dec° | Y−J; J−H | Two-class P(Q), uncalibrated | Mosaic visual screen |
| ---: | --- | --- | --- | --- | --- |
| 1 | -615262033606966661 | 61.526203, −60.696666 | 2.751; +0.118 | ≈0 (rounded) | Nearby extended source; association pending |
| 2 | -642780668608958826 | 64.278067, −60.895883 | 5.041; +0.123 | outside support | Bright linear streak nearby |
| 3 | -464108856519632999 | 46.410886, −51.963300 | 3.075; +0.144 | ≈0 (rounded) | Bright neighbor; deblending pending |
| 4 | -580697145407050206 | 58.069715, −40.705021 | 3.215; +0.146 | ≈0 (rounded) | Compact central NIR signal |
| 5 | -728937611292911060 | 72.893761, −29.291106 | 2.378; +0.169 | 1×10⁻⁸ | Central J/H PSF/ring structure |
| 6 | -464123285327831920 | 46.412329, −32.783192 | 4.753; +0.063 | outside support | Compact NISP source; flux check |
| 7 | -758338901467120575 | 75.833890, −46.712058 | 3.224; +0.064 | ≈0 (rounded) | Dark NISP stripe at source |
| 8 | -318088655326522487 | 31.808866, −32.652249 | 3.903; +0.191 | outside support | Structured noisy background |

The separately sorted **top three by conditional Q-vs-dwarf model probability *within only these 100 photometric ranks*** (not a global or calibrated quasar ranking) are: MER **-526096119485738264** (photometric rank 67, P(Q|Q-or-dwarf)=2.06×10⁻⁶), **2772682462555347715** (rank 87, 1.61×10⁻⁶), and **-51009972626798953** (rank 47, 1.03×10⁻⁶). These numbers are conditional to an arbitrary **0.1% quasar prior**, a Gaussian bivariate colour KDE of bandwidth **0.32 mag**, and only **24 quasar plus 14 dwarf labels**; a 1% prior is separately tabulated. The leave-one-out Brier score is **0.629861**, indicating a very poor classifier. **46/100** photometric ranks lie beyond four bandwidths from any labelled colour and have **no numerical model probability**. Omitted populations include galaxies, other stars, artefacts and blends. These decimal values are **not credible field posteriors or candidate confirmations**.

[Top 100 with labelled-model values and external receipts](output/ranked_conditional_top100.json) · [model diagnostics](output/conditional_model.json) · [conditional-supported subranking](output/conditional_supported_top100.json).

## Images, mid-IR, literature, spectroscopy — actual evidence and open gates

- **Real Euclid science pixels obtained:** IRSA SIA returned HTTP 502 for this run, but [public S3 MER tile lookup and bounded Astropy FITS section reads](pixel_review.py) succeeded for **VIS, Y, J and H for all top eight** (32 actual FITS mosaic cutouts, each 251×251 science pixels). [Per-band SHA-256 and WCS/pixel metadata](output/image_pixel_receipts.json) and [viewable contrast-stretched PNGs](output/pixel_previews/) are committed. Qualitative [visual-review flags](output/visual_review.json) include a bright streak near rank 2, unusual J-band annular/PSF structure at rank 5, and a **dark NISP stripe near rank 7**. **0/8 are exposure-level artefact-cleared:** multi-exposure persistence, ghosts, masks, edge geometry and quantitative PSF/deblending need review. Independently stretched PNGs are not cross-filter flux ratios.
- **CatWISE/unWISE and motion:** 16 ranked sources had IRSA TAP, alternate NOIRLab CatWISE2020 and unWISE Gator attempts. **Zero position-matched CatWISE W1/W2 colours and zero position-matched unWISE fluxes were obtained**; some requests failed, and catalogue-empty responses do **not** establish forced-flux upper limits. No significant-motion rejection was measured; **absence of a catalogue motion measurement is not evidence of zero motion**. [Primary external receipts](output/external_receipts.json) · [WISE alternate-route receipts](output/wise_mirror_receipts.json).
- **Literature and public spectra:** SIMBAD, Gaia, SDSS/DESI/ESO/MAST and other services were queried at the leading positions, but there are **service errors, unavailable associations and only metadata**, with **zero actual source spectra reviewed and zero new spectroscopic redshifts verified**. The missing-match statements are not negative spectroscopy. Detailed service/HTTP statuses are preserved in [external receipts](output/external_receipts.json) and [summary](output/external_summary.json).

## Expected z>7 objects and referee disposition

The original [Matsuoka-style luminosity-function integral](../euclid-dr1-wide-highz-vetting/output/qlf_expectation.json) predicts **18.710** z=7–9 quasars for a *hypothetical complete, uniformly sensitive* 1,900 deg² Wide release to J_AB≤23; a distinct Yang et al. area scaling predicts **7.6**. The **~8–19** expectation cannot be compared with the 1,959,717 photometric eligibility rows as a measured quasar yield; the effective unmasked area, injection completeness, depth/SED and purity are unmeasured. **No calibrated survey upper/lower bound follows.**

The independent [source receipt and scalar-control referee](referee.py) returns [INDEPENDENT_RECEIPT_PASS_PHOTOMETRIC_ONLY_NO_DISCOVERY](output/referee.json): 768 partitions, all published rows, source hashes, 21/24 controls, 4/14 dwarfs, rank-preserving top lists, conditional model support gating and 8/8 real image targets reconcile. [Referee decision](REFEREE.md): this is a reproducible *measured/inconclusive catalogue-control correction and photometric follow-up list only*. **No new z>7 quasar is confirmed. No remarkable-findings entry.**

Euclid ESA DR1-Foundation public Wide MER data via NASA/IPAC IRSA; [previous strict full-index study, PR #208](../euclid-dr1-wide-highz-vetting/README.md). These results and interpretation concern the dated 2026 control catalogues and cannot substitute for a future complete survey selection function.
