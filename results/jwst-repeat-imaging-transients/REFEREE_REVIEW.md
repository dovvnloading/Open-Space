# Independent referee review — JWST repeat imaging, PR #129

Review date: **2026-11-17 UTC**. This is a review of the existing `research/jwst-repeat-imaging-transients` branch, not a new search or a new PR. Neither source is claimed as a discovery or assigned a transient type.

## Numerical replays on deposited inputs

`tests/test_jwst_pr129_referee.py` uses a scratch copy of the checked-in scientific tables, executes **eight** of the original scripts, then compares their regenerated outputs with the frozen CSV/JSON results (numeric tolerance for serialization):

1. COSMOS-022 night/filter census: `cosmos-022-history/summarize.py` (nightly, medium-band changes, inventory hashes).
2. Pre-event MUSE narrow-line fit: `cosmos-022-validation/emission_lines.py` (line fluxes).
3. MUSE broad Hβ limits: `cosmos-022-validation/broad_hbeta.py`.
4. Full catalog-trigger dispositions: `repeat-catalog-followup/dispositions.py`.
5. COSMOS-022 host-dependent difference magnitude: `repeat-catalog-followup/candidate_photometry.py`.
6. GOODS-S blind-area, recovery and paired-detection summary: `goods-s-oct2026/blind_summary.py`.
7. GOODS-S host photometric redshifts and difference magnitudes: `goods-s-oct2026/nuclear.py`.
8. NEXUS same-epoch, matched-filter contrasts and injection recovery: `high-z-hosts/join_bands.py .`.

The independent unit assertions confirm the headline claims in PR #129 and the main result READMEs:

| Deposited result | Referee check |
|---|---|
| Repeat-catalog screen | 69 frozen catalogs with 317,286 source entries; 37 matched filter comparisons with 47,584 eligible entries; 463 single-filter triggers; 42 grouped candidates, all disposed; 100 image/position/filter checks, 14 groups reproducing in two filters at both apertures |
| COSMOS-022 | 105 date/filter samples across 20 dates and 20 filters; all eight medium-filter changes positive; native F162M medians 611→843 nJy and F182M 624→824 nJy, eight later exposures above nine earlier ones in both; all 34 native detector cores without SATURATED, DROPOUT, PERSISTENCE or OUTLIER bits; clean DO_NOT_USE-excluded 6+4 exposures per band; F162M aperture difference 236.179 ± 18.238 nJy; spectroscopic host z = 0.5067; host-conditional MΔ = −16.449 AB |
| GOODS-S 185 | 12 catalog filter/pointing pairs and 20,568 mutual associations; JADES DR2 host 103607, photometric z = 2.47 with [1.907, 2.514] interval; F277W aperture change 22.248 ± 1.154 nJy; host-conditional MΔ = −17.185; F277W independent-dither medians 50.61→71.49 nJy; relevant native peak DQ values zero |
| Blind GOODS-S | 2.3353/2.3307 arcmin² of tested area in F277W/F356W; 148/103 raw peaks, 105/80 morphology selected; 31 two-band coincidences, 28 flagged fixed-pixel defects; AB 27.5 recovery 212/294 (72.1%) and 246/324 (75.9%) |
| High-z NEXUS | 1,577 joined contrasts at 100 secure-z>5 host positions, **zero** jointly selected events; 5,743 control contrasts, three numerical crossings; F200W 1,867 contrasts / 120 positions, F444W 1,936 / 122; interpolation of recorded injections reproduces individual 90/50% depths 26.11/26.70 and 26.14/26.78 AB and joint 26.59/27.13 AB |
| MACS0416 and El Gordo | 0.190136 arcmin² and zero same-sign MACS paired events; 2,700 signed injection trials with interpolated 50% recovery AB 26.65; public CATS v4.1 critical-curve distance 0.410″ as deposited, with its distance-ratio computation checked independently; zero El Gordo joint events |

The repository's `Verify claims` GitHub Action also checks deposited numerical claims and links across the repository. The `Research checks` action initially failed to load two JWST test modules because its Python environment omitted pandas; the dependency was corrected and the full **Python 3.12 test run passed 292 tests**, including all six independent referee test methods and the existing immutable-catalog and NEXUS-stamp-cache tests. [Research checks run](https://github.com/dovvnloading/Open-Space/actions/runs/44020704940) · [Original full-tree verification and offline replay](https://github.com/dovvnloading/Open-Space/actions/runs/44019774408).

**Limits on this independent review:** This is an offline replay from the deposited compact inputs, *not* a second download and recalibration of every raw JWST FITS exposure. The original detector-level analyses relied on public remote FITS files, and some other READMEs document network-only replay steps. In the original offline verifier, three commands produced output drift after network attempts were refused; this does not demonstrate a scientific change, and these outputs have **not** been independently reconciled by raw-image reruns here. The GOODS-S README contained a literal `SITE PID` placeholder that failed automatic command replay; it has been replaced with executable examples. Numerical agreement with stored products is not proof that their original raw-photometry interpretation is correct.

## Input identity and receipts

Frozen STScI catalog lists record canonical HTTPS URLs, byte counts, source counts and SHA-256 digests in `repeat-catalog-followup/catalog_inputs.csv`, `repeat-catalogs/catalog_inputs.csv`, and the GOODS-S `8060_manifest.csv`/`12577_manifest.csv`. The independent tests verify receipt schema, digests and inventory consistency. COSMOS-022 history's `inventory_checksums.csv` SHA-256 digests were freshly recalculated against the three **deposited** source inventory CSV files and match. COSMOS detector `native_receipts.json` records 34 original exposures and byte-range hashes; GOODS detector, blind subtraction and MACS products have analogous byte-range receipts. NEXUS `f200w/mosaic_inputs.json` and `f444w/mosaic_inputs.json` list the original mosaics and SHA-256 values, with stamp-cache checks rejecting stale/changed inputs. The archival ESO MUSE spectrum records its original cube URL, raw byte count and SHA-256 in `cosmos-022-validation/muse_receipt.json`; the saved 1D spectrum and line fits are replayed here. The MACS CATS lens file hashes are in `macs0416/warhol_lens.json`. **Original remote FITS/MUSE bytes were not themselves re-downloaded in this review**, so their archived hashes are recorded and checked for format rather than asserted independently rehashed.

## TNS public-guest retry

At approximately **2026-11-17 05:00 UTC**, `tns_guest_retry.py` attempted **four** unauthenticated searches at 10-arcsecond radius: HTML guest search and `format=csv` download for each of COSMOS-022 (150.122033°, +2.323386208°) and GOODS-S 185 (53.188672°, −27.812554°).

| Position | Guest HTML | Guest CSV URL | Interpretation |
|---|---|---|---|
| COSMOS-022 | **HTTP 429** | **HTTP 200, HTML search page** (not a validated CSV) | No reliable automated cone result |
| GOODS-S 185 | **HTTP 200, HTML search page** | **HTTP 429** | No reliable automated cone result |

[GitHub Actions HTTP result and downloadable query receipt](https://github.com/dovvnloading/Open-Space/actions/runs/44020822897) (artifact `jwst-tns-guest-retry.json`). Queries are fully encoded by the source script, with timestamps, final URLs, HTTP status, media type and SHA-256 of any returned page. The service may throttle guests or require a session/API key; a successful HTML response is **not** a confirmed empty TNS search. Therefore the new official TNS attempt **did not clear either source as unpublished**. The pre-existing 305,064-row Fink TNS snapshot, published JWST catalogs, other broker mirrors and observing-team target lists show no positional counterpart **in the datasets inspected**, but cannot establish discovery priority.

## Referee scientific disposition

**COSMOS-022** and **GOODS-S 185** must remain **unclassified nuclear brightenings**. A low-level or obscured AGN fluctuating in brightness and an unresolved nuclear transient (including SN or TDE hypotheses) remain plausible. No luminosity in this PR is a fitted bolometric or transient peak magnitude: reported MΔ values are rest-band-equivalent, host-assumed *aperture differences*. The archival MUSE spectrum of COSMOS-022 predates the 2026 brightening and cannot classify it. GOODS-S 185's adjacent z = 5.2509 galaxy is not the candidate host.

**Required follow-up for each candidate:** (1) a classification-quality spectrum while it is elevated, looking for SN, broad/narrow AGN and TDE-like features against the host; (2) a **later JWST NIRCam epoch in matching filters**, with independent dithers and a consistent aperture/PSF prescription, to measure fading versus persistent or recurrent nuclear activity. No source in this PR reaches the [REMARKABLE] or discovery bar; do not add a remarkable-finding entry.
