# Independent-referee pass: DJA/MAST NIRSpec extreme-z preflight

**Date:** 2026-11-18 UTC/EST boundary. **Decision:** ACCEPT AS A PROVISIONAL DATASET-SCOPE / METHODS DEPOSIT ONLY; **REJECT any detection, any survey-complete null and any posterior-based physical conclusion**. This is a critical second pass over the deposited method and source ledger, not an independent observer or peer review. No spectrum bytes were available for an independent recomputation.

## Claim audit

| Test | Independent referee disposition | Why |
| --- | --- | --- |
| DJA version, shard counts, source URLs and publisher MD5 | **Scope frozen, raw verification pending** | Nine published files account for 80,367 *spectrum records*, not 80,367 independent galaxies; MD5 is publisher-reported, not observed |
| Newly public MAST selection | **Dashboard snapshot only** | The grouped observing-program rows do not establish individual public science product IDs, access rights, instrument exposure mode, source identity, or downloadable 2D spectra; the MAST API selection has not run |
| Per-spectrum independent fit | **NOT EXECUTED** | No complete `frozen_row_manifest.csv`, raw source SHA-256s, or generated `screen_rows.csv`. The implementation is reproducible in principle, but no astrophysical fit result is presented |
| Lyman-break + damping wing likelihood | **Exploratory method, not validated** | Only a step and simple exponential wing are approximated; line-spread function, IGM sight-line distribution, correlated extraction errors, sky residuals, throughput, photometric calibration and wavelength solution remain uncalibrated |
| Low-z interloper odds and redshift posterior | **NOT REPORTABLE** | Three low-z proxies and one artificial dwarf template do not span realistic brown-dwarf atmosphere models or emission-line galaxies; high-z class prior 0.02 is arbitrary; class evidence has a flat-amplitude proxy prior. Model posterior numbers, if calculated, must be labelled proxy and must **not** be quoted as calibrated `P(low-z|D)` |
| Break visible in 2D, contamination, bad pixels, slit overlaps | **NOT EXECUTED** | 1D combined-spectrum table is insufficient; native `s2d`/`cal`, masks, shutter maps and separated exposures are required |
| NIRCam photometry + catalog/alias literature priority | **NOT EXECUTED** | No NIRCam flux/limit files, calibrated SED convolution or object-specific ADS/arXiv/DJA/JADES/CEERS/UNCOVER/MoM/GLASS search receipts |
| Remarkable-finding gate | **PASS: no remarkable claim warranted** | No validated discovery-level evidence; do not add `[REMARKABLE]` or edit the notable log |
| PR #194 separation | **PASS** | The submitted files are limited to this independent NIRSpec campaign and its catalog/CI entry |

## Mandatory corrections before a scientific claim
1. Obtain every actual DJA v4.4 table and matching FITS product; preserve SHA-256, FITS schema receipts and complete source-level manifest. A 2D inventory must identify native instrument/program, slit and per-exposure units independently of the combined 1D table.
2. Confirm extraction-column orientation/wavelength units on real products. Check FITS `wave`, `flux`, `full_err`, `valid`, `npix` rather than trusting a description; hold out entire exposures/programs before tuning thresholds. Test computational and calibration feasibility on a subset before claiming full scan.
3. Estimate proper pixel covariance/background systematics and fit physical high-z/low-z/stellar template grids with explicit priors and prior-sensitivity, including grism/grating throughput and line confusion. Calibrate simulated injected-break recovery and false ranks on disjoint controls. Quote posterior intervals only when physically and empirically justified.
4. For every candidate, view the trace and break in both 2D pixels and independent exposures, audit bad pixels and slit contamination, compare NIRCam blue nondetections and red detections with convolved SEDs, and record all source-specific literature queries and priority decisions. Do not reuse a combined exposure as independent confirmation.
5. Keep the entire screened denominator, all ambiguous/rejected cases and blank-sky/contaminant controls, rather than reporting only a ranked subset. A null requires demonstrated sensitivity and coverage, not merely an empty candidate register.

## Released outcome

**Provisional archive and method setup; no newly unreported z > 12 galaxy established.** This is *not* evidence that none exist. A secure redshift, a calibrated probability for a low-z galaxy solution and a completeness limit are not available. Marking the PR ready/merging this documented preflight means accepting the reproducible scope and explicit uncertainty, **not** declaring the requested exhaustive search finished.
