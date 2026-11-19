# ACT depth-1 150 GHz pixel-search preregistration — before science FITS access

**Registered 2026-11-18 (America/New_York), after PR #192 was already merged.** This is a *post-merge branch experiment*, not a revision of or addition to the already-merged PR. PR #192's inventory-only conclusion remains intact. Do not count this document as a completed search.

## Immutable sampling rule

Source is the already committed gzip catalog `eligible-act-depth1.jsonl.gz`, whose compressed SHA-256 is `b6a0c020e3b03bc23aac66536ae1b31cb8889cbc3295c1d2acf16f8de3952503`. Test **all** rows with `tube_slot == "pa4"`, `frequency == "f150"`, and `"2018-09-10T00:00:00" <= start_time < "2018-09-13T00:00:00"` (catalog naive datetime convention pending FITS verification). No top-N truncation, selection by SNR, or cherry-picking; record count, map IDs, scan times and failures regardless of availability. One pre-fixed ROI: **ICRS RA 301.9952 deg, Dec +16.1652 deg, 0.30 deg great-circle radius**, at native 0.5-arcmin sampling, excluding incomplete PSF support / zero inverse variance / scan edges. This location deliberately tests the published 2018-09-11 Naess et al. 2021 transient, so it is a **targeted positive-control field** and not a statistically unbiased random survey field.

Reference candidates: *all* pa4/f150 frozen catalog rows with `"2018-08-21T00:00:00" <= start_time < "2018-10-04T00:00:00"`, **excluding every test-window map**. Require >=3 distinct other-night maps contributing valid science pixels to each reference/coadd location. Never use ACT DR6 or DR5 all-season coadds, nor a reference containing a test observation. For multiple test maps overlapping in time, treat correlated observations as one event cluster.

### Declared detection and falsification protocol

1. Retrieve authentic FITS primary T science, ivar, and mean-time pixels via HTTP Range (206 required) or bounded full download, recording per-product origin, HTTP metadata, byte-range receipts and FITS headers. Header-only reads do not count as science pixels; byte-range receipts are **not** full-file SHA-256.
2. Reproject each reference pixel to a fixed WCS with documented interpolation and ivar variance propagation; keep test scans entirely out of the reference. Preserve all failed reads.
3. Use measured ACT pa4/f150 optical/transfer-function beam (not an assumed Gaussian). Calibrate filtering with empirical off-source scan-noise, ≥3-night reference, independent split-epoch jackknife, sign-flip tail, and injected transient peaks. Screen at local SNR >=7 before empirical trials correction; any excess also needs a second independent array/visit and external vetting.
4. Predeclare injections into *actual* observed cutouts at 6, 8, 10, 15 sigma formal native white-noise levels, on a deterministic 5x5 valid sky grid, with random seed 192 and injected timestamps restricted to sampled observation times. A synthetic-only test does not establish empirical completeness.
5. Cross-check all survivors against published ACT transient tables, JPL Horizons ephemerides (with real pixel-level UTC timestamps), catalogs of variable radio/AGN sources, Gaia, VLASS/RACS, Swift/GCN and TNS. If any mandatory service fails, flag candidate **unvetted**; never call it a discovery.
6. Effective area-time at each flux/duration threshold requires per-pixel timestamp/visibility, real sensitivity and injection survival. No rate limit if any of these are unavailable. A published control recovery is not novel and cannot be relabelled as `[REMARKABLE]`.

This three-day and 0.30-degree bounded field was intentionally chosen before downloading the first science pixel. It tests pipeline feasibility, **not** the full ACT sky or 27,512-map inventory; extrapolation is forbidden. Report zero **usable** scans explicitly if products fail. The independent referee must treat incomplete controls or an unavailable measured beam as a stop, not as a null result.

References: [ACT public depth-1 documentation](https://phy-act1.princeton.edu/public/snaess/actpol/depth1/), [ACT published single-observation search](https://doi.org/10.3847/1538-4357/adcae0) and [Naess et al. 2021](https://arxiv.org/abs/2012.14347).
