# KM3-230213A: source-association audit (claim-first, frozen plan)

**Initial claim (2026-11-18, before novel-catalog screening):** No unique astrophysical source is yet established for KM3-230213A. A positional match, including an appealing flare, is *not* a detection-level counterpart. This investigation first reproduces the published literature and then screens disjoint additional public catalogues; all resulting hypotheses require a sky/background and trial-corrected chance calculation. Initial status: **IN PROGRESS; no discovery claim**.

## Fixed observing unit and localization
One event: 2023-02-13, MJD(TBD from paper)=59988.0533299; ICRS/J2000 RA=94.3 deg, Dec=-7.8 deg; circular containment radii 50%=1.2 deg, 68%=1.5 deg, 90%=2.2 deg, 99%=3.0 deg (Nature 2025, https://doi.org/10.1038/s41586-024-08543-1). Cone searches use great-circle distances, not RA/Dec box. The 3.0-deg region is the pre-registered eligibility mask; report inner radii separately. Circular limits are published approximations; no event-level HEALPix posterior is assumed.

## Frozen search and falsification plan (registered before additional screening)
- **Known/positive controls:** nature paper's 17 blazar candidates (4FGL-DR4, eROSITA+NVSS+WISE, RFC VLBI, 5BZCAT), published radio results (ASKAP 944 MHz, VLASS; including UGCA 127, Hebe, Narcissus), known blazar radio/X-ray flares, cosmogenic models and IceCube non-detections. Classify ALL as **KNOWN**, not new.
- **Optical transients / TNS / ZTF / ATLAS:** pre-registered temporal selection of first detection, discovery, or peak in UTC 2023-01-14 through 2023-03-15 (±30 days), separately slower TDE/nuclear events 2022-08-17 through 2023-08-12 (±180 days). Require original survey coverage, not a catalogue result alone; TNS discovery dates are not transient onset dates.
- **Short transients:** GRB / GCN and CHIME FRB epochs within ±1 day, with a separately recorded ±1-hour timing gate; distinguish published localization footprints from point centroids and reject spatially uninformative events.
- **X-ray:** Swift/BAT, Swift/XRT, MAXI time-resolved events ±7 days, with ±30-day sensitivity check. Instrument footprints and archived non-detection coverage matter.
- **Flaring AGN:** Fermi LAT light-curve repository qualifying weekly bins within ±30 days; secondary ±180 days. OVRO/VLBA time-resolved 15-GHz radio observations same ±30/±180d windows. Use original flux and uncertainty, never merely a 'variable' label.
- **Late radio / mid-IR:** VLASS epoch differences, RACS/ASKAP, and NEOWISE W1/W2 light curves (±180-day IR window, respecting sparse cadence), with catalogue epoch, beam, flux-error and selection-bias audits. These datasets cannot establish 2023 simultaneity from coadds alone.
- **Persistent alternatives:** nearby star-forming galaxies, clusters, TDE hosts / AGN; query within 3 deg; record redshift/distance and baseline density, not a temporal association.
- **Statistics:** per-catalogue (a) explicit per-object angular separation and within-contour membership, (b) area/flux/temporal-filter-adjusted source-density prior from a cited parent sample, (c) expected random associations and a *family-wise* sky-shifted control with ≥999 valid directions where public completeness permits. Fix declination in RA-shifts; exclude physical overlap, unobserved patches, and state that Galactic/survey gradients can defeat uniform RA nulls. Registered families: optical, GRB, FRB, X-ray, gamma, radio, MIR, persistent (8); multiple windows/trials corrected by max-score per shift across all attempted tests. No p-value if denominator or exposure footprint is absent: record **not estimable**, not zero, and do not call a candidate significant.
- **Decision rule:** only independent post-selection validation, stable evidence under sky shifts/completeness and properly adjusted global p, plus a plausible neutrino production model, could support a discovery claim. Otherwise null/inconclusive. Never mark [REMARKABLE] solely for positional match.

## Claim status and audit trail
- Pre-registered hypotheses: as above, **frozen at initial branch commit**.
- Completed: repo instructions, docs/REMARKABLE_FINDINGS.md criterion, prior IceCube study scope; published event position/time and contour numbers.
- Pending: independently transcribe 17-source table; download/freeze raw survey queries and hashes; run positional/time controls, check completeness, test finalists against independent data; referee audit.
- Related but distinct: [IceCube alert crossmatch](../neutrino-transient-crossmatch/README.md), whose results are *not* re-used as KM3NeT associations.
- **Data receipts, scripts, evidence, provenance and status** will be linked below as each is committed. No p-value is asserted before the denominator is known.
