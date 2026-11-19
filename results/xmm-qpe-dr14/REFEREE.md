# Independent methodological referee — XMM QPE DR14, PR #205 (2026-11-19)

**DECISION: KEEP DRAFT, DO NOT MERGE. No validated novel QPE; no population null; no [REMARKABLE].**

## Evidence I accept

1. Positive-control 0.2–2 keV background-subtracted EPIC-pn source curves recover the recurrent bright-peak pattern in all **six known-QPE active epochs**, with **29/29 published bright-event peak-count overlap**, but that is in-sample and not blind timestamp-by-timestamp sensitivity. RX J1301.9+2747 separately matches **8/10** fitted published times within 1200s. See outputs/CONTROL_RECOVERY.md and the unchanged control_thresholds.json.
2. The entire 33,488-entry soft/long DR14 X-ray parent was crossmatched against SDSS DR18 clean PhotoPrimary galaxy centers, 419/419 valid SQL chunks: **2,510 unique X-ray detections** with 2,519 possible galaxy counterparts within 2 arcsec. This is a reproducible **SDSS-footprint subset**, not a complete Legacy/SDSS galaxy survey.
3. The 4XMM-DR14 release was verified against the published FITS gzip MD5 **47c6bd113cb8510db418665472c81b48** following 60 HTTP 206 byte-range chunks. A 12-character PN_FLAG (F/T/-) implementation bug initially produced a false zero, was explicitly identified from the 4XMM quality-flag definitions and fixed. The corrected optical, EP_EXTENT, PN_ONTIME, PN flag and PN photon-count gate freezes **1,108 provisional detections / 897 X-ray sources / 742 XMM observations**. See outputs/archive_pn_point_freeze.json.
4. Statistical diagnostics on the six known positive controls: **240/324** synthetic paired rate injections matched both input peaks, **210/210** time permutations of eruption-containing control curves flagged recurrent patterns, **0/6** scaled source-background proxy flags. The time permutations have QPE outbursts, so **they are not null backgrounds**. The BACK proxy is **not an independent background region**. No survey-wide false-positive probability, FDR or independent completeness results from these numbers.
5. Three preselected observations (0872392901, 0306630101, 0200730401), six prospective pn detection rows, were requested via the existing PPS archive interface with 5-min/observation and 3GiB byte ceilings. **All three XSA connections closed without data; zero new-source pn light curves were analyzed.** Run receipts in outputs/archive_pilot_batch_01.json. The final safe six-control replay succeeded and keeps the gate false.


### Subsequent SDSS contaminant audit

After the 1,108-detection freeze, an independently recorded SDSS DR18 SpecObjAll/PhotoPrimary run completed without missing SQL batches: 205 unique spectroscopic host positions (419 spectroscopic rows), 119 detections classified as QSO/STAR/recognized BL Lac or coincident clean SDSS stars excluded, **989 surviving provisional detections / 795 unique SRCIDs / 688 obsIDs**. There were **zero** clean SDSS nearby stars within 1 arcsec, but **Gaia and full 4LAC/SIMBAD exclusions remain unverified**. A targeted connection probe again failed on the three pilot obsIDs and known control GSN 069, so archive PN access remains unavailable. [Receipt](outputs/archive_sdss_stellar_blazar_audit.json), [remaining subcohort](outputs/archive_sdss_after_exclusions.csv). No change to rejection disposition or validated new-QPE count.

## Blocking issues

- **No archive screen actually completed:** 0 of 742 provisional eligible pn observations successfully checked (three download failures are missing data, not non-eruptions); 0 vetted source flares. Therefore neither an archive-wide null nor a meaningful zero-candidate claim can be made.
- **Contaminants and host denominator unresolved:** the SDSS portion excludes much of Legacy Surveys footprint; Gaia stars, spectroscopic host-redshift/centroid astrometry and known blazars/AGN are not independently excluded. A morphological galaxy near an X-ray position need not be a confirmed nuclear QPE host.
- **Specificity and sensitivity unresolved:** no separately extracted full-field particle background or FBKTSR GTIs, source/field flare subtraction, positive stellar/TDE/AGN negative controls, injection of Poisson event lists, trial-count correction or held-out timestamp benchmark. A QPE-bearing curve time-shuffle is not a false-positive test.
- **Instrumental and external validation absent:** no SAS epatplot pile-up test or independent MOS confirmation, softening/temperature evolution, eRASS1/Swift/Chandra multi-epoch trend, literature/TNS/ATel QPE-vs-TDE-vs-AGN novelty audit for any new light-curve flag.

**Independent disposition:** refuse readiness and merge commit. Keep PR #205 draft until all 742 SDSS-selected pn observations and the LS sky outside SDSS are auditable, the background/stars/blazars and injection/FDR gates pass, and each independently flagged candidate is vetted. There are no new verified QPEs to label [REMARKABLE]. **Merge SHA: NONE.**

Detailed trial counts, source manifests, GitHub Actions run IDs and transport errors: outputs/ARCHIVE_SCREEN_AUDIT.md.
