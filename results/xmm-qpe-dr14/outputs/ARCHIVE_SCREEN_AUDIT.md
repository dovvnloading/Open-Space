# 4XMM-DR14 gated QPE host and pn archive-screen audit — 2026-11-19

**Disposition:** incomplete SDSS host-selection and failed PPS pn-screen pilot. PR #205 stays **draft**. No merge SHA, no validated new QPE, no archive null and no [REMARKABLE].

## Exact denominators

| Successive gate | Count | Interpretation |
| --- | ---: | --- |
| Full 4XMM-DR14 detections | 1,035,832 | Release 1.0 published FITS gzip MD5 47c6bd113cb8510db418665472c81b48 |
| Frozen soft/long X-ray prefilter | **33,488** | Existing TSERIES, SUM_FLAG<=1, EP_8_CTS>=250, camera ontime>=21600s, EP_HR3< -0.3, EP_HR4< -0.4 |
| SDSS DR18 clean primary galaxy at <=2 arcsec | **2,510** detections; **2,519** counterpart pairs | 419/419 SQL chunks completed for all 33,488 X-ray rows; 9 ambiguous counterpart detections |
| SDSS plus 4XMM pointlike and long PN quality gate | **1,108** provisional detections; **897** unique SRCIDs; **742** obsIDs | Unique host, EP_EXTENT<=6, EP_EXTENT_ML<4, PN_ONTIME>=21600, PN_FLAG serious bits 7/8/10/11/12 unset, PN_8_CTS>=250 |
| PPS new-source pilot attempted | **3 obsIDs**, **6** eligible source rows | Prespecified highest catalogue-ranked low-pileup subset; all three XSA connections aborted |
| **Successful pn observational screens** | **0 of 742 obsIDs**, **0 of 1,108 detections** | Failed downloads are missing observations, **not absence of QPEs** |
| Novel recurrent flags measured | **0 measured (sample unavailable)** | No per-source result; cannot call zero archive candidates |
| Validated new QPEs / remarkable labels | **0 / 0** | No discovery evidence, no novelty classification attempted |

Host crossmatch QA: X-ray parent SHA-256 c2793211c74f9a85992253fba76743a94bcdf2cffb2137f29a3037907d4e56a5; SDSS pair file SHA-256 b1495de1b3f28a7e3a7b21f73df30f96d43607a77db0231d445e1280eae0481f; final provisional point-pn manifest SHA-256 93475a72d00393c236b36065a60ef6fdb6f8c066c9e945c7e543ec47a8321357. All three input stages were recorded with source-code and API receipts.


## Subsequent SDSS-only stellar / broadline AGN rejection (not Gaia-complete)

A further bounded SDSS DR18 audit of the 1,108 initially point/pn-qualified detections completed all SpecObjAll and PhotoPrimary near-star batches. It identified 205 unique SDSS spectroscopic hosts (419 spectrum rows), **119** detection rows with class QSO or STAR or recognized BL Lac/blazar subclass to exclude, and **0** clean morphological SDSS stars within 1 arcsec. The resulting second-stage provisional host list has **989 detection rows, 795 unique X-ray SRCIDs and 688 distinct XMM obsIDs**. The baseline freeze of 1,108 is preserved unchanged for provenance; the smaller file is outputs/archive_sdss_after_exclusions.csv and complete query receipts are in outputs/archive_sdss_stellar_blazar_audit.json. This is an SDSS-only spectroscopic/stellar exclusion: host galaxies without spectra are unclassified, and proper-motion/parallax Gaia stars and SIMBAD/4LAC non-SDSS blazars have NOT been systematically excluded.

Four XSA no-download probes (the three pre-registered pilot observations and known-QPE control GSN 069) again failed with RemoteDisconnected. Therefore the **989 remaining provisional detections / 688 observations are still not screened at the pn curve level**, and the FDR and candidate vetting gates remain closed.

## Reproducibility and repairs

- NOIRLab LS DR10 TAP cone and upload timed out or returned HTTP 502; its anonymous QueryManager query endpoint returned HTTP 401. Beta endpoint certificate validation failed. A documented SDSS DR18 SQL alternative succeeded via SkyServer: all 419 chunks completed without missing IDs. Only SDSS-footprint host completeness is addressed; do not extrapolate to the whole LS sky.
- A direct 959 MB HEASARC 4XMM gzip stream exceeded its 18-minute cap, and the second hostname failed DNS. A later HTTP 206 probe proved reliable range support and a 1,005,466,736-byte compressed file. The range downloader retrieved 60 bounded 16 MiB chunks with 8 workers, assembled the file and independently checked published MD5.
- The first point-pn filter incorrectly cast the 12-character PN_FLAG into a number, excluding 1,645 sources. This zero was explicitly **invalid**. The 4XMM official catalogue quality flags are characters F/T/-; the corrected filter applies serious-artifact bits 7, 8, 10, 11 and 12 and retains **1,108** provisional eligible detections. See [official 4XMM flag definitions](https://xmmssc.irap.omp.eu/Catalogue/4XMM-DR14/col_flags.html).
- The pipeline restored the frozen six-positive-control no-network replay after the pilot; archive_screen_permitted remains **false**. Every runner job had bounded work and used output-only commits excluded by the path-filtered push workflow.

## Surrogate statistical tests and why they are not FPR

The fixed positive-control finder was applied to previously saved, source-position-verified NET EPIC-pn 0.2–2 keV 300-second curves; thresholds were not retuned. All six control epochs recover recurrent QPE activity. The 29/29 count overlap is **in-sample count-only**, not one-to-one timing recovery; RX has 8/10 fitted individual peaks within +/-1200s.

- Quiet-resampled positive-control baseline plus two Gaussian rate pulses, 450/850/1300-second widths and nominal 4, 8, 16-error amplitudes: **240/324** tests detected both injected peaks within 900s. Artificial rate injections on known controls are not event-list Poisson injections into the unknown survey population. No selection-function inference.
- Full permutations of good-bin *eruption-containing* positive-control curves: **210/210** emitted recurring flags. This is not a blank-sky/null experiment: genuine outburst bins remain in the rate distribution. No empirical archive false-positive rate or false discovery proportion can be estimated from these trials.
- Scaled source-extraction BACK1-3 substituted as a proxy series: **0/6** recurrent flags, but this is not independently extracted region photometry or FBKTSR/high-energy particle GTI. No particle-background screen and no stellar, TDE, persistent or flaring AGN negative sample.
- There is no multiple-testing correction across the 1,108-detection parent, no verified survey completeness and no non-detection exposure denominator. An archive occurrence-rate limit would be unjustified.

## PPS test and uncompleted vetting

Three pre-fixed (before PPS inspection) low-catalogue-pileup high-rank observations were attempted: **0872392901**, **0306630101**, **0200730401**, covering **six** provisional source rows in the frozen host list. All returned connection-aborted RemoteDisconnected from ESA PPS before any FITS product. No independent source FITS position, particle background, pile-up, stellar flare rejection, multi-epoch eRASS1/Swift/Chandra or literature/TNS/ATel QPE-vs-TDE-vs-AGN classification was possible. Any observed zero of flags is only a missing-data count.

## Independent referee decision and actions history

**Keep PR #205 draft, do not merge** until (a) every eligible PN observation is screened and the full LS/SDSS sky census/selection function is completed, (b) Gaia/star/blazar and other AGN contaminants are excluded, (c) off-source background, particle GTIs and pn epatplot plus spectral evolution pass, (d) genuine non-QPE negative controls and event-list injections calibrate FDR and completeness, and (e) candidate-specific eRASS1/Swift/Chandra, TNS/ATel and literature checks are complete. None is fully satisfied. No [REMARKABLE] claim; no docs/REMARKABLE_FINDINGS.md entry.

Evidence: [SDSS crossmatch run](https://github.com/dovvnloading/Open-Space/actions/runs/44493308016) and archive_sdss_receipt_00000_33488.json; [initial catalogue timeout](https://github.com/dovvnloading/Open-Space/actions/runs/44493881640); [shuffled and injected controls, HTTP range probe](https://github.com/dovvnloading/Open-Space/actions/runs/44497242211) and archive_statistical_checks.json; [first range download with invalid PN_FLAG zero](https://github.com/dovvnloading/Open-Space/actions/runs/44497525830); [corrected PN_FLAG freeze](https://github.com/dovvnloading/Open-Space/actions/runs/44498129634) and archive_pn_point_freeze.json / archive_pn_point_freeze.csv; [three failed PPS requests](https://github.com/dovvnloading/Open-Space/actions/runs/44499224827) and archive_pilot_batch_01.json; [final safe replay](https://github.com/dovvnloading/Open-Space/actions/runs/44499571847) and control_gate_replay.json. Control recovery: CONTROL_RECOVERY.md.
