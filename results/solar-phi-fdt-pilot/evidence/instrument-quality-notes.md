# SO/PHI FDT quality and timing metadata — 2026-11-17

Scope: instrument-team HTML, FITS headers and metadata tables only. No primary science image arrays or quicklook thumbnails were viewed. The [instrument team’s FDT table](https://www.uv.es/jublanro/phidata_fdt.html) was fetched successfully (6,469,587 bytes; SHA-256 `380b1b70e42d989344ce93441b893b83cbea0e5a7ba9459459bfe2c5a94575a8`; HTTP HEAD 200, Last-Modified 2026-10-21 14:37:59 GMT). The linked [quality-note table](https://www.uv.es/jublanro/phidata_fdt_notes.html?varnote=648130506) was fetched successfully (924,634 bytes; SHA-256 `99e06b11f76363c31d40cbdfa6d03466c56da6945fbcd44c2fc36bcfd470865b`; HEAD 200, Last-Modified 2026-10-21 13:08:10 GMT). Exact metadata receipt and extracted row strings are in `receipt.json`; HTML is retained in this directory.

The note-table description cell for each of the four exact DIDs says the product is an onboard average of three datasets, with lower noise. The release documentation independently says observations after 2025-02-22 generally use three on-board averages; this supports the four products despite the missing `ACCRAWIT` FITS keyword. This averaging is distinct from the `PHI_ACCU_accumulations=4` acquisition setting and from the 360 per-image time records; none makes 360 independent magnetograms.

| Product DID / observation start | Description-cell evidence | Interpretation |
| --- | --- | --- |
| 648130504 / Aug 13 10:25:09 | Three-dataset onboard average | No additional issue in this row’s description cell. |
| 648130506 / Aug 13 16:25:09 | Three-dataset average; ring-like pattern and stripes in Bmag | Relevant limitation: `B_LOS` is derived from Bmag and inclination; the local effect is not quantified. |
| 648150506 / Aug 15 16:25:09 | Three-dataset average | No additional issue in this row’s description cell. |
| 648150507 / Aug 15 19:25:09 | Three-dataset average | No additional issue in this row’s description cell. |

All four note cells also contain repeated strings including “Artefacts in azimuth and Bmag” and “Slight defocus.” The text visibly accumulates and repeats items across successive entries, whereas the separate description cells above are more specific. It is unsafe to assign every repeated note string independently to each exact exposure; absence from a description cell is not proof of pristine data. The main table links a group row for each date, and its hover text is group-level. Consult the exact note-table row and retain this attachment ambiguity. The team’s [FDT v03 release notes](https://www.mps.mpg.de/solar-physics/solar-orbiter-phi/data-releases) additionally describe dust and residual calibration artifacts, and the [data description](https://www.mps.mpg.de/solar-physics/solar-orbiter-phi/data) states that `blos` is derived from field strength and inclination.

All four primary headers report `ACCSTAT1=0`, `ACCSTAT2=64`, described only as “Accumulator Status.” Neither an `ACCRAWIT` header nor a matching table column was found in any extension. Searches of primary release documentation and exact-key web queries did not reveal the meaning of `ACCSTAT2=64`; do not interpret 64 as a quality pass or failure. The missing `ACCRAWIT` is a metadata inconsistency with the general release note, not evidence that no onboard average took place. A standard-L2 development run can be explicitly conditional on this unresolved status and documented artifact; finite pixels alone do not certify inversion quality.

All four headers use UTC and label `DATE-BEG`/`DATE-OBS` as start, `DATE-AVG` as average, `DATE-END` as end, and both `DATE_EAR` and `DATE_SUN` as corrected *start* times. These differences are from the respective `DATE-BEG`, not the midpoint:

| Scan start UTC | DATE_EAR − DATE-BEG (s) | DATE_SUN − DATE-BEG (s) |
| --- | ---: | ---: |
| Aug 13 10:25 | +345.571 | −160.008 |
| Aug 13 16:25 | +346.604 | −158.952 |
| Aug 15 16:25 | +353.477 | −151.899 |
| Aug 15 19:25 | +353.815 | −151.549 |

Use the appropriate scan midpoint offset if matching midpoint observations, retain observer ephemerides at observation time, and check HMI’s `T_REC` convention, UTC/TAI, and temporal response before comparing exposure windows. The [primary FDT–HMI comparison](https://doi.org/10.1051/0004-6361/202349096) accounts for light-travel delay; [HRT–HMI work](https://doi.org/10.1051/0004-6361/202245830) also accounts for UTC versus TAI. Four-minute spectral scans and HMI’s multi-sample temporal averaging prevent treating either as instantaneous. Geometric overlap and different line-of-sight projections remain separate constraints.
