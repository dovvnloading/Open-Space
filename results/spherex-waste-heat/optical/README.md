# Optical dimming in the Hephaistos fields

**Seven of the eight stars with ZTF DR24 r-band coverage have no measured night at least 0.20 mag fainter than baseline. D has one discrepant night, 2019 September 19.** B and C have no r-band points in the queried cones. All ten archive requests succeeded, including G and I. The catalog returned 3,934 measurements and the quality and matching cuts left 2,701 star-nights across eight targets. The per-field and per-night measurements are in [r_lightcurves.csv](r_lightcurves.csv) and [nightly_summary.csv](nightly_summary.csv).

## Completed measurements

| Hephaistos target | Archive points | Unique qualifying nights | First–last MJD | 95th-percentile fading (mag) | Largest nightly fading (mag) | Nights ≥0.20 mag |
|---|---:|---:|---|---:|---:|---:|
| A | 167 | 163 | 58219–60829 | 0.043 | 0.071 | 0 |
| B | 0 | 0 | — | — | — | — |
| C | 0 | 0 | — | — | — | — |
| D | 645 | 495 | 58268–60967 | 0.044 | 0.304 | 1 |
| E | 577 | 415 | 58337–60935 | 0.030 | 0.117 | 0 |
| F | 215 | 201 | 58370–60969 | 0.034 | 0.072 | 0 |
| G | 466 | 366 | 58288–60969 | 0.029 | 0.055 | 0 |
| H | 440 | 341 | 58280–60965 | 0.046 | 0.099 | 0 |
| I | 616 | 328 | 58217–60968 | 0.049 | 0.109 | 0 |
| J | 808 | 392 | 58217–60968 | 0.021 | 0.106 | 0 |

D measured r = 17.969 ± 0.028 in field 497 on MJD 58745.172. Five hours later field 1544 recorded r = 18.404 ± 0.044 and 18.640 ± 0.052 in two 30-second exposures separated by 39 seconds. That field's nightly median was 0.574 mag below its normal brightness. Combining the independent field medians for the same date gives 0.304 mag and a local score of 8.0. D returned to r = 17.928 the next observed night. Its infrared source has independently been spatially separated into a foreground M dwarf and a background obscured galaxy. [JWST/MIRI follow-up](https://arxiv.org/abs/2607.09460)

Both faint exposures pass the public catalog flags, but their ZTF PSF-fit χ values, 1.925 and 1.832, exceed all 101 other retained measurements of D in field 1544 (median 0.688, maximum 1.776). Their centroids also sit 0.50 and 0.73 arcsec from the propagated Gaia position, compared with a median 0.11 arcsec on the other dates in this field. These diagnostics motivate a check of the science pixels before assigning the fading to D.

## Measurement and replay

`ztf.py` propagates each Gaia DR3 position to J2022 and retrieves the r-band DR24 detections inside 2.5 arcsec. Each exposure must have zero catalog flags, magnitude uncertainty ≤0.10 mag, and separation ≤1.5 arcsec from Gaia propagated to that exposure. The two southern targets B and C returned header-only tables. The archive stores distinct object IDs for the same star observed in different fields or detector quadrants. [IRSA API](https://irsa.ipac.caltech.edu/docs/program_interface/ztf_lightcurve_api.html) · [DR24 release](https://irsa.ipac.caltech.edu/data/ZTF/docs/releases/dr24/ztf_release_notes_dr24.pdf)

`nightly.py` verifies the archived response SHA-256, groups observations by integer UTC MJD, and uses every object ID with at least five distinct nights. It computes a nightly median per ID, subtracts its time-series median, takes the median across contributing IDs for each date, and centers that star's combined time series. A positive residual means fainter. The local score divides a residual by the larger of the median reported exposure error on that date and 1.4826 times the median absolute deviation across the star's nights. The score describes the excursion and does not estimate its global false-alarm probability. Recalculate the table with `python nightly.py` from this directory. [Query receipts](query_receipts.json) · [Raw API responses](lightcurve_inputs.zip) · [Night-by-night measurements](nightly.csv)

## Pending checks and limits

The r-band archive selects detections. Rare eclipses, faint states below the detection threshold, events between the sampled times and dimming averaged out within a night remain possible. The recorded passbands and dates do not establish continuous optical stability or a survey-wide occultation occurrence rate. B and C have no optical constraint from these queries. The same-field science-image check for D is in progress.
