# Native W2 tests of the literature-unmatched wide-area trajectory alerts

These sources passed the six-epoch catalogue triage. None of the five targets below produced a repeated positive motion residual at its reported flux and velocity in both scan seasons. The image tests compare new, held-out 2021–2024 NEOWISE visits to the early mission, with seven visits for J192049 and eight each for the other targets (312 native exposures total). Detection of a stationary mid-infrared source or a chance catalogue alignment is insufficient for assigning distance, temperature or a brown-dwarf type.

| CatWISE source | W2 Vega | W1−W2 | Moving / best alternative in expanded half-year sample | Native motion contrast, seasons 1 / 2 | Same position after adding a catalogue-flux mover, seasons 1 / 2 |
|---|---:|---:|---:|---:|---:|
| J010638.19+122908.3 | 15.274 | 2.850 | 2 / 1 (13 catalogues) | +0.89 / +2.24 | +4.93 / +5.31 |
| J213348.06+391040.9 | 14.917 | 3.091 | 3 / 3 (16 catalogues) | −0.06 / +1.28 | +6.28 / +9.31 |
| J192049.57+111222.6 | 15.042 | >3.875 | 10 / 11 (16 catalogues) | +0.51 / +0.03 | +3.67 / +2.23 |
| J152220.04−512432.2 | 15.134 | 2.018 | 5 / 2 (16 catalogues) | +0.68 / −0.97 | +4.84 / +1.98 |
| J175039.41−182034.4 | 15.135 | 2.395 | 4 / 5 (16 catalogues) | +0.04 / +0.20 | +5.57 / +3.22 |

J010638 has a weak positive difference in one season. It has only two expanded-sample centroid matches, with both in the same scan direction. J213348 has as many nominal matches along a rotated path as along its reported motion. The centroid cluster for J192049 also fits a static source. Its six-pass moving and stationary tracks each match all six visits when given the same freedom in reference position; the expanded sample matches 12 on each optimized track. J152220 crosses a crowded stellar field. Its best stationary location within 3″ matches six expanded-sample epochs and its optimized moving path matches seven. Only one of five nominal moving matches falls after CatWISE. The direct images show static features along the path and no positive moving signal in the later frames.

J175039 supplied only one scan direction in its primary triage. In the expanded sample the best alternative has five matches and the optimized stationary and moving models both have seven. Both direct-image seasons are consistent with zero motion contrast.

## Measurements and reproduction

In `assess.py`, the expanded W2 catalogue test chooses the highest-detection-count catalogue in **every** 182.625-day bin, as in the deep-field pilot. This can pick different scan directions and epochs from the six-window primary test. A match uses the same 1.5″ radius, S/N ≥4 and quality factor >0.85. The reported best alternative is the best of a fixed position and three rotated same-speed motions. Each tile directory contains detailed detection records, association scores and input URLs with digests. The algorithm samples the same observational archive, so the six-pass and expanded tallies are correlated.

For each source, `images.py` queries the final NEOWISE W2 frame inventory. It uses up to the first four visits before MJD 57400 and last four after MJD 59215, dividing them by pass parity. For each visit it selects up to eight W2 exposures by quality and proximity of the pointing centre. It projects the native pixels onto a 1.375″ grid with each exposure's WCS, scales by its photometric zero point and median-combines by epoch and scan season. The orange circle in each figure marks the early predicted position and cyan the late position.

The moving-source statistic compares the late-minus-early image with a positive-late minus positive-early Gaussian PSF profile. Its units are the median-subtracted contrast divided by 1.4826 × MAD over 52 nearby displaced model paths. The control apertures can overlap. These are empirical contrasts, without Gaussian significance calibration. To measure local sensitivity the code injects a mover of the target's catalogue W2 and motion at its actual position in **every native frame**, reprojects, combines and measures again. The model assumes constant brightness and a circular 6.4″ W2 PSF. An injected score near two indicates limited rejecting power in that season. These local tests characterize the specific paths and fluxes supplied.

Use NumPy, SciPy, Matplotlib and Astropy. Run `python assess.py SOURCE` for the expanded epoch set and `python images.py SOURCE` for the native-image tests. `--catalog` overrides the default `../wide/selected.csv`. The image script imports the native FITS strip reader from the pilot J0352 image audit. Source directories preserve the query URL, native metadata, exposure URL/header/strip SHA-256 digests, image panel and numerical results. Where `frames.tbl` was omitted for space, the script retrieves the original query specified in `frames.url`. Use the source name exactly as spelled in the CSV (the southern name uses an ASCII minus there).

[NEOWISE image archive](https://registry.opendata.aws/wise-neowiser/) · [unTimely catalogue DOI](https://doi.org/10.26131/IRSA580)
