# LoDeSS image validation and release provenance

## Reproducing the measurements

Use Python 3.12 with NumPy, SciPy, pandas, Astropy, lxml and Matplotlib. Commands below run from this directory. Set `DATA` to a local data directory and `MATCHES` to the directory containing the comparison pipeline's `match_lotss.csv.gz`, `match_nvss.csv.gz` and `match_vlssr.csv.gz`. The images selected here include each of the 13 good fields and the mediocre-quality repeat pointings P136+30, P210+80 and P264+86. Allow approximately 7 GB for the FITS inputs. Successful downloads use temporary files and replace completed products atomically.

```bash
DATA=./data
MATCHES=../comparison/results
python acquire_lodess.py "$DATA" --images P36+35 P57+47 P128+36 P130+60 P136+36 P139+79 P151+49 P171+49 P181+55 P219+50 P238+32 P246+74 P249+80 P136+30 P210+80 P264+86
python map_fields.py "$DATA"
python local_controls.py "$DATA" "$MATCHES" --targets 9927 2114 8029 8092 13292 6907 2928 7012 18228 8889 11986 7001 13277 4926 19498 11686
python inspect_rows.py "$DATA" 9927 2114 8029 8092 13292 6907 2928 7012 18228 8889 11986 7001 13277 4926 19498 11686
python resolution_check.py "$DATA" "$MATCHES" 6907 2928
python target_images.py "$DATA"
python channel_checks.py "$DATA" "$DATA/channels"
```

The archive manifest gives the byte count, SHA-256 and direct URL for both official FITS tables, all 40 Stokes I WCS headers, and 16 complete Stokes I fields. The official images combine the full 15–30 MHz bandwidth. Of the 40 headers, 31 give a frequency-axis reference value of 21.971 MHz, and nine give 22.264 MHz. The survey paper reports an effective central frequency of 23 MHz. These values describe a single continuum flux measurement. [DR1 release](https://lofar-surveys.org/lodess.html) · [Cordun et al. 2026](https://arxiv.org/html/2610.04698v1)

## Field identity and a release discrepancy

The source table includes 21,465 rows and the Gaussian table 33,890. Neither provides a field name per row. The source identifier resets between pointings. We split the Gaussian table into 40 contiguous sequences by Gaussian ID and assigned each sequence by forward projection from its sky coordinates through each of the 40 actual field WCS headers. Each winning median pixel residual is below 0.000001 pixel. Every alternative field exceeds 2,325 pixels. We verified the 32 source blocks by matching the source and island identifiers of each single-component source to its parent's Gaussian catalog position within 1 arcsec. All 17,811 such sources matched in one field, with zero concordant alternative assignments. The component catalogue retains repeated observations. Sources in the published deduplicated catalogue can preferentially reflect the brighter detection in overlap regions.

Using Table 3 of Cordun et al., the source table contains 11,008 good, 8,822 mediocre, and 1,635 poor-field rows. The last group comprises 707 in P100+30, 464 in P17+40 and 464 in P317+27. This conflicts with the release webpage's blanket statement on poor-field exclusion. The row-to-field crosswalk in `evidence/source_fields.csv.gz` and block-level verification in `field_membership.csv` allow this check without guessing which field was closest.

## Photometry and local checks

`RadioImage.measure(ra_deg, dec_deg)` measures signed Jy at a fixed position by projecting the background-subtracted pixels onto a Gaussian with the image's restoring FWHM. Its local statistical error is 1.4826 × MAD of the same PSF estimate at up to 40 positions 6–12 restoring beams away. Reported aperture, residuals and peak displacement diagnose departures from the point-source model. Error estimates exclude the survey's additional 39% systematic recommended for compact, bright, isolated sources in good fields. Fainter sources and field edges can suffer larger losses. A 60-arcsec-radius aperture in the separate resolution check tests sensitivity to diffuse emission and neighbouring flux. An aperture's pixel noise is not its integrated uncertainty.

In 1,008 compact, isolated good-field references, a straight spectrum through 73.8, 144 and 1,400 MHz generally overpredicts the measured 23 MHz flux. The reference selection, regression and local radii are in `local_controls.py`. These ratios combine spectral curvature, calibration, selection and epoch differences. They are insufficient to derive a source-specific scale factor. The 6- and 20-arcsec LoTSS products share visibilities, and repeated LoDeSS fields share an instrument and calibration pipeline.

## J112022.1+502936

All four WCS-preserving FITS stamps and their image mosaic are in `evidence/j112022/`. In the best pointing the fixed PSF returned **0.0868 ± 0.0199 Jy** at the LoLSS position, before flux-scale systematics. Two additional pointings returned 0.177 ± 0.102 Jy and 0.002 ± 0.131 Jy at off-axis distances of 8.61° and 11.78°. Neither gives a well-measured repeat flux. A 60-arcsec-radius aperture measured 0.069 Jy at 23 MHz and 0.387 Jy at 54 MHz after convolution to 45 arcsec. At 144 MHz the same aperture and beam measured 0.378 Jy using the 6-arcsec input and 0.337 Jy using the 20-arcsec input. The 54- and 144-MHz images show one dominant compact component. The nearest qualifying reference has a factor 7.20 predicted/measured ratio at 23 MHz. The seven references within two degrees have a median of 2.26 (interquartile range 2.09–2.88). Angular smoothing alone preserves the brighter-frequency aperture fluxes, while the absolute flux lost locally at 23 MHz remains uncertain. [LoLSS DR1](https://lofar-surveys.org/lolss.html) · [LoTSS DR3](https://lofar-surveys.org/dr3.html)

## 44–64 MHz channel checks

The published LoLSS DR1 page and release paper provide six channel source and Gaussian catalogues. Their public field image links expose the full-band Stokes I and V mosaics. The ASTRON TAP table `lolss.mosaic` returned only the preliminary continuum mosaic on 2026-11-16. The channel FITS catalogue headers identify the original image naming pattern. The corresponding source pointing `p170+50-mosaic-chan00.fits` returned HTTP 404 in the published Stokes I directory. We could not remeasure six independent fixed-aperture channel fluxes from these public products. [LoLSS DR1 products](https://lofar-surveys.org/lolss.html) · [Release paper](https://arxiv.org/abs/2301.12724)

Instead, `channel_checks.py` compares two measurements from the same catalogue fits. Its unweighted log-log slopes for J112022 are +0.638 for integrated flux and +1.109 for peak brightness. The fitted major axis decreases from 18.40 to 16.43 arcsec, and the source classification changes to multi-component at 48 MHz alone. Among 58 compact, isolated high-S/N comparison sources within 2°, three have integrated slopes at least as positive. Their median difference between peak and integrated slopes is +0.263, with a 16th–84th percentile range of −0.004 to +0.486. The target's difference (+0.471) lies in that range. The statistics are descriptive for detected, quality-selected sources and do not represent a preregistered discovery significance. The catalogue comparison cannot substitute for a common-aperture remeasurement of the six channel images.

## Expanded LoLSS pilot

The expanded search lowered the integrated 54-MHz threshold to 60 mJy and allowed an integrated-to-peak ratio below 1.65. It retained full-band single-component sources with no neighbour within 60 arcsec, within 3.5° of a good LoDeSS field centre. Each of six channel matches had to lie within 5 arcsec and have peak S/N at least eight. Of 780 parent sources, 743 passed. We ranked the spectral slope after subtracting the published indicative +0.27 bias and followed up the leading 24. This ranking has no correction for the 743 trials or source-specific calibration errors.

J115948.8+582020 ranked first. Multiple papers already identified it as a peaked compact steep-spectrum source; the fitted turnover is about 270 MHz. LoDeSS row 8833 lies 1.49 arcsec away. Its 23-MHz forced flux, 264 ± 25 mJy statistical, appears to rise again below the 44-MHz point of 62.6 ± 10.0 mJy. A 39% fractional uncertainty gives a combined 23-MHz uncertainty of 106 mJy. The only repeat observation gives −9 ± 68 mJy at 10.17° off-axis, and the nearby lobes of the unrelated X-shaped radio galaxy CGCG 292-057 complicate separation. The 23-MHz aperture sum decreases from 250 to 137 mJy when its radius increases from 30 to 90 arcsec. The public image supports a follow-up target, while the current measurements do not independently establish a secondary emission component. [Misra et al. 2023](https://doi.org/10.1093/mnras/stad1533)

J144527.0+500921 yielded 19.8 ± 12.2 mJy at 23 MHz. Its six-channel spectral index was +1.34 ± 0.43 before bias correction. The residual score of 2.51 across this search, weak 23-MHz measurement, and lack of a VLSSr detection leave it a lower-priority lead. The 743 ranks, 24 multi-survey crossmatches, field image and aperture sweep are in `evidence/expanded`.

```bash
python ../expanded_sample.py "$DATA" "$DATA/channels" /path/to/comparison/raw "$DATA/expanded"
python expanded_images.py "$DATA" "$DATA/expanded"
```

## J100005.91+812702.7 and overlap deduplication

The released source catalogue includes row 9460, 0.83 arcsec from the LoTSS position. It comes from the brighter P210+80 detection, a mediocre pointing 9.22° off axis. Filtering the published source table for good fields first hides the positional match. All three overlapping observations are present in the Gaussian catalogue within 2.1 arcsec. The good P139+79 counterpart is at Gaussian row 8034. The coordinates, identifiers, integrated and peak fluxes, fit errors and field names are in `evidence/J100005/components.csv`. The single retained source matches the largest measured integrated flux, consistent with the published overlap rule.

| Field | UTC date | Field quality | Off-axis distance | Forced flux (Jy), statistical error | 60-arcsec-radius aperture (Jy) |
|---|---|---|---:|---:|---:|
| P139+79 | 2021-12-02 | Good | 3.02° | 0.671 ± 0.025 | 0.701 |
| P210+80 | 2021-12-30 | Mediocre | 9.22° | 1.059 ± 0.054 | 0.808 |
| P264+86 | 2022-05-03 | Mediocre | 10.48° | 0.435 ± 0.063 | 0.135 |

Each image has a peak within 3.4 arcsec of the LoTSS position. The first pointing also has consistent 45- and 60-arcsec-radius apertures (0.663 and 0.701 Jy), with a small PSF-fit residual. The very far off-axis observations show larger PSF residuals and sensitivity to aperture size. The three epochs establish positional repeatability. Amplitude differences are confounded with beam and calibration effects. The 144-MHz images have one dominant compact component. Identical 60-arcsec-radius apertures yield 3.125 Jy from the 6-arcsec image and 3.031 Jy from the 20-arcsec image, a 3.1% difference relative to the latter.

```bash
python J100005_images.py "$DATA" "$DATA/J100005"
```
