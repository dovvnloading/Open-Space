# Comparison surveys and reverse search

The matching unit is a zero-based row of the official LoDeSS source FITS table. Overlapping pointings can record the same physical object. Each survey match retains the nearest centroid, angular separation and a count of possible neighbours. Primary measurements come from NASA HEASARC's native catalogues and ASTRON's LoTSS DR3 TAP service. The table lists the number of LoDeSS rows with a match among all 21,465 rows.

| Survey | Frequency (MHz) | Search radius | Matched rows | Photometry |
|---|---:|---:|---:|---|
| Revised 8C | 38 | 105″ | 2,840 | Integrated. For unresolved sources the peak/SNR approximates local noise. Keep the extent and P/C flags. |
| VLSSr | 73.8 | 45″ | 14,556 | Peak Jy/beam and local noise. Compare to integrated values only for unresolved sources. Axes in this archive are the convolved sizes. |
| GLEAM EGC v2 | 76–227 in 20 subbands | 60″ | 3,631 | Independent survey with correlated subbands. Store fit error, error on the common scale, internal-scale error and local rms. The wideband entry overlaps the subbands. |
| LoTSS DR3 | 144 | 25″ | 18,840 | Integrated and peak fluxes, with counts and sums inside 25, 45, 75, 135 and 180″ for blend checks. |
| TGSS ADR1 | 150 | 25″ | 18,000 | Integrated with native error. Retain mosaic fit size, peak and structure code. |
| WENSS | 325 main, 352 polar | 45″ | 16,079 | Integrated. The error column is local rms. Inspect fit-problem and source-type flags, especially groups and their separate components. |
| NVSS | 1,400 | 30″ | 19,209 | Integrated with catalogue uncertainty. Retain fitted size and residual flags. |

The common flux columns use Jy, except VLSSr peak brightness (Jy/beam). The native shape and auxiliary flux columns retain HEASARC units. The code uses great-circle distances. Unmatched rows carry empty flux fields. Catalogue radii and counts alone cannot distinguish detector limits from incomplete coverage. GLEAM's extragalactic catalogue ends around +30° and has sky masks. Its absolute scale error is 13% at the northern edge. [GLEAM](https://doi.org/10.1093/mnras/stw2337) · [LoTSS DR3](https://lofar-surveys.org/dr3.html) · [VLSSr](https://heasarc.gsfc.nasa.gov/w3browse/all/vlssr.html) · [8C survey](https://www.mrao.cam.ac.uk/projects/surveys/8C/)

## Reverse search and repeat observations

VLSSr supplied a parent independently of LoDeSS detection. Select sources brighter than 1 Jy/beam with both convolved axes below 90″ within 4° of each of the 13 good-quality field centres. A match required a same-field LoDeSS centroid within 60″. This yielded 632 field/source pairs (622 unique VLSSr identifiers), of which 155 pairs (153 objects) lacked a LoDeSS catalogue centroid in the corresponding field. Forced point-source fluxes and empirical local errors are available for all 155. There were 35 measurements below three times their local uncertainty. Requiring an independent catalogue match with flux above 20 mJy (NVSS) or 80 mJy (LoTSS or TGSS) retained 18 different objects. All 18 were in P136+36. The low-flux appearances come from this one pointing and are spatially and instrumentally correlated.

Every one of the 18 was detected at its original VLSSr position in P136+30, with S/N from 8.92 to 53.40. Median S/N was 1.04 in P136+36, 22.43 in P136+30 and 16.37 in P128+36. The three pointings were observed on 2022-04-14, 2021-12-02 and 2022-11-24, respectively. For VLSSr J090432.3+353906 the forced brightnesses in the same order were 0.035 ± 0.043, 0.952 ± 0.039 and 1.250 ± 0.026 Jy, against a VLSSr peak brightness of 3.744 Jy/beam. The errors here describe local noise only. Field-dependent losses undermine an extreme spectral cutoff inferred from P136+36. The other 17 low-S/N cases had no counterpart meeting the cuts. For example, several near 3C 239 have no corresponding bright emission in NVSS, TGSS or LoTSS. Their archival VLSSr entries alone are inadequate to establish an astrophysical spectrum.

The signed image flux is stored and should be used directly in a flux likelihood. The one-sided quantity called `limit95_jy` is a screening statistic, max(0, measured flux) + 1.645 × empirical noise. It describes apparent flux in a single image under Gaussian thermal scatter. Field-dependent attenuation and PSF errors prevent using it as a bound on intrinsic emission. The evidence tables retain per-field parent entries, all candidate image measurements, external matches and repeat observations.

## Independent coverage of ILTJ112022.08+502936.2

The NRAO VLSSr postage server returned a 1° FITS cutout from field 1100+550, observed 2003-09-20 at 73.7939 MHz. Forced photometry at the LoLSS position gave 0.408 ± 0.084 Jy before CLEAN-bias correction using the published 75″ restoring beam. The aperture estimate is 0.356 Jy and the brightness peaks 7.14″ from the target. Forty nearby fixed-aperture controls provide the empirical statistical error. The survey's point-source CLEAN-bias prescription adds about 0.047 Jy here. This earlier VLA observation independently anchors the spectrum around 74 MHz. [VLSSr image service](https://www.cv.nrao.edu/vlss/VLSSpostage.shtml) · [Lane et al. 2014](https://doi.org/10.1093/mnras/stu256)

The VCSS epoch-one bright catalogue gives a 340 MHz flux for VCSS1 J112022.1+502935, offset 0.36″. The raw integrated value is 220.05 ± 4.46 mJy, and the corrected integrated value is 251.72 mJy. The quoted fit error is on the uncorrected measurement. The fitted component is single, 14.5″ × 13.6″, with S/N 85.46. This constrains the higher-frequency spectrum. [VCSS release](https://cdsarc.cds.unistra.fr/viz-bin/cat/VIII/111)

The actual NASA LAMBDA 36.528 MHz OVRO-LWA map has target sky brightness contrast of −925 K relative to a local annulus. The normalized median absolute deviation of 48 equivalent nearby measurements is 1,372 K (about 1.84 Jy/beam). The published thermal noise and estimated confusion are 0.799 and about 1.64 Jy/beam. Other LoTSS detections contribute approximately 0.320 Jy at 144 MHz in a Gaussian beam with 18.48′ FWHM, compared with the target's 0.380 Jy. This is a 144 MHz confusion diagnostic. The target's spectrum at 36.528 MHz cannot be separated in this map. [Map and beam](https://lambda.gsfc.nasa.gov/product/foreground/fg_ovrolwa_radio_maps_info.html) · [Eastwood et al.](https://arxiv.org/abs/1711.00466)

LoLSS DR1 re-ran the calibration of the partially calibrated preliminary-release data. The two catalogues share observations, so their independent processing cannot serve as a second observing epoch. The VLA 4MASS programme became the 74 MHz VLSS survey. [LoLSS DR1 §3](https://doi.org/10.1051/0004-6361/202245389) · [NRAO VLSS](https://www.cv.nrao.edu/4mass/4MASSlist.shtml)

## NVSS J100005+812702

The reverse search recovered this source in the P139+79 image even though the deduplicated source catalogue has no P139+79 centroid within 60″. The full DR1 catalogue contains it as zero-based row 9460 in P210+80, 0.83″ from the LoTSS position. The Gaussian catalogue includes components from all three overlapping pointings. LoDeSS measures 0.671 ± 0.025 Jy by forced PSF in good P139+79, 3.02° from field centre. P210+80 and P264+86 give 1.060 ± 0.054 Jy and 0.435 ± 0.063 Jy at 9.22° and 10.48° in mediocre fields. All three errors here describe local noise. Image apertures of radius 60″ give 0.701, 0.808 and 0.135 Jy respectively. The third field has a large residual and a strongly distorted PSF. The SED uses the best-pointed epoch with a 39% gain error and substitutes the brightest repeat in a separate fit. [Pointing images and complete component associations](../release/evidence/J100005/)

The NRAO VLSSr catalogue gives a 2.191 Jy/beam peak with measured axes 88.44″ × 75.00″. In the restored 75″ image, a point-source-normalized 90″-radius aperture at the LoTSS position gives 2.273 Jy. Forty control apertures give 0.203 Jy scatter. The local point-source CLEAN correction of +0.080 Jy gives 2.353 Jy, with an additional 15% calibration allowance. Corrected fluxes in 75″ to 105″ radii span 2.322–2.365 Jy. The more conservative catalogue Gaussian-area estimate, corrected the same way, gives 2.663 Jy. A separate fit uses this higher alternative. The primary VLSSr model data point is the image integral. [VLSSr calibration](https://doi.org/10.1093/mnras/stu256) · [Apertures](evidence/j100005_vlssr_apertures.csv)

At 144 MHz the LoTSS DR3 source has morphology code M and deconvolved major axis 4.54″. Its catalogue total is 3.159 Jy and the sum of other catalogue sources within 180″ is 0.00088 Jy. Images at 6″ and 20″ give 3.125 and 3.031 Jy in 60″ apertures. At common 75″ resolution these same apertures give 3.115 and 3.017 Jy. TGSS measures 3.489 ± 0.349 Jy. The 352 MHz WENSS polar entry gives 2.541 Jy and NVSS measures 0.863 Jy. The 1988 Miyun field independently measured an unresolved 2.13 Jy/beam at 232 MHz. Its S/N of 41.8 implies 0.051 Jy/beam local noise. The Miyun sensitivity fit allows 16.6% calibration scatter and a further fit uses 20%. [LoTSS aperture test](evidence/j100005_lotss_apertures.csv) · [WENSS polar frequency](https://support.astron.nl/LOFARImagingCookbook/bbs.html) · [Miyun survey](https://doi.org/10.1051/aas:1997314)

| Smooth spectral fit | ΔAIC, power law minus curved fit | Peak (MHz) | 95% profile interval (MHz) |
|---|---:|---:|---:|
| Six surveys, including LoDeSS | 30.80 | 187 | 128–221 |
| Five surveys, 23 MHz held out | 20.71 | 162 | 80–219 |
| Four surveys, all LOFAR held out | 17.94 | 164 | 76–223 |
| Add Miyun | 28.15 | 181 | 118–218 |
| Add Miyun, hold out 23 MHz | 18.61 | 152 | 67–212 |
| Replace 23 MHz with brightest repeat | 30.40 | 169 | 115–202 |
| Add 20% independent scatter per epoch | 12.99 | 198 | 129–249 |

The intervals profile the log-parabolic peak after refitting its amplitude and curvature, using Δ(−2 log L) = 3.84. They assume a smooth spectrum and the assigned calibration model. The comparison is source-specific and has no search-wide false-alarm calibration. The 6C and 7C catalogues report 2.54–3.55 Jy near 151 MHz. This historical spread, the scatter between LoDeSS pointings and the strength of the known compact centimetric component limit claims about intrinsic variability and the precise location or cause of the turnover. [Fit results](evidence/j100005_fits.csv) · [Input fluxes](evidence/j100005_photometry.csv) · [Historical SPECFIND entries](evidence/VIII_104_1000.tsv)

No 38 MHz 8C source falls within 17.5′. The closest source is the 27.5 Jy 8C 1001+817. The Cambridge archive describes the published images as microfiche and distributes a corrected source list. Branson's independent 81.5 MHz polar catalogue has a 10.3 Jy source 15.9′ away, NB81.14. Its nominal depth excludes areas near bright sources. Neither catalogue supplies a reliable target flux limit here. GLEAM does not cover this polar position. [Revised Rees 38 MHz survey](https://www.mrao.cam.ac.uk/projects/surveys/8C/) · [Branson, Table III and completeness discussion](https://articles.adsabs.harvard.edu/pdf/1967MNRAS.135..149B)

## Reproduction

Run from this directory with Python 3.12, NumPy, SciPy, pandas, Astropy and curl. First follow the adjacent release instructions to obtain the LoDeSS catalogues, field mapping and the images used below. Set DATA to that directory. The HEASARC inputs are pinned by SHA-256 in archives.csv. The ASTRON retrieval batches 500 uploaded positions in one query and returns all catalogue components in a 3′ cone per source. The script validates its CSV columns and the uploaded identifiers. Delete the corresponding cache if the source FITS row order changes.

```bash
DATA=../release/data
export PYTHONPATH="../release:..${PYTHONPATH:+:$PYTHONPATH}"
python fetch_archives.py raw
python fetch_lotss.py "$DATA/sources.fits" raw
python match_catalogues.py "$DATA/sources.fits" raw results
python reverse_vlssr.py "$DATA" .
python fetch_lotss.py results/reverse_positions.fits reverse_raw
python match_catalogues.py results/reverse_positions.fits reverse_raw reverse_results --surveys lotss
python match_catalogues.py results/reverse_positions.fits raw reverse_results --surveys gmrtas150m nvss rrs8c38mhz wenss gleamegcat
python inspect_reverse.py .
python reverse_controls.py "$DATA" .
python measure_vlssr.py raw/vlss_1120.fits results/vlss_1120.json
```

The compressed files for the complete parent and crossmatches can be regenerated with these commands. The local-control field offsets are ≤10° and the fixed empirical PSF aperture must be fully covered for a repeat to be accepted. To reproduce the OVRO-LWA contrast, use the Mollweide map listed below and ASTRON's 144 MHz source catalogue.

```bash
curl -fLsS https://lambda.gsfc.nasa.gov/data/foregrounds/ovro/images/ovro_mollweide_lwa_sky_map_36.528MHz.fits -o raw/ovro_36.fits
curl -fsSG https://vo.astron.nl/tap/sync --data-urlencode REQUEST=doQuery --data-urlencode LANG=ADQL --data-urlencode FORMAT=csv --data-urlencode "QUERY=SELECT source_name,ra,dec,total_flux FROM lotss_dr3.main_sources WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',170.092166,50.493238,0.6))" -o results/ovro_neighbours_1120.csv
python check_ovro.py raw/ovro_36.fits results/ovro_neighbours_1120.csv results/ovro_1120.json
```

The source-specific script uses the archived reverse-counterpart table under `evidence/`. It retrieves the external images on first run, records their checksums and creates the aperture checks, input photometry, model fits and plot. The two main LoDeSS fields are required. The separate release analysis contains the third epoch from P264+86.

```bash
python ../release/acquire_lodess.py "$DATA" --images P139+79 P210+80
python j100005.py "$DATA" raw evidence
```
