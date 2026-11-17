## DECam and Legacy Surveys

**No two-night moving point source passes the search at the trial Phan velocities.** Two public 90-second DES DECam r exposures, 605285 (2016-12-30) and 605571 (2016-12-31), observed the same CCD 21.47 hours apart. The two-lobe significance peaked at 5.23 in a dense sky-position and displacement search with threshold 5.5. At the 3,617 discrete model tracks, 3,424 passed the masks and none passed the detection criterion. The models, not a binned sky density or physical posterior, define these denominators.

### What the optical model tests exclude

Each cell gives the number of artificial planets recovered out of **400** trial draws, in geometric albedo order **0.10 / 0.30 / 0.75**. Every draw includes masked and outside-search outcomes in the denominator. The full 42-case table includes 95% Wilson intervals, realized r magnitudes, eligibility counts, and additional envelope models. [model_recovery.json](model_recovery.json)

| Assumed mass and radius | Published centroid draws | Centroid stress draws |
|---|---|---|
| 5 M⊕, 2.71 R⊕ | 379 / 389 / 389 | 371 / 379 / 379 |
| 5 M⊕, 2.94 R⊕ | 385 / 389 / 389 | 379 / 379 / 379 |
| 10 M⊕, 3.70 R⊕ | 389 / 389 / 389 | 379 / 379 / 379 |

The trials draw uniformly with replacement from the bound nodes of a separate 13 × 32 × 56 orbit grid. The stress experiment assigns equal probability to the 25 pairs formed from the central IRAS direction and the four endpoints of its reported 18″ × 7″ error ellipse (position angle 79°), and the AKARI central direction plus four cardinal 40″ offsets. These weights and the AKARI offset size are sensitivity choices. Neither a probability law for the actual candidate nor a confidence region for AKARI astrometry is available. The Wilson intervals quantify trial-draw precision for each stipulated design. They omit uncertainties in radii, albedos, spectra, or the astrometric model.

For each draw, the calculation uses the object's actual propagated solar distance and observer distance at the DES epoch, with light time. The reflected flux uses full phase, solar Sloan r = −26.93 AB, and the stated radius and geometric albedo. It assumes the adopted solar r calibration is applicable in the DECam r band; an object's spectrum can change that transformation. Each model is planted in both registered science images on its own trajectory with the local empirical PSF, Poisson shot noise and a 3.83 e⁻/ADU gain. The second visit uses the relative photometric scale from field stars. The code recomputes data-dependent masks on the planted patches and runs the complete two-lobe detection and localization rule. Saturation is not reached by these modeled planets. [planet_injections.py](planet_injections.py)

### Motion search and validation

The search uses FITS midpoint times, TPV astrometry, pipeline DQ masks, local sky estimates and a residual astrometric translation from the field stars. Each epoch is convolved with the other exposure's empirical PSF. The difference is filtered with the resulting common PSF. The search compares a positive residual at the earlier predicted position with a negative residual at the later one. Both lobes must be present, each must exceed 1.75 times the dipole noise divided by √2, and the lobe flux ratio must be between 0.25 and 4. Randomly displaced control dipoles supply the local robust noise scale. A threshold of 5.5 is applied to the summed signed response. This is a search threshold in units of empirical noise, not a Gaussian global false-alarm significance. Fluxes are calibrated using the Legacy DR11 r template only as a field-star reference; the moving-source difference uses independent science frames. The 31 local stars give robust single-stamp zeropoint scatters of 0.075 and 0.037 mag. [difference_result.json](difference_result.json)

The rectangular starting-position search extends ~76″ beyond the nominal model bounds, clipped to the CCD interior, and evaluates 10 distinct rounded two-epoch motion vectors. At least one motion passed the geometric and detector masks at 1,342,312 of 1,431,771 starting pixels. Five first-pass features surviving the two-lobe flux test were reviewed in the original pixels. Four lie on narrow curved or branched tracks in the first exposure, and one lies next to a compact sharp spike. Their later predicted positions lack a matching point source. A 3 × 3 median sharp-feature mask with a large fractional-gradient requirement removes these particle events; PSF injections include a bright-mover check on the resulting masks. No candidate reached minor-planet identity or variable-star vetting. The audit montage shows the first-exposure position in cyan and the second in red. [First-pass feature coordinates](pre_morphology_peaks.json) · [Original pixel stamps](dense_review.jpg)

An independent sensitivity check used 400 off-grid trajectories, each drawn from a finer grid and shifted within a uniform-area 12″ disk. The blind bank search recovered 387 at r=23, 373 at r=23.25, 235 at r=23.5 and 18 at r=24. Injections used the real background, empirical PSF, source shot noise and a 2.5-pixel localization criterion. The fixed-magnitude runs retain the measured field masks, while the physical-model runs also reapply those masks at each planet brightness. One PSF per epoch, catalogue centroid treatment, and the absence of a calibrated posterior limit the interpretation. [dense_injections.json](dense_injections.json)

### 2024–2026 imaging

| Public imaging tested | Coverage and result | Injections |
|---|---|---|
| 2024-11-12, DECam g, 2 × 10 s | Compared two separate frames to Legacy DR11 g on CCD 50. Of 3,617 nodes, 2,200 passed both masks, and none passed joint S/N ≥5.5 and individual S/N ≥2. Peak joint S/N 3.57. The 62-second baseline does not resolve the expected motion. | Of 320 random nodes, 198 were eligible. Recovered 198 at g=21, 170 at g=21.5, 14 at g=22. Photometric scale from four local stars in each exposure. |
| 2025-07-22, DECam N395, 720 s, proposal 2023B-646244 | Tested two CCDs of the 395 nm MAGIC exposure with empirical matched filters, DQ masks and a pre-existing Legacy g counterpart mask. Of 3,617 trial nodes, 1,466 passed all cuts. None passed 5.5. Peak S/N 3.50. | Of 320 draws, 129 were eligible. Recovered 101 at a total flux of 80 ADU, and all 129 at 150 ADU. The frame names a Gaia g calibration; the calibrated narrowband AB scale and empirical-PSF injections are now available in the [independent Gaia XP analysis](../n395_cal/README.md). |
| 2026-09-13, 10-28, 11-15, DECam r, 12 × 900 s | The archive lists four dedicated exposures per night, with object names `p9-candidate` and `p9-1hour-20261114`. Their release dates are 2028-03-05, 2028-04-19 and 2028-05-07. An anonymous retrieval on a November frame returned HTTP 403 and the release date. | No public pixel limit. |

DR11 g observations end 2024-09-11, before the 2024 November exposures. Their independent template removes static sky, with a deep-source mask around significant stationary objects. The pair's shorter exposure and calibration sample limit its depth. The N395 frame supplies one epoch. A subsequent narrowband calibration uses the published MAGIC system-response curve, nine Gaia XP standards and realistic stellar PSFs. It measures 126/129 eligible injections recovered at N395_AB 23.25 and 91/129 at 23.50. [Standards, injection results and spectral assumptions](../n395_cal/README.md) [Legacy DR11 observing dates](https://www.legacysurvey.org/dr11/description/) · [MAGIC survey](https://arxiv.org/abs/2605.26581) · [2024 recovery](modern_g_result.json) · [2025 recovery](modern_n395_result.json) · [2026 receipts](dedicated_2026.csv)

### Reproduction and outstanding checks

Run with Python 3.12, NumPy, SciPy, Astropy, PyERFA, and Matplotlib. The repository's sibling `orbits/` directory holds the propagation code. `download.py` verifies SHA-256 before using FITS inputs; the two dynamic archive queries for exposure inventory are snapshotted as gzip files. `inventory.py` propagates the Phan grid to 226 distinct exposure midpoints returned from a broad NOIRLab rectangular boresight search. The 389 file rows include repeated reductions. Of those 226 midpoints, 42 have positive CCD intersections demonstrated by independently queried legacy Data Lab polygons; blank cells denote missing matching polygon data. The newer 2024/25 WCS intersections are measured separately. A boresight match alone does not establish coverage. Other unprocessed historical visits and broad-band configurations may yield better constraints. [Exposure inventory](decam_epochs.csv) · [Archive queries](archive_queries.json) · [Image hashes](decam_inputs.json) · [2024/25 focal plane intersections](modern_footprints.json)

```bash
python download.py --group des
python decam_difference.py
python download.py --group all
python modern_g.py
python n395.py
python inventory.py
```

This work uses the NOIRLab Astro Data Archive and Astro Data Lab, DES imaging and Legacy Surveys coadds. The optical survey products can reuse underlying DECam exposures. Their sensitivities should not be multiplied as independent probabilities. [NOIRLab archive](https://astroarchive.noirlab.edu/) · [DES data access](https://www.darkenergysurvey.org/the-des-project/data-access/) · [Legacy Surveys](https://www.legacysurvey.org/)
