# SPHEREx moving-source search

None of the 2,551 eligible nominal trajectories (70.53% of 3,617) passed the three-season test in any of the six SPHEREx detectors. We measured 572 DR1 and QR3 images with observation epochs from June 2025 through July 2026. The union of 901 eligible translated controls also had zero triggers. A flat-spectrum injected source was recovered at AB 18.0 on 234/234 eligible tracks, at 19.0 on 194/234, at 19.5 on 158/234 and at 20.0 on 13/234. At 20.5 none passed. For the brightest tested reflected-light model (10 Earth masses, 5.09 Earth radii, geometric albedo 0.75), detections on 151 of 320 injected trajectories imply estimated nominal-grid sensitivity of 45.5% after applying the exact measured coverage. No albedo-0.1 model was recovered.

| Detector | Wavelength range (µm, approximate) | Native exposures | Eligible grid states | Faintest tested AB magnitude with ≥90% recovery among eligible injections |
|---|---:|---:|---:|---:|
| D1 | 0.75–1.12 | 114 | 1,249 | 19.0 |
| D2 | 1.10–1.65 | 82 | 730 | 19.0 |
| D3 | 1.62–2.43 | 93 | 570 | 19.5 |
| D4 | 2.42–3.82 | 107 | 272 | 19.5 |
| D5 | 3.81–4.42 | 89 | 591 | 18.5 |
| D6 | 4.41–5.00 | 87 | 1,093 | 18.0 |

The query returned 671 exposures across six detectors (421 DR1, 250 QR3). Time, sky coverage, flag requirements and available reference images differ by band. Ninety-nine listed observations did not supply a usable cutout in these requests. Each band's `images.json` records its acquisition status. Positions and controls use common IDs. `combine.py` tests whether any same-detector three-season search passed, counting a shared trajectory only once. `combined_injections.csv` gives measured rates at 14 constant magnitudes and for all 30 planet radius/albedo cases. The near-IR physical test approximates the planet as a grey reflector. Cold planets with strong molecular absorption or thermal atmospheric windows can have different spectra.

### D1 measurement detail

D1 had 1,249 eligible trajectories (34.53%). Native image injections recovered 124 of 124 eligible flat-spectrum objects at AB 19.0. Recovery declined to 105 of 124 at AB 19.5 and 4 of 124 at AB 20.0. The denominator refers to the fixed grid fitted to the two published positions. It carries no posterior probability and excludes astrometric perturbations and alternative coadd weighting.

## Completed image measurements

The IRSA Simple Image Access cones had angular radius 0.40° about ICRS (35.04°, −49.67°), a 0.90 µm detector selection, and a 5,000-row limit. DR1 returned 71 D1 observations and QR3 returned 50. The native Level 2 image cutouts were 0.78° squares at this fixed center. 114 hashes represent usable downloads (46 June 2025, 23 December 2025 and 45 June–July 2026). The seven additional reported footprints gave one HTTP 422 and six HTTP 503 responses on the cutout endpoint after retries. Each used file contains calibrated surface brightness, flags, variance, zodiacal brightness and a spatially varying ePSF. D1 covers approximately 0.75–1.12 µm in narrow slices with a linear-variable filter. [SPHEREx public data](https://irsa.ipac.caltech.edu/Missions/spherex.html) · [SPHEREx data products](https://irsa.ipac.caltech.edu/data/SPHEREx/docs/)

Each trajectory was propagated to the midpoint of every exposure. The WCS with distortion gives detector coordinates. The statistic uses the native IMAGE minus ZODI brightness, converted from MJy/sr into µJy per pixel using the projection-plane pixel area. Pixels require FLAGS = 0, finite data and positive variance. A locally normalized difference of 0.8- and 3-pixel Gaussian smoothings subtracts background. Each smoothing requires at least 70% of its kernel support to be valid. At each trial position, the median in reference images from the other seasons subtracts persistent structure. The reference set uses up to four exposures from each other season, prioritizing coverage and matching the archived ePSF central wavelength at the median track position. Each science measurement needs at least two usable reference images. Each season requires at least three usable science measurements. The final statistic uses the median in each season.

We translated 1,500 randomly drawn orbit tracks by a fixed sky offset per track between 8′ and 13′. Of these 437 met the same coverage and masking rules. The median and 1.4826 × MAD across these controls calibrate the statistic separately by season. The detection requirement was a statistic above +2.5 control scales in **every** season and a mean above +4. None of the 437 eligible controls passed. For the target grid, the largest minimum statistic was +0.993 with a mean of +1.538. Both thresholds were fixed before the search. Shared images and nearby trajectories make these trial statistics correlated.

| Constant AB magnitude in D1 | Flux density (µJy) | Recovered among eligible injections | Recovery among all 320 injected grid points |
|---|---:|---:|---:|
| 19.0 | 91.21 | 124/124 (100.0%) | 124/320 (38.8%) |
| 19.5 | 57.55 | 105/124 (84.7%) | 105/320 (32.8%) |
| 20.0 | 36.31 | 4/124 (3.2%) | 4/320 (1.3%) |
| 20.5 | 22.91 | 0/124 | 0/320 |

The 320 injections are uniformly drawn, without replacement, from the nominal grid (seed 20261117). The eligible fraction among them fluctuates around the exact 1,249/3,617 coverage. The estimated whole-grid sensitivity at AB 19.0 is 34.5% when conditional recovery is multiplied by known coverage. The binomial 95% interval for the measured 124/124 conditional recovery is 97.1–100%; at AB 19.5 it is 77.1–90.5%. These describe the numerical sampling of tracks given these images. They do not account for ePSF and spectral-model systematics.

At each track and epoch, injections use the nearest position-dependent archived ePSF, sampled at the source's detector pixel phase and normalized to a total flux. We add the source to native science stamps and refilter them with their original mask. Linearity gives the response at each flux before the nonlinear image and season medians. This propagates self-subtraction when the source lies at a reference location. Masks and empirical noise are held fixed and extra source shot noise is neglected. Astrometric phase, measured backgrounds and flagged pixels vary across trials. The replay uses the subset of original trials whose true positions were injected. It measures recovery at the true track. A blind search could also trigger on a neighboring trial.

## Planet models

`grey_model_injections.csv` repeats the native-pixel test for the 5, 8 and 10 Earth-mass model radii in the accompanying brightness grid. It uses a constant geometric albedo of 0.1, 0.3 or 0.75 throughout D1, full phase, and a solar blackbody of 5,772 K and radius 695,700 km. The local ePSF wavelength, heliocentric object distance and geocentric range set each trial exposure's flux. For an object radius R (in au), the flux is the monochromatic solar flux at 1 au multiplied by albedo × R²/(heliocentric distance² × geocentric range²), all distances in au. The 37–48 K interior thermal model contributes negligible thermal continuum in D1. Deep absorption and clouds can change the near-infrared albedo. Objects with the same mass span different radii and intrinsic temperatures. The model test reports the outcome of each injection as its own flux changes with wavelength and orbit. [Fortney et al.](https://arxiv.org/abs/1604.07424) · [Linder and Mordasini](https://arxiv.org/abs/1602.07465)

Across all 320 injected tracks the compact 5 Earth-mass, 2.71 Earth-radius model recovered 0, 0 and 2 tracks for albedos 0.1, 0.3 and 0.75. An extended 5 Earth-mass, 4.12 Earth-radius planet recovered 0, 0 and 59. An extended 10 Earth-mass, 5.09 Earth-radius planet recovered 0, 14 and 112. The latter gives a raw 112/320 (35.0%) sample recovery overall and 112/124 (90.3%) among eligible injections. Applying that conditional recovery to the exact grid coverage gives 31.2% estimated sensitivity for that one bright model. All ten models at albedo 0.1 gave zero recoveries. The full model range of trial-median D1 brightness was AB 18.52 to 23.51.

## Reproduction and remaining coverage

Run `python search.py --orbit-dir /path/to/orbits --cache /path/to/cache` with NumPy, SciPy, Astropy and PyERFA installed. The orbit directory must contain `propagate.py`, `orbit_input.json` and `brightness.csv`. The script fetches successful entries from `images.json` and checks each SHA-256 before measurement. The two compressed CSV files contain the raw SIA responses. `query_receipts.json` identifies the exact requests, acquisition time and SHA-256 of the decompressed response. `eligible_grid_indices.json` supplies the zero-based indices into the nominal orbit grid. The SPHEREx spacecraft's displacement from the geocenter is much smaller than one survey pixel at 500–700 au. The narrowband ePSF central wavelength approximates the local spectral response. A detailed planet atmosphere and per-channel passband integral would be needed for spectrally structured predictions.

For D2–D6, add `--manifest D2/images.json --out D2`, replacing D2 with the detector. Run `python combine.py` once all detector products are present. For each detector `injection_trials.csv.gz` records the zero-based orbit ID and every pass/fail result. `all_query_receipts.json` contains the all-wavelength search requests and the SHA-256 of the adjacent compressed CSV files' decompressed responses.

The measured nominal-grid fractions refer to exactly those trial trajectories. In a separate spacing check, 15,841 bound continuously drawn nominal orbits had median 5.72″ maximum seasonal displacement from the closest nominal grid track and a 95th percentile of 11.95″. Only 55.2% lay within one 6.15″ native pixel at all nine representative epochs. With a sampled 18″ × 7″ IRAS centroid ellipse and a 40″ AKARI disk, the corresponding fraction dropped to 16.1%. `orbit_spacing.py` and `orbit_spacing.json` document the offsets. Model colors, grid spacing, flagged pixels and reference coverage all affect the completeness of a blind continuous-orbit search.
