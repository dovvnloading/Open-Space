# DECam N395 calibration and moving-source limits

The July 22, 2025 N395 frame contains no source above the original threshold of 5.5 on the 1,466 usable nominal Phan tracks. Artificial planets with the shapes of observed Gaia stars recovered on 126 of 129 eligible tracks (97.7%) at **N395 AB 23.25** and 91/129 (70.5%) at AB 23.50. Of 3,617 nominal trajectories, 40.53% pass the geometric and quality cuts. Combining the exact measured coverage with the conditional injections gives 39.6% nominal-grid sensitivity at AB 23.25. A single epoch supplies a brightness constraint and cannot establish motion.

| N395 AB | Among 129 eligible | CCD 15 (116 trials) | CCD 21 (13 trials) |
|---|---:|---:|---:|
| 23.00 | 129 | 116 | 13 |
| 23.25 | 126 | 114 | 12 |
| 23.50 | 91 | 87 | 4 |
| 23.75 | 34 | 30 | 4 |
| 24.00 | 8 | 7 | 1 |
| 24.50 | 0 | 0 | 0 |

## Throughput and zeropoint

We extracted the MAGIC team's N395 instrumental response, including telescope and mirror response, from the dashed vector path in the authors' original Figure 2 artwork. Its transmission-weighted wavelength is 395.39 nm and its FWHM is 9.90 nm. For each Gaia DR3 sampled XP spectrum we computed the photon-counted passband integral, dividing by the equivalent response to a flat 3,631 Jy source. The TAP cone rectangle contains 63 stars with continuous spectra, of which 21 have the bright-source sampled spectra used here. Seventeen fall inside the two CCDs with adequate border margins. Nine pass the variability, RUWE, spectrum S/N, image S/N, saturation, DQ and centroid checks. The astrometry includes Gaia proper motions to the frame date. AIP's TAP response contains signed 343-point spectra and their sampled errors. The measurement table retains the flux, local image uncertainty and reason for excluding each rejected star. [MAGIC Figure 2](https://arxiv.org/abs/2605.26581) · [Gaia DR3 sampled spectra](https://doi.org/10.17876/gaia/dr.3/52) · [Gaia external BP/RP calibration](https://doi.org/10.1051/0004-6361/202243880)

| CCD | Accepted standards | 13-pixel-aperture AB zeropoint, per ADU in the 720 s frame | Standard-to-standard RMS (mag) | Approximate random error on median (mag) |
|---|---:|---:|---:|---:|
| 15 | 4 | 28.627 | 0.055 | 0.060 |
| 21 | 5 | 28.626 | 0.023 | 0.041 |

The calibration uses the 13-pixel radius consistently for stars and injected flux. Random errors include the larger of the measured star-to-star MAD and the median per-star spectral/photometric uncertainty, scaled for the sample median. Globally correlated errors in the physical XP flux scale, atmospheric color terms at the observed airmass of 1.09 and XP resolution through the Ca H and K absorption lines remain. No color transformation was fit with this small standard set. Replacing the published response with a 10.01 nm rectangular passband changes the median synthetic magnitude by 0.020 and 0.019 mag on the two CCDs. That difference tests bandpass choice only.

The original scoring kernels enclosed 74.5% and 93.4% of their volume within two pixels on CCDs 15 and 21. Recognized Gaia stars enclose approximately 52–54% at the same radius. Direct linear conversion of the original artificial ADU threshold using a single detector-independent scale would misrepresent a stellar image. We retained the original sky subtraction, mask, empirical detection kernel and 5.5 threshold, and injected sky-subtracted profiles of the retained same-CCD Gaia stars into the original science pixels. These include the true predicted 0.15–0.21-pixel motion during the exposure, fractional pixel phase, Poisson shot noise at the CCD gain and a fresh sharp-feature veto following injection. Existing fixed-source and bad-pixel masks are retained. The replay reproduces the real-data maximum S/N of 3.499 and both original CCD eligibility counts.

## Planet brightness and recovery

The absolute CALSPEC solar spectrum integrated through the measured passband gives apparent N395 AB −25.598 at 1 au. We use a fully illuminated, grey reflecting disk within N395 and apply each orbital state's heliocentric and geocentric distances separately. The albedos are geometric and refer specifically to 395 nm. A planet's methane, haze, clouds or other absorption can make them differ from its broad visible albedo. The 37–48 K thermal models contribute negligible thermal continuum here. [CALSPEC solar spectrum](https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/calspec)

| Physical model | N395 AB on the nominal grid | Recoveries among 129 eligible |
|---|---:|---:|
| 5 Earth masses, 2.71 Earth radii, albedo 0.10 | 23.53–25.04 | 23 |
| Same mass and radius, albedo 0.30 | 22.34–23.85 | 120 |
| Same mass and radius, albedo 0.75 | 21.34–22.85 | 129 |
| 10 Earth masses, 5.09 Earth radii, albedo 0.10 | 22.16–23.67 | 128 |
| Same mass and radius, albedo 0.75 | 19.97–21.48 | 129 |

For the compact, low-albedo five-Earth-mass model the empirical coverage-scaled whole-grid sensitivity is 7.2%. For the extended ten-Earth-mass model at the same albedo it is 40.2%. These are finite-grid results for this date, model and field mask. They include neither a continuous-orbit prior nor the complementary coverage of other epochs. All 30 radius and albedo combinations are in `model_recovery.csv`.

## Exact DES–Euclid misses

We reused the saved 6D state vectors of the dark 5 Earth-mass, 2.71 Earth-radius, albedo-0.10 planets that both DES and Euclid missed in their matched experiment. The three published-centroid draws predicted N395 AB 24.87–25.04. Only one survived the original N395 mask, with real S/N +0.57, and its artificial planet missed on all 64 repeated Poisson and empirical-PSF trials. The twelve centroid-perturbed misses predicted AB 23.86–24.91. Six were eligible, four lay in the Legacy stationary-object mask and two fell outside the tested CCDs or their edge margins. Their eligible raw pixel S/N ranged from −1.80 to +1.02. None of the six first seeded injections passed. For stress draw 389 (AB 23.86) three of 64 repeats passed, and none passed for the other five. This visit added zero recoveries on the predefined first realization of those 15 specific missed states. Their limited sensitivity and masks account for the remaining uncertainty.

The per-draw N395 position, magnitude, selection, real-data S/N and artificial-source response are in `paired_misses.csv`. Append `--paired-root /path/to/research/result` to the injection command to reproduce them. The paired root contains `archives/planet_draws.csv.gz` and `orbits/euclid/paired_trials.csv.gz`. The original two-survey experiments and this narrow-band check are linked by draw number and saved orbit rather than independent per-survey efficiencies.

## Reproduction

The catalogue and XP response are bundled as compressed raw responses. The input SHA-256 in `receipts.json` applies to their decompressed bytes. That file also identifies the authors' source archive, the optical input FITS and mask, the two Legacy reference images and the CALSPEC solar spectrum. To recreate `bandpass.csv`, pass the downloaded arXiv source archive to `python bandpass.py PATH_TO_ARCHIVE` (PyMuPDF and NumPy). The published PDF vector is sampled directly; the response outside its available wavelength interval is taken to be zero.

Use the repository's optical download script to populate the four native image paths, then run the commands below with NumPy, SciPy, Astropy and PyERFA. The solar FITS is fetched if absent. The two optical input scripts import the same detection routines and orbit propagator as the parent analysis. `standards.csv` reports per-star physical flux in Jy and its conservative uncertainty, aperture fluxes, matched amplitudes, quality cuts and residual astrometry. `injections.csv.gz` and `planet_trials.csv.gz` retain individual injected track IDs, calibrator IDs and observed outcomes.

```sh
python calibrate.py --optical-dir /path/to/optical --orbit-dir /path/to/orbits --cache /path/to/optical/receipts
python planet_injections.py --optical-dir /path/to/optical --orbit-dir /path/to/orbits --cache /path/to/optical/receipts
```
