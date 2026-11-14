# DECam recovers the OGLE oscillation on the image of 2dFS 757

The DELVE DR3 Magellanic Clouds catalog places source 0129m727.52197 within 0.01 arcseconds of the position of 2dFS 757. In **156 unforced DECam exposures on 26 nights**, its resolved photometry recovers **1.107 ± 0.102 times** the amplitude of the 4.88575-hour OGLE-IV waveform at its predicted phase. A scale of one means the same semiamplitude as OGLE's 9.83 millimagnitudes. Omitting each night in turn gives scales from 1.030 to 1.158. The analysis used this telescope only after freezing the period and phase with OGLE. [Fitted scales, uncertainties and all controls](result.json)

| Band | Distinct exposures | Nights | Scale ± conditional night-cluster SE |
|---|---:|---:|---:|
| g | 21 | 12 | 0.88 ± 0.43 |
| r | 48 | 23 | 0.85 ± 0.17 |
| i | 37 | 7 | 1.32 ± 0.09 |
| z | 50 | 5 | 1.27 ± 0.10 |

The DECam field is crowded. The nearest catalog entry, 0129m727.143190, lies 1.62 arcseconds away and is 3.98 r magnitudes fainter. Its complete disappearance could dim the combined light by about 0.027 mag if the aperture admitted all its flux, enough to account for the measured excursion. Only six of its own measurements passed the quality and comparison cuts. The independently fitted scales for neighbors at 1.91 and 2.96 arcseconds were 0.19 ± 0.46 and 0.08 ± 0.42. The target's PSF-fit centroid response was +0.33 ± 0.17 arcsec per magnitude east and −0.12 ± 0.09 north. Deblending, cadence and seeing affect these estimates. They favor the target image while leaving the closest companions unresolved as possible contributors. The small numbers of distinct nights in i and z also limit those bands' error estimates.

## Reproduce

Run `OPENBLAS_NUM_THREADS=1 python measure.py` from this directory. It uses the neighboring directory's OGLE-only fit, all DELVE rows associated with the local source IDs on target exposures, and the original field catalog. We omit forced photometry and select magnitudes from 5 to 24 with positive errors below 0.08 mag. Overlapping tile entries with the same exposure number and object count once, keeping the smallest uncertainty. For each object and filter, its median measurement sets its baseline. Each image uses the median offset of at least three of the nearest 15 comparison objects with r = 14–18.7 and separation above 8 arcseconds. We account for exposure midpoint and heliocentric light time, fit a separate offset and orbital template scale per band, and project the short waveform with a common scale. The release includes every unselected near-neighbor entry in the field table; missing catalog magnitudes are 99.99.

Photometry comes from the DELVE DR3 Magellanic Clouds release at NOIRLab Astro Data Lab. The files retain the selected exposure IDs, positions, original magnitude errors and exposure times. [Original queries and SHA-256 digests](sources.json) · [DELVE survey](https://delve-survey.github.io/)
