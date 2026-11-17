# Uranus, JWST program 8975

**No new satellite passed two-visit validation in either filter.** The measured F162M recovery depth is **3 μJy (AB 22.71)** for observable circular, equatorial orbits from 79,500 to 80,500 km. At that flux, a separate native-detector test recovered 40 of 40 off-grid injections. Its exact two-sided 95% binomial interval is 91.2–100%. The field and masks admitted 80 of 120 uniformly drawn radius/phase combinations. Coverage, sampling error and systematic PSF uncertainty apply separately.

| Real observation pair | Grid templates | Eligible for two visits | Positive templates | Negative templates |
|---|---:|---:|---:|---:|
| F162M, October 6, 2025 | 63,720 | 10,970 | 10 | 4 |
| F140M, October 8, 2026 | 63,720 | 7,261 | 0 | 169 |

Each pair consists of four distinct on-planet CAL exposures in each visit, on eight short-wave detectors across both modules. F162M used visits 08975001001 and 08975003001. F140M used 08975005001 and 08975007001. The 500 km by 1° grid tested all phases at 52,000–140,000 km. Only 45 of 177 radius bins had any eligible templates in each band, spanning 54,500–86,000 km with internal gaps. Each filter has its own per-radius counts in [F162M](f162m_grid/radial_bins.csv) and [F140M](f140m_grid/radial_bins.csv). Neighboring templates share pixels.

The ten F162M positives all followed the same narrow radial feature at 55,000 km, across ten consecutive degrees of phase. The strongest template gave visit scores 6.41 and 7.24. Modeling the azimuthally extended halo with a median in 0.65-pixel-wide projected equatorial annuli reduced those scores to 2.88 and 4.41. All 720 half-degree phases then failed. Belinda, held out of the moon mask in a separate positive-control reduction, gave 8.99 and 16.08 before halo subtraction and 6.54 and 12.66 afterward. None of its 36 displaced comparison tracks passed. The residual fails the source validation criterion. [Ring test](f162m_ring_check.json) · [Belinda control](control_714.json) · [Original positive templates](f162m_grid/passing.csv)

## Measured completeness in the Perdita–Puck gap

The denser local search used five radii at 250 km intervals and 480 phases per radius at 0.75°. Other test particles had radii and phases drawn continuously and independently of this bank. We required at least two usable exposures in each visit. Observable paths were selected using only the science mask and template support. A matching detection had to pass the full local blind search and lie within 1.5 pixels of truth throughout both visits. All injections were made individually before glare modeling. A preexisting matching peak counted as a miss.

| Planted flux | AB magnitude | F162M recovered | F140M recovered |
|---:|---:|---:|---:|
| 0.1 μJy | 26.40 | — | 0/36 |
| 0.3 μJy | 25.21 | 0/40 | 0/36 |
| 1 μJy | 23.90 | 25/40 | 0/36 |
| 3 μJy | 22.71 | 40/40 | 0/36 |
| 10 μJy | 21.40 | 40/40 | — |

The F140M field admitted 55 of 120 random paths. Its second visit has substantially stronger structured detector noise, evident also in the 169 inverted-image passes. Its tested flux levels gave no two-visit completeness depth. In F162M, the 3 μJy flux corresponds to illustrative diameters **66, 47 and 33 km** for geometric albedos **0.04, 0.08 and 0.16** in the same band. The computation uses the CALSPEC solar spectrum, SVO F162M throughput, a phase factor of unity, and mean Sun/Uranus and JWST/Uranus distances of 19.505 and 18.814 AU. Diameter scales as the inverse square root of assumed phase factor and albedo. [F162M recovery](f162m_plants/report.json) · [F140M recovery](f140m_plants/report.json) · [Native injection record](f162m_native/report.json) · [Size calculation](f162m_size.json)

The native challenge independently integrated a subpixel PSF in each original detector frame, perturbed the planted position by Gaussian offsets of 0.3 native pixel RMS per coordinate per exposure, and transformed added flux through the SIP WCS and bilinear resampling. It used the linearity of the calibrated resampling operator to add that contribution to the real image. The search used the original orbital bank, threshold, variance and masks. Measured efficiencies are conditional on the chosen circular population, PSF, mask, and fixed calibrated noise. Both injection methods use the same underlying radial encircled-energy profile. They do not measure sensitivity for eccentric/inclined moons, other radii, unrecovered angular phases, or other filters.

## Reproduction

The two input manifests contain original MAST exposure filenames, SHA-256 hashes, detector centers, exposure timing and DQ policy. Use the public calibrated files, Python with Astropy, NumPy and SciPy, and the ephemerides described in [GEOMETRY.md](../GEOMETRY.md). The native science mask accepts DQ values 0 and 4 alone, and valid SCI and ERR. The known moons, rings and planetary disk are hard spatial exclusions. The search also requires 80% of the finite trailed PSF stamp on supported pixels. Limb fits supply corrections as large as 1.83 pixels in F162M. [F162M inputs](f162m_inputs.json) · [F140M inputs](f140m_inputs.json) · [F162M limb fits](f162m_centers.csv) · [F140M limb fits](f140m_centers.csv)

From the parent result directory, the following reproduces the F162M run after obtaining the ephemeris and individual CAL files. Use the distinct sets named in the manifest in place of the shell glob where the directory also contains off-target-calibration detectors. The preparation script applies the center table to select the on-planet detector for each exposure.

```sh
python limb_center.py centers.csv cache/jw0897500[13]*_cal.fits --radius-px 60
python prepare.py ephem.npz cube.npz cache/jw0897500[13]*_cal.fits --centers centers.csv --size 421 --samples 33 --psf-sigma .9 --mask-pad 5 --allow-corrected-jumps
python jwst8975/pilot_inputs.py cube.npz work
python jwst8975/broad_search.py work/science.npz broad
python jwst8975/coverage.py work/science.npz work/plants.npz work/covered.npz coverage.csv
python search.py work/science.npz work/tracks.npz measured --plants work/covered.npz --flux 3e-7 1e-6 3e-6 1e-5 --n-plants 40 --cache-templates --match-px 1.5 --psf-support .8
python jwst8975/native_plants.py work/science.npz work/tracks.npz work/covered.npz cube.npz.inputs.json cache work/Encircled_Energy_SW_ETCv2.txt native --n 40 --flux 3e-6
python jwst8975/reflected_size.py work/science.npz calibration size.json --band F162M --flux-jy 3e-6
```

For F140M change the exposure glob to observations 005 and 007, use 65 path samples and an initial sigma of 0.85, then inject 36 at each flux shown. In both filters the pilot replaces the initial Gaussian with a pixel-integrated radial model from STScI's encircled-energy table, retaining the measured fraction inside the finite stamp (0.844 for F162M). The score is inverse-variance trailed PSF photometry with broad clipped glare subtraction and an empirical radial noise correction. A passing template exceeds score 5 independently in two visits and 7 overall while keeping a score of 2 after removing the most influential exposure per visit. Nearby blank controls, inverted-image passes and the ring test quantify several sources of false positives. A grid score has no direct global false-alarm probability. [STScI NIRCam PSFs](https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-performance/nircam-point-spread-functions) · [NIRCam flux calibration](https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-performance/nircam-absolute-flux-calibration-and-zeropoints)
