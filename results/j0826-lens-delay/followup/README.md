# J0826+7002: independent time-delay fits and tests that do **not** yet confirm a precision delay

This is a **new, separately reproducible follow-up** to [the original ZTF–Gaia echo study](../README.md). Only files in `followup/` and dedicated validation workflows were added. The original light curves, analysis code, and other findings were not changed.

## What survives, and what does not

**The leading hypothesis remains that northern quasar image B precedes southern image A by approximately seven months, not by a precise measured 210.0 days.** A joint Gaussian-process (GP) fit with intrinsic quasar variability shared between both images and separate smooth, image-specific extrinsic processes has a favorable delay near 210–230 days. A common intrinsic cubic-spline model with separate slow microlensing also favors a delayed solution over zero for some predeclared roughness and knot spacings. The signal is **not model-independent**: changing intrinsic spline flexibility or the extrinsic timescale changes the preferred delay, and later held-out observations do not reliably prefer the exact 210-day hypothesis to zero delay.

The established lens has source redshift **z=1.618**, two spectroscopically identified quasar images separated by **5.83 arcsec**, and multiple candidate foreground lens galaxies of **unknown secure redshift**. This precludes a precision time-delay distance, mass or Hubble-constant measurement. [Discovery spectroscopy and field description](https://academic.oup.com/mnras/article/520/3/3305/6948349).

### A joint, undetrended GP rather than high-pass correlation

[`gaussian_process.py`](gaussian_process.py) uses 14-day disjoint medians of the *original* ZTF g-band nights. An OU intrinsic process (q(t)) is shared between images, with image A evaluated at (t-\Delta), while each image has an independent long-timescale squared-exponential extrinsic process. The two mean magnitudes are integrated out through generalized least squares, and each photometric error also receives an independently fitted variance floor. No prewhitening is applied. The independent-images comparator deletes the cross-image intrinsic covariance. [`robustness.py`](robustness.py) varies 10-, 14-, 20- and 30-day aggregation and microlensing timescales of 180, 400, 800 and 1,200 days.

A representative 14-day fit gives a fixed-nuisance profile minimum at **about 220 days**, with a strong conditional preference for a delay over zero and one-year aliases. This *profile is not a Bayesian posterior*: nuisance parameters were optimized conditional on a chosen delay; alternative noise priors, correlated calibration errors and flexible microlensing can change the result. In particular, 14-day bins are not 239 independent points of information about the delay.

### Future-epoch holdouts are mixed, not universal confirmation

For a stringent split, [`gaussian_process.py`](gaussian_process.py) learns the nuisance hyperparameters on data **before** a specified MJD, using the 210-day hypothesis, and scores all later measurements by their conditional predictive density without refitting the nuisance parameters or means. The below differences are `held-out NLL at zero delay minus held-out NLL at +210 days`; positive values favor the delayed model. They are model-specific *predictive score differences*, not significance levels or Bayes factors.

| Training ends before MJD | ZTF 14-day bins trained/tested | NLL(0) − NLL(210) | Interpretation |
| --- | ---: | ---: | --- |
| 59600 | 129 / 110 | **+21.14** | The locked early model prefers the delayed hypothesis |
| 60000 | 163 / 76 | **+5.50** | Reduced advantage |
| 60200 | 171 / 68 | **−0.26** | The most recent held-out block does **not** favor exactly +210 over zero |

At MJD 60200, nearby lags (+215 to +235 days) perform similarly or modestly better than +210, so this tests the precision and stability of the proposed lag rather than conclusively rejecting a delayed echo.

The separate [`seasonal_template.py`](seasonal_template.py) uses a lagged nonparametric template of B and a slow spline for A's extrinsic variation, with each A observing season held out. It compares competing lags on **identical A test nights**, but permits the B template to use *all* B years, including B exposures temporally corresponding to the held-out A season. This is **cross-image reconstruction, not a prospective forecast**. With eight seasons and 418 common A nights, pooled RMS errors are:

| Slow extrinsic spline knot spacing | +210 d | 0 d | +365 d |
| --- | ---: | ---: | ---: |
| 350 d | 0.2148 mag | 0.2664 mag | 0.2400 mag |
| 550 d | 0.1836 mag | 0.1918 mag | 0.2331 mag |
| 850 d | 0.1636 mag | **0.1515 mag** | 0.2421 mag |

The zero-delay template sometimes obtains a **negative fitted response coefficient**, which cannot represent positive macro-lens magnification. But this unconstrained comparator remains a useful rival *statistical* explanation for noisy correlated trends. The reversal when microlensing is smoothed more strongly demonstrates that an apparent narrow optimum would understate systematic uncertainty.

### Joint splines: a serious model-flexibility check

[`free_knots.py`](free_knots.py) models one shared intrinsic cubic spline and an independent slowly varying A-only extrinsic spline, using identical numbers of coefficients and roughness regularizers at each lag. It sweeps intrinsic knot separations of 90–280 days and extrinsic separations of 450–1,000 days, and in a separate test allows the intrinsic knots to move ±22% of their nominal spacing (at most 50 days). Depending on the knot spacing, some fits prefer delays roughly **260–300 days** rather than 210. This disagreement means **do not claim a calibrated ±5-day, ±15-day or 95% confidence interval**. Current defensible reporting is *a candidate positive delay of order 200–300 days, with non-negligible aliases and extrinsic-model uncertainty*, not a precision cosmographic measurement.

### Attempted third, wholly independent survey: not a confirmation

The public [PTF light-curve catalog](https://irsa.ipac.caltech.edu/Missions/ptf.html) does cover this northern sky field. But the [27 extracted PTF epoch records](ptf-epoch-query.csv) consist of **17 measurements assigned to one object near bright image A and 10 to a different nearby field object**. **None provide a separately resolved B light curve.** All have `photcalflag=0`. Thus PTF cannot yield a third telescope-based time-delay estimate. The frozen query, exact IDs and extraction/provenance checks live in [`ptf_probe.py`](ptf_probe.py), [`ptf_direct.py`](ptf_direct.py) and [`test_followup.py`](test_followup.py). The parent study's Pan-STARRS1 DR2 detections have just **one same-filter epoch pair within 12 days** of +210 days across the two images, and pooled disparate passbands actually favor an unrelated alias. They must not be cited as third confirmation.

The independent cross-era observational support is still **ZTF and Gaia**, not three surveys. A forced-PSF extraction from actual PTF images at **both** published image positions, with image-blending tests and reliable calibration, would be a distinct project and must not be described as completed here.

### What the delay would imply *if confirmed*: SIS sanity check only

[`sis_sanity.py`](sis_sanity.py) adopts the *published* quasar redshift, two-image 5.83 arcsec separation and a toy singular-isothermal-sphere lens in an assumed flat `H0=70 km/s/Mpc, Ωm=0.3` cosmology. It sweeps **assumed**, not measured, lens-galaxy redshifts. At an **illustrative z_l=0.5**, the SIS has an Einstein radius of **2.915 arcsec** and a maximum double-image delay scale of **1,530.5 days**. Producing 210 days corresponds to `β/θ_E ≈ 0.137` and a toy `σ ≈ 415 km/s`. The SIS's predicted macro flux ratio is only about **1.32**, versus a Gaia G-band A/B flux ratio near **6.1**. This discrepancy is unsurprising when multiple foreground galaxies, a possible blended third quasar image, galaxy-light contamination and chromatic microlensing are ignored.

**These numbers are conditional consistency checks, not measurements of lens redshift, mass, dispersion or H0.** The unknown deflector redshift and mass-sheet/external-shear degeneracies prevent robust inversion of the delay.

## Specific falsification tests and research status

A persuasive confirmation needs (1) independently reduced, contemporaneous and *separately resolved* multi-band light curves, (2) multiple distinct intrinsic quasar features seen first at B and then in A at a stable, common delay, (3) out-of-era prediction with lag and nuisance priors frozen before future measurements, and (4) model-averaged uncertainty that includes extrinsic variability, instrumental correlations and annual sampling aliases.

The 210-day interpretation would fail if sufficiently sampled new light curves consistently require **zero, negative or annual-alias offsets**, if B's well-resolved flares do **not** recur in A in the predicted window, or if realistic microlensing/survey-calibration models reproduce the apparent echo **without any shared delayed signal** and outperform it on independent held-out epochs. At present, the last holdout is **not positive evidence for the exact 210-day model**, and a third successful independent survey is **missing**.

### Reproduce and limits

Run on Python 3.10 or 3.12 with NumPy and SciPy:
```sh
python -m unittest discover -s results/j0826-lens-delay/followup -p 'test_followup.py' -v
python results/j0826-lens-delay/followup/gaussian_process.py --quick
python results/j0826-lens-delay/followup/robustness.py --quick
python results/j0826-lens-delay/followup/free_knots.py --quick
python results/j0826-lens-delay/followup/seasonal_template.py --quick
python results/j0826-lens-delay/followup/sis_sanity.py
```
All input ZTF/Gaia/PS1 nights in this study come directly from the immutable upstream tables in the parent directory. The dedicated [GitHub J0826 numerical workflow](../../../.github/workflows/j0826-model-validation.yml) runs `test_followup.py` on two Python versions. The PTF audit workflow retains the original public retrievals in a short-lived Actions artifact; [`ptf-epoch-query.csv`](ptf-epoch-query.csv) is the frozen inspectable source extract. No uncertainty interval has been calibrated by injection-recovery or a fully specified hierarchical posterior; numerical precision in any machine output does not change this limitation.
