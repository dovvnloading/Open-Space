# OGLE-2026-BLG-0058: follow-up advanced on 23 November 2026

**Disposition: candidate isolated black-hole lens; mass unmeasured.** The useful new result of this pass is an observing-time correction: propagating the deposited annual-parallax chains places the median maximum of the **unblended source-centroid deflection in spring 2027**, with broad uncertainty. The earlier July estimate used a straight-line approximation. Spring deserves priority for a feasibility assessment. A separate solar-elongation check places this target outside JWST's field of regard in July.

This packet is based on repository commit `972cd5169eea2b83d7407df84eb60c167b3c084e`. It adds reproducibility checks and forecasts; it contains **no new astronomical observations, measured angular Einstein radius, or confirmed black hole**. No observing request, paper, ATel or email was submitted in this pass.

## Question, observing unit and decisive tests

The question is whether a long microlensing event with a small fitted parallax is caused by a dark, high-mass lens, rather than a slowly moving ordinary star or a model/systematics alternative. The existing observing unit is an OGLE I-band photometry row, with HJD, magnitude, uncertainty, seeing and sky level. The retained data contain 1,900 rows from HJD 2457791.88602 to 2461345.51161. These are the same photons used in the existing analyses.

Two measurements would change the interpretation. First, independently calibrated 2027 photometry can test the frozen light-curve forecasts and the sensitivity of the parallax to the seven late-2026 points. Second, a joint fit of multi-epoch astrometry and photometry can measure the angular Einstein radius, source motion, lens–source relative parallax and blend contribution. A low inferred mass, detectable luminous lens, or a better source-orbit/blend explanation would weaken the black-hole interpretation. A sufficiently massive dark lens with these alternatives tested would strengthen it. See the existing [result and limitations](../../results/ogle-dark-lenses/README.md).

## Checks completed in this pass

| Check | Result | Evidence and limit |
| --- | --- | --- |
| Input integrity | All five selected science-input SHA-256 hashes match `sources.json` | Photometry, frozen forecast, saved chains, posterior summary and deposited independent fit; [checks.json](ogle-evidence/checks.json) |
| Independent fitting code replay | 1,900 points; rescaled parallax improvement Δχ² = 24.879393; positive and negative geometries recovered | [fits.json](ogle-evidence/fits.json), [stdout](ogle-evidence/fit-replay.txt). All three χ² values match the deposited output; the largest parameter difference is 0.00000585 day in tE. This is a numerical reproduction, not a new observation. |
| Baseline replay | Output matches the deposited baseline report exactly | [baseline-replay.txt](ogle-evidence/baseline-replay.txt). The actual 2017–2024 season-median range is **I = 18.112–18.129**; the old prose omitted the 2024 minimum. |
| Frozen forecast replay | All 24 parsed CSV rows match at the deposited 0.001-mag precision | [forecast-replay.csv](ogle-evidence/forecast-replay.csv), [stdout](ogle-evidence/forecast-replay.txt). No posterior was refitted or regenerated for this test. |
| Calendar audit | All eight frozen forecast dates have the correct JD-at-midnight labels; HJD 2461215.7 corresponds approximately to **24 June**, not 21 June 2026 | HJD-to-calendar comparisons here do not remove the heliocentric correction or change time scale. The several-minute distinction does not explain the previous two-day mismatch. The alert's unblended fit and the deposited parallax fit have different peak dates. |
| Astrometric timing | The frozen annual-parallax trajectories give spring median maximum dates | Details below. This revises a forecast, not the event's scientific disposition. |
| Public status check | The returned OGLE event page lists the same coordinates and a last epoch of HJD 2461345.51161, dated 1 November 2026 | [Official event page](https://ogle.astrouw.edu.pl/ogle4/ews/blg-0058.html). The web result can be cached; it does not establish that the live photometry bytes are unchanged. |
| Fresh photometry retrieval | **Not completed:** the photometry endpoint returned HTTP 503 | [Receipt](ogle-evidence/live-photometry-receipt.json). No updated photometry was substituted into the frozen fit. |
| Prior-report check | Exact-name web searches returned the OGLE/KMTNet event pages and no separate black-hole report in the returned results | [Search record](ogle-evidence/web-receipts.json). This is not exhaustive ADS, ATel, proposal, or unpublished-work clearance. |

The fit optimizer emitted warnings at singular trial flux solutions; the selected solutions reproduce the archived finite fit values. The no-parallax fit has negative blend flux and therefore an undefined blend magnitude (`NaN`), as in the archived output. This pass did not rerun the full 5,519-event fit screen, new MCMC chains, Galactic population integration, or all systematic tests. The remnant-model probability of 0.83–0.93 remains the existing model-dependent result; equal weighting of the geometries gives 0.6–0.7. Neither is an astrometric mass measurement.

## New calculation: when the source deflection is largest

For an isolated point source and dark point lens, the centroid shift relative to the unlensed source is

\[
\delta(u)=\theta_E\frac{u}{u^2+2},\qquad
\frac{d(\delta/\theta_E)}{du}=\frac{2-u^2}{(u^2+2)^2}.
\]

Thus a rising crossing of `u = sqrt(2)` gives the maximum `theta_E / (2 sqrt(2))`. The trajectory must include the same annual-parallax terms used in the light-curve fit. This is the source-centroid expression in [Sahu et al. (2022), §2.3 and equation 6](https://oro.open.ac.uk/84469/1/Sahu_2022_ApJ_933_83.pdf); an unresolved neighboring source requires an additional blend model.

The calculation uses **all 4,800 saved rows per geometry** in `posteriors/2026-BLG-0058-chains.npz`, without Galactic reweighting. For every sample the source starts below `sqrt(2)` at 1 November 2026 and has exactly one rising crossing and no falling crossing on the 800-day, one-day grid. Linear interpolation locates each crossing; `numpy.quantile(method='linear')` gives the intervals. The script retains the original geocentric reference epoch. [Code](ogle_checks.py) · [Full numerical output](ogle-evidence/checks.json)

An independent review of this pass caught a sampling problem before publication: applying a second `[::4]` stride to the stored array retained only two original walker trajectories. In `candidate.py`, emcee's step-major/walker-minor flattened array of 40 walkers was already saved with `[::5]`, so the stored file contains eight walkers. Removing the extra stride moved the positive-geometry median maximum from 28 April to 6 May and broadened the upper interval. All dates below use every stored row. The rows remain correlated, and the saved file does not contain the full original 40-walker ensemble; no new convergence claim is made. The original frozen photometric forecast is preserved for its future-data test and replayed with its original selection. A future posterior audit should regenerate/check the full ensemble separately and label any revised forecast explicitly.

| Photometric geometry | Median date of source-centroid maximum | Central 68% interval | Central 95% interval |
| --- | --- | --- | --- |
| u0 > 0 | 6 May 2027 | 29 March–13 June | 1 March–21 July |
| u0 < 0 | 11 April 2027 | 15 March–13 May | 21 February–18 June |

Dates are rounded calendar labels for the HJD calculation. The independent best-fit parameters instead give 16 May and 11 April: a best-fit date is not the median of the posterior-derived dates. These intervals are conditional on the original photometric model and error assumptions, and do not include new systematic or source-orbit models. They also describe the Earth-based trajectory; an observatory-specific forward model is required for actual observations.

For an **assumed** θE = 3.5 mas, the theoretical unblended maximum is 1.24 mas. Selected evaluations of the independently fitted annual-parallax trajectories are:

| Date | Source deflection for u0 > 0 / u0 < 0 (mas) | Same lensing term after weighting by the OGLE fitted source fraction (mas) |
| --- | ---: | ---: |
| 1 April 2027 | 1.229 / 1.237 | 0.890 / 0.995 |
| 15 August 2027 | 1.214 / 1.195 | 0.856 / 0.946 |
| 1 September 2027 | 1.205 / 1.184 | 0.847 / 0.936 |
| 1 September 2028 | 0.913 / 0.872 | 0.629 / 0.681 |

The last column is only a dilution illustration using `A Fs / (A Fs + Fb)`. It is **not a complete photocenter prediction**: the changing weight of an offset blend, different filter flux ratios, proper motion and source parallax also matter. Nor does the same assumed θE in both columns represent each geometry's Galactic mass posterior. In particular, fading back to photometric baseline in 2028 does not make the astrometric signal negligible.

## Observing feasibility and proposed measurement package

The target is RA **268.1587083°**, Dec **−31.7933611°** (J2000; 17:52:38.09, −31:47:36.1). An Astropy geocentric apparent-Sun calculation sampled at midnight UTC and applying the documented 85°–135° solar-elongation bounds gives approximate JWST planning windows of **15 March–5 May** and **6 August–27 September 2027**. Solar elongations on 1 and 15 July are 166.9° and 155.1°, outside those bounds. [STScI field-of-regard description](https://jwst-docs.stsci.edu/methods-and-roadmaps/jwst-moving-target-observations/jwst-moving-target-supporting-technical-information/moving-target-field-of-regard)

These original approximate dates are not a GTVT or APT result. The subsequent GTVT check and a date-alignment issue are documented below. Guide stars, target-specific detector placement, saturation, scheduling and other constraints remain unevaluated. Complete APT and instrument feasibility checks for this coordinate and the requested setup. STScI explicitly distinguishes a visibility estimate from schedulability in its [GTVT guidance](https://jwst-docs.stsci.edu/jwst-other-tools/jwst-target-visibility-tools/jwst-general-target-visibility-tool-help). HST and ground-based adaptive optics require their own feasibility checks; JWST restrictions should not be transferred to them.

### Subsequent GTVT 1.2 audit

A GTVT 1.2 run later in this pass used the captured live Horizons `JWST_merged` ephemeris, centered on the Sun (`DE441`), with LT-corrected ICRF vectors. Review found that this installed package's initializer replaces each position coordinate by `(coordinate.shift() - coordinate) * 1 + coordinate`: the previous row's vector, while retaining the current date. All 364 nonmissing output vectors match the preceding Horizons date within 2.33 × 10⁻¹⁰ km; the first row has missing positions. This is a one-sample shift in the calculation, not a new astronomical effect.

The [audit script](ogle_gtvt_audit.py) preserves the [unmodified output](ogle-evidence/gtvt-original.csv) and [Horizons response](ogle-evidence/gtvt-horizons-ephemeris.txt). In a separate object it joins the original vectors by exact `JDTDB`, then reruns the existing GTVT fixed-target methods. It makes no network request and does not patch the installed package. [Audit, versions and hashes](ogle-evidence/gtvt-audit.json) · [Date-aligned output](ogle-evidence/gtvt-date-aligned.csv)

| Calculation | First to last accepted spring sample | First to last accepted late-summer sample |
| --- | --- | --- |
| Unmodified GTVT 1.2 output | 16 March–5 May 2027; 51 samples | 8 August–28 September 2027; 52 samples |
| Explicitly date-aligned variant | 15 March–4 May 2027; 51 samples | 7 August–27 September 2027; 52 samples |

These are **midnight TDB sample labels**, not exact boundary instants or whole-day availability; midnight TDB is about 69 seconds earlier in UTC, with the preceding UTC calendar date. The audit records both scales. GTVT uses strict 84.8°–135° bounds; substituting 85°–135° produces the same daily windows for this target. An independent vector-dot-product calculation agrees with the corrected separations and all field-of-regard flags. All 31 July samples remain excluded, spanning solar elongations of 140.75°–167.07°. The corrected table is a documented variant using GTVT methods, not an unmodified official result. It supports spring and late-summer planning; APT, guide-star, ETC, instrument and scheduling checks remain pending.

The proposed request to a specialist team is:

1. **A spring 2027 astrometric sequence**, preferably with separated epochs across the accessible window, accompanied by two-color imaging to resolve or constrain the blend and identify suitable reference stars.
2. **A late-summer/autumn 2027 sequence** to measure evolution and curvature. The conditional source deflection is still close to its maximum during this period, so missing the spring maximum does not eliminate the target's value.
3. **Additional seasons through and beyond 2028**, with the full lensing trajectory retained in the fit. Their purpose is to separate lensing from the unlensed source's position, proper motion and annual parallax; a late image should not simply be declared the unlensed reference.
4. **An instrument-specific precision budget**, using source/neighbor colors, dithers, distortion calibration, repeatability and reference-star systematics. A per-epoch scale around 0.1–0.2 mas is a useful feasibility target for the 0.35–1.4 mas source-shift range, not a demonstrated sensitivity for this field. Exposure times cannot be assigned from OGLE I alone.

The public [JWST program 10476](https://www.stsci.edu/jwst-program-info/download/jwst/pdf/10476/) is a relevant design example for a different event: it combines multiple filters and repeated seasons with a stable reference-star field. It does not provide telescope time or validated exposure settings for OGLE-2026-BLG-0058.

## Frozen photometric test for 2027

| Date | No-parallax I | u0 > 0 I, central 95% interval | u0 < 0 I, central 95% interval |
| --- | ---: | ---: | ---: |
| 15 February | 17.885 | 17.934 [17.917, 17.950] | 17.945 [17.924, 17.966] |
| 1 April | 17.950 | 17.985 [17.973, 17.997] | 17.995 [17.979, 18.011] |
| 1 June | 18.010 | 18.028 [18.020, 18.037] | 18.035 [18.023, 18.045] |

These are unchanged model predictions from the deposited [forecast](../../results/ogle-dark-lenses/data/forecast.csv). They exclude the noise of a future exposure; the repeated no-parallax value is a best-fit curve, not zero uncertainty on a physical prediction. February is a forecasting date, not a guarantee of ground-based visibility or low airmass.

For the next data release, retain the original chains and evaluate the new calibrated points against them before refitting. Preserve time system, observing site, filter, uncertainty, seeing, sky level, image quality and calibration-star measurements. Use independent field stars to constrain zero-point and seeing terms, and include correlations between exposures from a night. A 0.01-mag level for a multi-night mean is a sensible calibration objective given the April separation of 0.035–0.045 mag; repeated points cannot average away a common zero-point error. Report predictive agreement or disagreement for both parallax solutions, then fit the combined data separately. An independent telescope is especially useful because reprocessing OGLE measurements does not remove shared observational systematics.

KMTNet's [2026 policy page](https://kmtnet.kasi.re.kr/ulens/event/2026/) confirms the event cross-identification and the general release rule of 1 July in the following year. It permits preparation using posted curves during the proprietary period but restricts submission, arXiv and public reproduction. This pass did not use KMTNet photometry. Confirm the applicable release status before including those measurements in a public joint analysis.

## What would count as a mass result

Use a joint astrometric/photometric posterior for `M = theta_E / (8.144 pi_E)`, including covariance, both geometries, blend/lens light and source motion. At the existing median parallaxes, a 3-solar-mass lens corresponds to θE = **0.94 mas** and **1.38 mas** for the positive and negative geometries. Substituting each existing one-sided 95% parallax upper value gives **1.19 mas** and **1.75 mas**. These are planning scales only: combining separate marginal bounds is not a joint confidence statement or a black-hole confirmation criterion.

The most useful rejection tests remain a slow luminous lens, source orbital motion (xallarap), chromatic or seeing-dependent blending, and astrometric reference-frame systematics. A coherent displacement measured in one epoch cannot by itself distinguish these from lensing. The next specialist review should assess both scientific alternatives and practical feasibility before any telescope request is finalized.

## Reproduce this pass

Requirements used for the photometric replay and original timing calculation: Python 3.12.14, NumPy/SciPy from the workspace runtime, Astropy 8.0.3, emcee 3.1.6, PyERFA 2.0.1.5 and astropy-iers-data. Exact NumPy/SciPy versions are recorded in [checks.json](ogle-evidence/checks.json). From the repository root, with these packages importable:

```sh
python results/ogle-dark-lenses/check/pspl.py results/ogle-dark-lenses/data/ews/2026/blg-0058.dat docs/followup-2026-11-23/ogle-evidence/fits.json
python results/ogle-dark-lenses/check/baseline.py results/ogle-dark-lenses/data/ews/2026/blg-0058.dat > docs/followup-2026-11-23/ogle-evidence/baseline-replay.txt
python results/ogle-dark-lenses/forecast.py results/ogle-dark-lenses/data/ews 2026-BLG-0058 results/ogle-dark-lenses/posteriors/2026-BLG-0058-chains.npz docs/followup-2026-11-23/ogle-evidence/forecast-replay.csv
python docs/followup-2026-11-23/ogle_checks.py
```

For the separate offline GTVT audit, use `jwst_gtvt==1.2` and the dependency versions recorded in [gtvt-audit.json](ogle-evidence/gtvt-audit.json), then run `python docs/followup-2026-11-23/ogle_gtvt_audit.py`. Its captured ephemeris and original CSV are included, so this replay does not query Horizons.

The fitting/forecast functions read frozen inputs and write separate replay outputs; the original posterior chains and predictions remain intact. The changed README, paper and message drafts correct the observing-time statement, calendar mismatch and baseline-range wording. Author review, complete literature checks, final observing feasibility, collaborator responses and any submission are still pending.
