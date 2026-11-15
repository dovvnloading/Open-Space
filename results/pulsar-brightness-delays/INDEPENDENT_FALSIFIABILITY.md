# Follow-up P: B1822−09 blind-frequency challenge and frozen falsifiers

**New additive result, 14 November 2026.** This follow-up analyzes arrival-time channels **not used to fit the original PR #21 chromatic slope**, while preserving the original 52-epoch CSV and analysis verbatim. It does not claim an independent telescope discovery or identify the physical mechanism.

## New observation: the channel holdout works

The base result fits a brightness-dependent effective DM using receiver pairs 0–3 and 4–7. Under a purely frequency-independent varying cold-plasma column, the effective DM slope against simultaneous flux density should be *the same for every pair* after a shared time trend. The original data found −0.0118507 and −0.0245615 pc cm⁻³/mJy respectively.

The **new falsifiable prediction** is obtained before looking at the other four arrival-time channels. If arrival time is locally described by

```
arrival_time(nu, epoch) = achromatic(epoch) + K * [A(epoch)*nu^(-2) + B(epoch)*nu^(-4)]
effective_DM(i,j) = A + B*(nu_i^(-2) + nu_j^(-2))
K = 4148.808 s MHz^2 / (pc cm^(-3)); nu is in MHz
```

then **per visit** the two original outer-subband pairs determine A and B. For each untested pair, calculate its actual `nu_i^-2+nu_j^-2` from that visit's channel centers, interpolate between pairs 0–3 and 4–7, then regress predicted and actual effective DM against simultaneous flux with a separate linear time trend. The pairing is frozen; no coefficients are optimized to the cross-band values. A constant receiver-band template offset can change intercepts, but not the tested flux slopes.

| Pair (0=low frequency, 7=high) | Measured | Predicted from 0–3 and 4–7 | Absolute residual | Residual / measured |
| --- | ---: | ---: | ---: | ---: |
| 0–7 (shares endpoints) | −0.0156874 | −0.0160917 | +0.0004043 | 2.58% |
| **1–6 (disjoint TOA channels)** | **−0.0175098** | **−0.0182221** | **+0.0007123** | **4.07%** |
| **2–5 (disjoint TOA channels)** | **−0.0182725** | **−0.0196329** | **+0.0013604** | **7.45%** |

All DM/flux slopes are pc cm⁻³/mJy. The independent-channel predictions are much closer to the measurements than a constant effective-DM/flux slope. The five estimates also have an ordered progression across the band. A two-term *smooth spectral model* passes this cross-frequency holdout. **It does not prove that the extra term is plasma dispersion:** intrinsic and instrumental profile evolution can produce the same function.

## Robustness tests performed on this fixed sample

The new script makes its numerical requirements executable, using `visits.csv` *as published in PR #21* and asserting its Git blob hash before analysis.

- **Four-visit moving circular-block bootstrap** (4,000 samples, seed 1715) on the original full-visit rows. All three cross-band prediction-error slope intervals include zero. The original upper-minus-lower, a new inner-minus-outer (1–6 minus 0–7), and center-minus-outer (2–5 minus 0–7) intervals are negative. The inner-minus-outer effect is **−0.0018224 pc cm⁻³/mJy**, with approximately **[−0.00249, −0.00126]** in a separately seeded equivalent block experiment.
- **Exact circular-shift placebo:** shift the 52-point contemporaneous flux series by each of its 51 nonzero visit offsets without changing the measured delays or time coordinates. **0 of 51** shuffled alignments equals or exceeds the full original upper-minus-lower, inner-minus-outer, or center-minus-outer slope in absolute value. This is a stationarity-sensitive diagnostic, not a discovery p-value, because the pulsar was originally selected using these 52 visits.
- **Time halves:** the new inner-minus-outer slope is −0.001881 in the first 26 visits and −0.001798 in the last 26. The original upper-minus-lower effect is −0.01104 and −0.01369 in the two halves. The signs are stable.
- **Early→late frozen extrapolation:** train only on channels 0–3 and 4–7 in the *first* 26 visits. Predict the later 26 visits for untrained pairs. For 1–6, the frozen prediction is **−0.0164403** versus measured **−0.0182562** (9.95% difference); for 2–5 it is **−0.0176629** versus **−0.0191101** (7.57% difference). The law approximately transfers, but a stationary amplitude is **not exact**. Both checks remain retrospective: the original source was selected using all epochs.
- **Flux extremes:** after removing the five brightest and five faintest visits, the inner-minus-outer slope remains negative (−0.0023024); the original chromatic difference remains negative (−0.01800).
- **Instrument proxies:** adding duration, annual sine/cosine, L-band midpoint and bandwidth as nuisance regressors retains a negative inner-minus-outer slope (approximately −0.00195). Remove any one visit and the inner-minus-outer slope remains negative (about −0.00203 to −0.00169).
- **Cross-check in another facility:** the three Parkes single-pulse files and their broad individual bootstrap intervals were already examined in PR #21. They are not a successful independent replication and are not rebranded here.

The `cross_frequency_holdout.py` command computes all these values, prints JSON, and aborts if the predeclared checks do not survive. The inter-epoch bootstrap is **conditional on selection**; its negative confidence intervals do not correct for 476-object discovery screening, multi-metric exploration or long-range autocorrelation.

## Mechanisms that remain viable

1. **True time-variable cold unmagnetized plasma column only** predicts a common effective DM across all five frequency pairs. The strong difference between band halves *rejects that complete explanation for the measured delay under the current timing template*, not the possibility that plasma supplies part of it.
2. **Pulse-profile state mixing / mode switching** is a leading explanation. B1822−09 has established correlated changes among the main pulse, precursor and interpulse, and extensive pulse-shape variation. A fixed template's fiducial point can slide with mode fraction and frequency, creating apparent DM and brightness changes without a corresponding change in electron column. [Single-pulse/mode morphology](https://doi.org/10.1051/0004-6361/202451384); [long-term spin-down and profile variability](https://doi.org/10.1093/mnras/staf427).
3. **Chromatic propagation and scattering** can cause genuinely frequency-dependent arrival times. Multipath sampling of inhomogeneous electron columns causes frequency-dependent apparent DM; pulse broadening often has a steep chromatic behavior. A `nu^-4` expansion is a fitting basis, **not** a uniquely inferred scattering exponent. [Cordes, Shannon & Stinebring 2016](https://doi.org/10.3847/0004-637x/817/1/16); [historical profile-evolution degeneracy](https://academic.oup.com/mnras/article/377/2/677/1036774).
4. **Calibration, flux selection and radiometer/template bias**: a higher S/N visit can change TOA uncertainty, frequency-dependent centroiding or fitted-mode weights. Shared calibration, correlated channel errors, RFI masking, imperfect clock/barycentering, narrow-band scintillation and using only retained detections remain available explanations. The two held-out pairs are disjoint in TOA **channels**, not in visit, observing system, source or reduction pipeline.
5. **Pure linear time trend / seasonal or frequency-center drift**: the supplied measured trends, chronological halves and adjusted nuisance regressions reduce these particular explanations. They do not rule out all nonlinear red noise, chromatic shape covariates or long-term instrumental drift.
6. **One epoch drives the result**: leave-one-epoch-out and 10-visit end trimming preserve the sign. The magnitude changes under trimming, so numerical amplitude should not be treated as fixed physical law.
7. **A fitted high-order dispersion term is unique**: false. Other smooth chromatic functions can interpolate these five non-independent pair estimates, especially over the narrow L-band range. The cross-frequency test supports *predictive smoothness*, not identification of a physical microscopic index.

## Preregistered falsification and next observations

- **Current-sample null `frequency-independent DM`:** all effective-DM slopes must agree if a single cold-plasma column is sufficient. The negative band-half difference and independent channel gradients reject this strict null *conditional on the arrival-time reduction*.
- **Current-sample new spectral prediction:** fail if a disjoint pair 1–6 or 2–5 has a frequency-model residual slope whose 4-visit-block 95% interval excludes zero, or if its absolute residual exceeds 20% of its observed slope. Both pass now; do not tune these thresholds after seeing new observations.
- **Chronological transport stress:** early-only training should predict a later disjoint pair to within 12% of the observed later slope. Both pass narrowly (9.95%, 7.57%), not a calibration-free physical validation.
- **Prospective independent replication:** acquire *new* wideband baseband or folded pulse portraits, preserve all complete visits, blind the later dates and frequencies, derive TOAs with state-specific frequency-evolving profile templates, then independently fit flux slopes in separated frequency bands. A chromatic slope difference consistent with zero with tight uncertainty, or vanishing after state-specific template correction, would **falsify the persistence of the timing effect under a profile-controlled model** while favoring a template-induced origin in the original data. A persistent positive/negative spectral pattern after controlled profile correction and matched calibrators would challenge a profile-only explanation and strengthen a propagation hypothesis.
- **Mechanism-specific preregistration:** profile hypothesis predicts covariance with independently inferred mode weights and a substantial reduction in fixed-template bias when mode-conditioned templates are used. Scattering hypothesis predicts broadening/asymmetric delay and chromatic scattering proxy that tracks the extra delay. A physically changing common electron column predicts `nu^-2` timing with a consistent DM in every independent subband. Measure these *before* comparing the test set.
- **Controls:** observe a stable nearby pulsar and stable timing/calibration source with matched backend, observing cadence and gain/RFI state; split tuning into training and a genuinely withheld new campaign. Recompute a whole-survey multiple-testing null rather than interpreting a post-selection single-pulsar p-value.

## Reproduce with no downloads

Run from repo root:

```sh
python results/pulsar-brightness-delays/cross_frequency_holdout.py > /tmp/B1822-heldout.json
```

Only Python standard library is required. The program verifies the unchanged PR #21 `visits.csv` by Git blob hash, reproduces all 52 epochs and tests, and exits nonzero if the results fail. The original PR #21 result files are untouched. This new cross-frequency prediction is **a new analysis of archival data**, not a claim of peer-reviewed discovery or independent observational priority.
