# Held-out and adversarial test protocol: WISE J0458+6434 AB

**Scope:** follow-up to PR #41, not a reclassification of the original 2017/2020/2023/2025 and 144 MHz discoveries. Do not interpret a forced measurement, catalogue absence, or fitted Stokes V outlier as independent confirmation.

## Pre-specified competing hypotheses

1. **M (moving radio dwarf)**: the dominant 2–4 GHz emission follows the fixed infrared proper-motion trajectory. Each epoch has an unconstrained flux because bursts/variability are possible.
2. **S (static background)**: the same sightline contains one stationary, potentially variable synchrotron source. Fix its location independently at the **LoTSS DR3** source centroid to prevent refitting to VLASS.
3. **M+S (blend)**: a moving centimetre-wave emitter and stationary metre-wave/centimetre contaminant coexist. Fit two convolved restoring-beam templates simultaneously, with unrelated epoch fluxes; guard against near-collinearity.
4. **A (artefact/selection)**: the weak Quick-Look (QL) peaks or fitted positions are impacted by phase errors/ghosts or the selection of an apparent source from 1,526 trial ultracool-dwarf positions. A single field-offset rank is not a global FAP.

## Locked tests

- Use the actual four FITS epochs **separately**, actual observing dates, known IR trajectory, and 144 MHz catalog coordinate. Hold the 144 MHz coordinate fixed throughout the image fits.
- Refit native QL pixels in a *common per-epoch 19×19-pixel patch* with the known elliptical restoring beam, constant background and two background gradients: moving-only, LoTSS-stationary-only and moving+stationary. Report amplitudes, patch SSE, and separation at each date; do not treat pixel SSE as independent-noise chi-squared because synthesized-beam pixels are correlated.
- Repeat a negative control with co-translated template pairs in noise-only nearby patches to estimate **empirical** two-source amplitude uncertainties and covariance. Do not claim that a barely resolved pair is deblended in 2023.
- Compare the measured LoTSS fitted/deconvolved size and position angle with the **maximum** proper-motion smear possible across dated local input pointings. Account for uncertainty in unknown pointing weights, resolution and ionospheric smearing.
- Compute cross-band spectral-index *scenarios*, not a measured simultaneous spectral index, including a steep-spectrum static contribution at 3 GHz. Derive luminosity and brightness temperature under explicit emitting-area assumptions; do not infer ECM or field strength from Stokes I alone.
- Test constancy of the four GHz forced amplitudes, including the 2020 non-detection, while stating the pre-selection caveat. Recompute at multiple patch radii to check fitting-window dependence.
- Separate **no coverage**, **covered but too shallow**, **catalog nondetection**, and **unqueried archive**, and provide primary survey citations.

## Falsifiable follow-up

First priority: multi-hour, **phase-referenced, simultaneous** VLA S-band and L-band or C-band full-Stokes imaging and a LOFAR HBA time-resolved Stokes-I/V observation, with sub-arcsecond centroid comparisons to the propagated IR coordinate and one or more stationary in-field astrometric controls. A 144 MHz detection fixed on the 2023 LOFAR location while the GHz source continues northeast would falsify a single common radio emitter, even if the GHz moving source is genuinely the dwarf. A fixed GHz source would falsify the dominant moving GHz interpretation. Repeat low-frequency observations in at least two well separated dates; the 6-arcsec survey mosaic alone cannot measure 0.37 arcsec/year motion.

## Provenance

Inputs remain in `results/wise-j0458-radio/` and `results/wise-j0458-radio/independent-radio/` as committed in PR #41; this follow-up only adds files. Relevant authoritative descriptions: NRAO VLASS QL guide https://science.nrao.edu/vlass/data-access/vlass-epoch-1-quick-look-users-guide ; LoTSS DR3 https://lofar-surveys.org/dr3.html ; LoTSS DR2 https://lofar-surveys.org/dr2_release.html ; TGSS ADR1 https://doi.org/10.1051/0004-6361/201628536 .
