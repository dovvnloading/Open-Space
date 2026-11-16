# Research checkpoint: SPHEREx moving-object spectra

**Bottom line:** a target-centered spectral-excess lead in comet C/2025 R2 (SWAN) repeats on December 22, December 25 and a later January 2, 2026 sequence. This is not yet a validated molecular identification, new discovery or unusual-composition result. Calibration-aware testing exposed limitations that prevent promoting it as one.

## Review corrections: read before interpreting fit statistics

**January is not a pristine holdout for the later model-family comparisons.** Its Gaussian-model failure was inspected before CO2, Lorentzian and triangle comparisons were developed. Although the CO2 temperature grid was selected using December values only, subsequent family comparison adaptively reused January. Call these later-epoch shape-transfer checks, not independent held-out validation. Earlier report wording and JSON keys containing “heldout” record the historical workflow and must be read with this correction. A genuinely untouched epoch or independent instrument is still required.

**Spectral weighting convention remains unresolved.** The original scripts assume photon-weighted flat-Fnu calibration, using integral(R F_lambda lambda)/integral(R/lambda). The official explanatory supplement v2.1 Eq. 22 instead states Fnu_pred = integral(t S_k)/integral(t). Under the interpretation that delivered NormTransmission equals t and S_k is Fnu, convert the PSG model with Fnu proportional to lambda² F_lambda before averaging. An independent reviewer obtained December chi-square 6.910 and January 4.749 for the same 100 K model, versus 6.473/4.879 under the original weighting; the four-temperature ranking survives. Numerical stability between these interpretations does not establish which convention is correct. Original results below remain conditional, not calibrated measurements.

The [official supplement](https://irsa.ipac.caltech.edu/data/SPHEREx/docs/SPHEREx_Expsupp.pdf), [R7 spectral product guide](https://irsa.ipac.caltech.edu/data/SPHEREx/docs/SPHERExSpectralFiles_202607.pdf), and [Hui et al. calibration paper](https://arxiv.org/abs/2602.09139) were checked. Their public descriptions do not resolve the mapping sufficiently for this analysis. No discovery, production rate, molecular abundance or temperature inference is warranted.

## Completed

| Target | Selected records | Retrieved | Interpretation |
|---|---:|---:|---|
| 29P | 90 | 82 | Known positive control, not novel chemistry |
| 174P/Echeclus | 8 | 5 | Insufficient pilot sensitivity, not gas absence |
| C/2024 E1 (Wierzchos) | 15 | 9 | Known JWST CO2 comparison |
| C/2025 R2 (SWAN) | 57 | 56 | Repeat spectral-excess lead needing further validation |

Every missing download was retained as an HTTP failure, not a nondetection. For SWAN, 27 strict-mask and 46 OUTLIER-inclusive diagnostic apertures passed at radius 2.5 pixels. The first nine records were exploratory; the remaining 48 were selected before inspection. Full JSON, final CSV, source URLs and hashes are preserved.

The work includes official spectral-response calibration, aperture/background sensitivity, midpoint agreement within 44 microseconds, target-position and ePSF checks, AllWISE contamination checks, and empirical forward modeling. Four focused unit tests pass, and the scientific scripts compile. No browser suite or whole-repository CI was invoked.

## What survived, and what failed

The December 25 excess survives three background annuli and is centered close to the ephemeris. January 2 separately samples the rising feature again. Simple smooth continua cannot fit these measurements under the formal noise assumptions.

However, the six-point December Gaussian fit has six parameters and zero degrees of freedom. Its frozen shape fails to transfer to January 2: chi-square 44.8 for two degrees of freedom, or 91.6 with a larger aperture. Adding an assumed 10% uncertainty floor changes the answer, but that floor is not a measured systematic. A single 4.260-micron delta line also fails to reproduce the red wing. These failed approximations motivated the subsequent physical-template tests below; they do not themselves disprove a molecular band.

Point-source morphology is inadequate in bright feature channels. Gaussian broadening helps, but the continuum also broadens; this is not evidence of resolved gas. October has approximately three pixels of motion trailing and strong contamination, so it cannot be compared naively against December. September 11, 2026 measurements are weak, but distance and aperture geometry changes prevent interpreting that alone as physical fading.

All December pilot feature apertures require explicitly permitting the OUTLIER flag. Another core flag rejects the recovered December 22 blue-side image even under that diagnostic mask. No other bad flags were overridden. Historical catalog neighbors and bright-star wings remain deblending concerns; broadband catalog flux is not contemporaneous narrowband photometry.

## Subsequent physical-template and field-control tests

A documented NASA PSG non-LTE CO2 fluorescence template at 100 K, selected among 20/60/100/200 K using December data, fits the small-aperture January 2, 2026 sequence much better: chi-square **4.879 for two residual degrees of freedom**, with amplitude and constant refitted. December chi-square is 6.47 with two residual degrees of freedom before accounting for temperature-grid selection. This is shape compatibility, not a temperature measurement, production-rate retrieval or unique chemical identification.

Spatially weighted detector responses using point and steady-coma profiles leave January chi-square at 4.86–4.90. The larger 3.5-pixel aperture remains in tension at about 8.35/2. Particular fixed CO and H2O templates, and December-fitted Gaussian, Lorentzian and triangle shapes extrapolate worse, but the alternatives are not exhaustive. Shared radiometric assumptions and source-dependent systematic errors remain.

Translated-track field controls give January a studentized contrast exceeding all 65 valid controls (descriptive upper rank 1/66, not a global discovery significance). Their robust scatter is about 1.10 times formal noise. September 11, 2026 is an ordinary faint control, with 16 of 44 control contrasts exceeding the target. Neither test calibrates bright-source model errors.

PSF-convolved rho^-1 profiles improve the bright images, with feature amplitudes reasonably stable over fitted radii. Continuum profiles improve too, and residual structure persists; this does not demonstrate a gas coma distinct from dust. Approximate light-weighted response integration is not a full response-before-optical-convolution or radiometric retrieval.

New primary reports: `PSG_PHYSICAL_SHAPE_REPORT.md`, `COMA_PROFILE_TEST.md`, `EXTENDED_RESPONSE_AND_ALTERNATIVES.md`, `ALTERNATIVES_AND_INDEPENDENT_DATA.md`, and `FIELD_NOISE_AND_PSF_TEST.md`. Exact input configurations, spectra, hashes and scripts are preserved.

## Required next steps

Advance the tested generic fluorescence/approximate spatial-response models into a jointly calibrated coma, trail and neighboring-star retrieval; calibrate source-dependent uncertainties through injection/recovery; and seek independent observations. No independently confirmed contemporaneous infrared spectrum was retrieved. LCO optical observing activity in the relevant season is documented, but the archive query did not return calibrated frames. The bounded literature search found no indexed SWAN CO2/SPHEREx paper, but that does not prove priority. Do not claim unusual CO2/H2O, production rates or a first detection.

Review `SWAN_CANDIDATE.md`, `SWAN_SYSTEMATICS.md`, `RESPONSE_VALIDATION.md`, `EMPIRICAL_FIT_REPORT.md` and the independent spatial reviews. Pilot files are deliberately retained intermediate checkpoints, not the final 57-record sample.
