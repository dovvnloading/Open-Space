# Research checkpoint: SPHEREx moving-object spectra

**Bottom line:** a target-centered spectral-excess lead in comet C/2025 R2 (SWAN) repeats on December 22, December 25 and a held-out January 2 sequence. This is not yet a validated molecular identification, new discovery or unusual-composition result. Calibration-aware testing exposed limitations that prevent promoting it as one.

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

The December 25 excess survives three background annuli and is centered close to the ephemeris. January 2 independently samples the rising feature again. Simple smooth continua cannot fit these measurements under the formal noise assumptions.

However, the six-point December Gaussian fit has six parameters and zero degrees of freedom. Its frozen shape fails to transfer to January 2: chi-square 44.8 for two degrees of freedom, or 91.6 with a larger aperture. Adding an assumed 10% uncertainty floor changes the answer, but that floor is not a measured systematic. A single 4.260-micron delta line also fails to reproduce the red wing. Proper molecular-band and extended-source modeling is required before identifying the feature.

Point-source morphology is inadequate in bright feature channels. Gaussian broadening helps, but the continuum also broadens; this is not evidence of resolved gas. October has approximately three pixels of motion trailing and strong contamination, so it cannot be compared naively against December. September measurements are weak, but distance and aperture geometry changes prevent interpreting that alone as physical fading.

All December pilot feature apertures require explicitly permitting the OUTLIER flag. Another core flag rejects the recovered December 22 blue-side image even under that diagnostic mask. No other bad flags were overridden. Historical catalog neighbors and bright-star wings remain deblending concerns; broadband catalog flux is not contemporaneous narrowband photometry.

## Required next steps

Use a physical fluorescence template through the aperture-weighted response; jointly model coma, trail and stars; calibrate uncertainties through injection/recovery; and seek independent observations. The bounded literature search found no indexed SWAN CO2/SPHEREx paper, but that does not prove priority. Do not claim unusual CO2/H2O, production rates or a first detection.

Review `SWAN_CANDIDATE.md`, `SWAN_SYSTEMATICS.md`, `RESPONSE_VALIDATION.md`, `EMPIRICAL_FIT_REPORT.md` and the independent spatial reviews. Pilot files are deliberately retained intermediate checkpoints, not the final 57-record sample.
