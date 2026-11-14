## KELT-South test design

Set before downloading epoch photometry, after inspecting the four archive metadata rows. Use the two TFA orientations as the primary data and repeat on the underlying raw curves. The duplicate reductions share exposures.

Keep the Gaia shape and phase and the previous TESS quadratic timing scale fixed at −0.80164096. Compare zero delay and the same magnitude with reversed sign. Fit an inner period correction and linear or quadratic phase drift as nuisance terms within KELT, with a separate Fourier waveform per orientation. Hold out entire nights deterministically across each survey season. Compare flux-space prediction error, give night-block uncertainty in paired differences, and report phase-shifted controls and the freely estimated scale as diagnostics. Do not interpret a null outcome as exclusion without a sensitivity estimate. Use the supplied archive barycentric times and report any clock uncertainty.
