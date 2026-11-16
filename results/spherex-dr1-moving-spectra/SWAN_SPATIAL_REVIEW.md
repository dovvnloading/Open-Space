# SWAN independent spatial-code review

This is an independent implementation on the SAME image evidence, not independent observations. Script and JSON accompany this note. Exact FITS URLs/hashes retained. No repository files were modified.

Used nearest detector-position EPSF row, oversampling5 and 33x33 support. Fitted only radius2.5 native pixels (amplitude + constant background + centroid offset bounded0.5px). Outlier bit19 explicitly permitted; other invalid bits excluded. This is a spatial diagnostic, not calibrated comet flux. A finite EPSF cannot model extended coma, and correlated variance/calibration uncertainty is not included in chi-square.

- Dec25 4.2517um feature is centered within0.043 pixel by fitted ePSF centroid, but point-source model is rejected by formal reduced chi-square6.70 (15 dof).
- Same sequence 4.2155 continuum has chi-square1.46 (17 dof); 4.2882 feature-side has3.62 (15 dof). These differences merit coma/extended morphology modelling, not a claim of resolved gas.
- Dec22 4.2928 feature-side has chi-square4.90 (16 dof) and fitted offset0.306px. Continuum4.3294 has1.03 (15 dof).
- Invalid core of Dec22 4.2566 remains rejected; do not relax bit28.
- October4.0431 image is NOT pointlike at the predicted center: fit hits0.5px bound, chi-square~846, positive centroid~1.09px away, consistent with strong trail/extended structure and possible field contamination. Do not compare October raw flux to December as a physical change without motion-aware modelling.

Next falsification: fit convolved coma/trail profile with spatially varying background; recover feature/continuum contrast under common morphology and synthetic source injection. Check stationary catalogs and source-following centroid independently. All December diagnostic signals depend on permitting OUTLIER flags. No discovery claim is established.
