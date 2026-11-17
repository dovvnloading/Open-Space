# Blind spectral selection for tile 7119

Select every WISE-eligible Gaia dwarf proxy with W1<12 whose epoch-2025.9 proper-motion position is in order-5 tile 7119. The cohort uses astrometry and the W1 quality anchor without selecting W3/W4 color. Record attempted and usable denominators separately.

Use the DR1 cube, 19″ source aperture, 36–60″ sky annulus, identical blank apertures, and quality masks in the sibling spectra analysis. Test the ratio of median Fν over 4.42–4.90 µm to 1.51–1.72 µm, with 3.16–3.62 µm as context. Require at least half the nominal channels in each tested band (minimum two). Compare each star to up to 15 other usable stars in this tile within 0.5 mag BP−RP and 1.5 mag G. Require ≥4 controls and a positive median ratio for the ensemble. Use 1.4826×MAD of the ratios with a 5% photospheric scatter floor and finite-sample median uncertainty. Add aperture noise measured from ≥8 corresponding sky positions. Require H and long-band S/N≥5, ratio ≥1.25 times the peer median, and ≥5σ excess. A strict-mask trigger must survive a repeat including OUTLIER-flagged pixels before interpretation. Report any repeat-only trigger separately.

Inject 100–600 K Planck spectra at L_excess/L_star of 0, 0.001, 0.01, 0.03 and 0.1 into measured spectra within 3σ of the peer locus with a Gaia temperature. Normalize the stellar model to its measured H continuum, add wavelength-dependent flux, then rerun band medians and detection. This tests fixed calibrated spectral fluxes and omits source finding, blending and temporal splits.

This specification records the choices frozen before inspecting cohort photometry. Source-level confusion, dust exclusion and all-sky completeness require separate evidence.
