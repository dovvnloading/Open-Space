# HD 21997 and HD 202628

HD 21997 has no F444W threshold crossing in the 2.5–5.5 arcsec search annulus. HD 202628 has two, at 3.91 and 4.13 arcsec. Both disappear when the identical science integrations are reprocessed with a different reference star. The primary reference produces a deep negative feature whose shoulders pass the locally background-subtracted matched filter. Neither peak establishes a point source.

| Field | F444W date | Science integrations | Reference dithers | Search coverage | Positive / negative crossings | 95% recovery at 2.65 / 5.35 arcsec |
|---|---|---:|---:|---:|---:|---|
| HD 21997 | 2024-02-05 | 5 | 9 | 99.67% | 0 / 0 | 1.58 / 0.74 µJy |
| HD 202628 | 2024-06-02 | 9 | 5 | 99.64% | 2 / 8 | 2.48 / 0.93 µJy |

Each field contributes 2,376 single-source injection trials, with six radii, 36 angles and 11 brightnesses. The published parent code uses independent subsets of five F444W standard-star dithers for the detection and injected PSFs, a relative mask-transmission map, small-sample annular thresholds and a negative-peak control. Its explicit injection before reference fitting agrees with the local calculation to better than 2e−13 in amplitude for both targets. The curves in [measurements](measurements/) measure F444W point-source recovery in a single detector visit. The approximate contrast denominator comes from an isophotal blackbody extrapolation from 2MASS Ks (2.159 µm, 666.7 Jy zero point) to 4.433 µm. The input temperatures of 8300 and 5750 K give 0.720 and 1.748 Jy respectively. Stellar SED, variability and disc excess affect the relative contrasts. Absolute flux recovery does not use this denominator. [2MASS calibration](https://www.ipac.caltech.edu/2mass/releases/allsky/doc/sec6_4a.html)

In the alternate-reference experiment, HD 25945 replaces all five HD 190580 dithers. At the two HD 202628 peak pixels, S/N changes from 5.44 to 0.31 and 8.16 to −0.50. Within two detector pixels of either peak there are zero crossings with HD 25945. The absolute residual at each original peak was negative. [Reference test](measurements/HD202628_reference_test.csv)

## Publication screen

The search on 2026-11-17 covered target names, JWST and companion/coronagraph/point-source terms, the GO 2538 and 3989 programs, and the earlier review's object-linked titles. Search coverage depends on indexing and does not establish priority.

- HD 21997 (HIP 16449, HR 1082). Subaru/HiCIAO SEEDS section 4.3.2 reports zero candidates in the January 2013 data. Rameau et al. surveyed dusty systems with VLT/NaCo. The 2026 Spirit of Lyot conference already reported JWST/NIRCam and MIRI disc imaging for HD 21997 with a broad dust distribution. A dust detection or excess cannot establish a point source. Other object-linked studies include the Meshkat et al. debris survey, ISPY, Gemini/GPI polarimetry, the MLOCI methodology, a two-epoch study of β Pic and gas imaging of HD 131835. These entries alone do not identify a companion in this field. [SEEDS](https://doi.org/10.3847/1538-3881/153/3/106) · [Rameau](https://ui.adsabs.harvard.edu/abs/2013A%26A...553A..60R/abstract) · [Bendahan-West conference program](https://conference.ipac.caltech.edu/SpiritofLyot6/page/agenda) · [ALMA](https://arxiv.org/abs/1310.5069)
- HD 202628 (HIP 105184). Ferrer-Chávez et al. 2026, the JWST-related SIMBAD lead 2026ApJ..1010L..19F, mentions HD 202628 only in §7 as an eccentric-ring comparison. Their detected planet orbits RX J0534.0−0221. HST/STIS had already imaged the ring. Faramaz et al. 2019 Table 4 gives ALMA sources S1–S4. S1 was discussed as possible circumplanetary material with a background alternative. Its existence is previously reported and its interpretation was unresolved in that study. ISPY lists this target in Table D.1. An earlier transit-feasibility study considered a background object that the projected disc would cross. [Ferrer-Chávez](https://arxiv.org/abs/2609.20748) · [HST](https://arxiv.org/abs/1206.2078) · [ALMA + HST](https://arxiv.org/abs/1909.04162) · [ISPY](https://ui.adsabs.harvard.edu/abs/2020A%26A...635A.162L/abstract) · [Background transit prediction](https://doi.org/10.1093/mnras/stt2471)

## Public archive and reproduction

The 0.012-degree MAST footprint cones and frozen product-list responses are in [receipts](receipts/). The NIRCam association members explicitly distinguish one deep science visit per field from PSF reference dithers. For HD 21997, a single public F1140C Level 3 mosaic mixes September 2023 and January 2024. The auxiliary search uses separate PSF-subtracted visit cubes before mosaicking. For HD 202628, both F1140C rolls are from June 2023. GO 9828's October 2026 HD 202628 F2550W images have MAST release dates of 2027-10-17 and do not contribute a public epoch to this test.

From the repository root set the variables below to the field directory containing this README, the merged analysis and a writable cache. Use Python with NumPy, SciPy and Astropy.

```sh
field_dir=$(dirname "$(find results/jwst-five-field-followup -name reduce_f444w.py -print -quit)")
merged=results/jwst-nearby-companions
cache_dir=./jwst-two-field-cache
out_dir=./jwst-two-field-measurements
python "$merged/nircam/nircam_fetch.py" --manifest "$field_dir/inputs.csv" --inputs "$cache_dir"
python "$field_dir/reduce_f444w.py" "$merged" "$cache_dir" "$out_dir" HD21997
python "$field_dir/reduce_f444w.py" "$merged" "$cache_dir" "$out_dir" HD202628
python "$field_dir/check_reference.py" "$merged" "$cache_dir" "$out_dir"
python "$field_dir/auxiliary.py" "$merged" "$cache_dir" "$out_dir"
```

The blind F444W tables treat every threshold crossing as a trigger. Motion comparisons at two wavelengths require an independently recovered centroid and relative registration. A missing counterpart at another wavelength cannot determine motion. Atmospheric grid fits to signed flux proxies on subtraction residuals have no physical interpretation.
