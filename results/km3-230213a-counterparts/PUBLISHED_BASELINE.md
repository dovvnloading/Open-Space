# Published baseline: every source is KNOWN

[known_17.csv](known_17.csv) transcribes the 17 source names, positions, rounded separations and three published redshifts from Lincetto et al., PoS(ICRC2025)1100 Table 1 (https://pos.sissa.it/501/1100/pdf). These sources were already in the Nature 2025 (https://doi.org/10.1038/s41586-024-08543-1) and subsequent study, Astronomical Journal 172, 139 (September 2026; https://doi.org/10.3847/1538-3881/ae7a4e). The #13–#17 entries appear at the end of the proceedings table; numbering is not globally sorted by separation. #2 is a gamma-only source whose nature may be confused with mismodelled Galactic diffuse emission.

| Literature report | Assessment |
| --- | --- |
| #8 PMN J0606−0724 (z=1.277) prominent OVRO 15 GHz flare, ~5 d from neutrino; 0.26% chance in ~3000 OVRO light curves | **KNOWN**, p=0.0026 PRE-TRIAL from published work only; NOT a global p-value; collaboration rejected conclusive association |
| #6 0605−085 (z=0.87) year-long LAT gamma flare with peak ~180 d before event | **KNOWN**; post-trial association p not established |
| #1 MRC 0614−083 X-ray increase between eROSITA visits and Swift ToO in late 2024 | **KNOWN**; multi-year gap prevents measurement of timing at neutrino epoch; p not established |
| ASKAP+VLASS 1052 radio sources within 1.5 deg, 10 blazar-like, UGCA 127/Phaedra, WISEA J061715.89−075455.4/Hebe, EMU J062248−072246/Narcissus | **KNOWN** from Filipović et al. (https://arxiv.org/abs/2503.09108); no new source attribution |
| Cosmogenic/UHECR origin and steady blazar population | **KNOWN HYPOTHESES**, not a source identification. Cosmogenic https://arxiv.org/abs/2502.08508 ; blazar population https://doi.org/10.1088/1475-7516/2026/03/033 |
| Nature flux-versus-IceCube/Auger comparison | **KNOWN**; model-dependent joint configuration ~0.5% for assumed E^-2 diffuse flux, NOT an association p-value |
| IceCube 2026 15-year steady/flaring/time-window study | **KNOWN NULL**; no excess, flux limits: https://arxiv.org/abs/2609.10931 |

**Already tested, not unexplored:** collaboration also examined OVRO, RATAN-600, VLBA, Fermi-LAT, WISE/NEOWISE, ZTF, ATLAS, eROSITA/Swift, SVOM pointings, and a VLBA VLASS-selected follow-up. Reprocessing these without newly released input or a frozen alternative selection is not independent. A literature-only match to 2026 studies is not novel.

**Statistical guardrail:** The 17 records are preselected for angular proximity: sky shifting *these* 17 would circularly give inflated significance. The chance probability for *each* one is **not estimable from this selected table**, except the externally published 0.0026 unadjusted temporal p of #8. Reliable chance estimates require full parent-catalogue density, exposure, coverage and independent controls. Uniform whole-sky Poisson density is an illustrative prior only (4FGL-DR4 ~7195/41253 sq.deg and 5BZCAT 3561/41253 sq.deg); these do NOT represent Galactic-plane-calibrated source probabilities. The duplicate catalogues may refer to the same physical objects.

**Receipt:** downloaded/read 2026-11-18: https://pos.sissa.it/501/1100/pdf ; https://pmc.ncbi.nlm.nih.gov/articles/PMC11821517/ ; https://doi.org/10.3847/1538-3881/ae7a4e . Proceedings PDF text was parsed; screenshot attempt returned HTTP 403, so manually crosscheck when endpoint is available. Coordinates are J2000 degrees. The analysis script checks their spherical separations.
