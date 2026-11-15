# Supplemental band-by-band and sampling analysis

This file is an independent additive test; it does not modify the original PR #24 period searches.

At the **fixed 3.829488-hour frequency** and its harmonics, a simultaneous multi-frequency fit of 246 **exposure-paired ZTF g-band difference measurements** yields **A(f)=0.0484 mag** and **A(2f)=0.0037 mag**. The g-band second harmonic is *not detected*: month-block resampling of its fitted amplitude gives a 95% empirical distribution extending to ~**0.020 mag**, and a 20,000-trial one-sided energy-style month sign-flip test of 2f after removing f/2, f and 3f gives **p≈0.86**. This upper range is a distribution of the **fitted** amplitude, not a survey-calibrated upper bound on astrophysical 2f amplitude.

In contrast, TESS's broad red band measures 2f/f amplitude ratios around **0.225 (sector 80)** and **0.165 (sector 92)** in the fully separate raw-image linear reduction. A red-enhanced orbital/ellipsoidal component is a possible physical interpretation, but it is **not established**: ZTF g has sparse sampling and substantially different, changing blend fractions from the broadband TESS pixels, and the g-band non-detection is compatible with multiple causes. No strong color-dependent waveform classification is made.

The Pan-STARRS DR2 detection table contains, after the original PR quality cuts, **27 g detections on 11 nights, 10 r on five nights, 13 i on six nights, 14 z on nine nights, and five y on three nights**, all before the 2017 eruption. Separate fixed-phase period fits to these filter fragments would be dominated by sampling/systematics, not evidence of a pre-eruption clock. ASAS-SN's **249 V-band rows span 2015–2018**, largely before or during the eruptive flares; there are **no V rows from 2019 onward** in the retrieved master-list source. **ASAS-SN g** does overlap the late decline but is dominated by correlated aperture blending and fails a significant phase test (see the main follow-up README).

The phase-locked red waveform's **first harmonic** is therefore a *new* detection in TESS, while the physically important **multicolor harmonic ratio remains unresolved**. Higher-resolution, simultaneous blue/red fast photometry is a decisive next measurement.

To replay the g-band harmonic significance and pre-eruption filter coverage:
python results/v612-sct-clock/followup/band_waveform.py --check --output /tmp/v612-bands.json

All results are conditional on the target having been identified previously in the ZTF discovery screen.
