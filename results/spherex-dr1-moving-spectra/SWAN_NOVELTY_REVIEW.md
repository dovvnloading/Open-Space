# Independent novelty and validity review: C/2025 R2 (SWAN)

Search date 2026-11-16 UTC. Queries combined exact designation, SWAN25B alias, SPHEREx, spectroscopy, infrared, CO2 and carbon dioxide. No directly indexed prior SWAN SPHEREx/CO2 spectrum was found. This is a bounded negative search, NOT proof of novelty, priority or absence of unpublished/in-press work.

## Verified prior observations

1. University of Jena primary observing report: https://www2.astro.uni-jena.de/Users/markus/C2025R2.html
   Optical spectroscopy on 2025-10-18, C2, forbidden oxygen and NH2 features; no detected sodium. The observed optical wavelength range does not test 4.26-micron CO2 emission.
2. TRAPPIST ATel 17538, published 2025-12-10: https://orbi.uliege.be/bitstream/2268/338585/1/ATel%2317538.pdf
   SWAN measurements 2025-11-27 at r_h=1.61au, Delta=.90au: Q(OH)=7.52 +/-2.92e26/s, Q(CN)=1.19+/-.07e25/s, Q(C2)=1.41+/-.10e25/s, Q(C3)=3.09+/-.64e24/s. These are optical narrowband Haser-model estimates, not direct infrared parent CO2 measurements. They are 25 days before December22 SPHEREx pilot and cannot be divided directly into a December22 CO2 estimate without activity-evolution modeling.
3. Levine et al., Assessing Potential Visibility of Comets Like C/2025 R2 by NEO Surveyor: https://ntrs.nasa.gov/citations/20250009793
   Visibility simulations, not an observed SWAN spectrum. Does not establish chemistry.

## Wierzchos is a positive control, not a new CO2 discovery

First JWST spectrum of distant activity in long-period comet C/2024 E1 (Wierzchos), DOI https://doi.org/10.1093/mnrasl/slaf046 (https://arxiv.org/abs/2503.14071), already reports CO2-driven activity at7au, water ice, and absent detectable water/CO gas. A July2025 SPHEREx 4.26-micron feature is therefore expected known chemistry unless a separately demonstrated new temporal/spatial/compositional result is obtained.

## Required validation of SWAN pilot

Engineering reports December22 feature in adjacent samples, but values are preliminary aperture diagnostics. Do not describe these as validated gas emission, unusually CO2-rich chemistry or discovery yet.
- Recover missing blue-flank image and independent repeated sequences.
- Integrate continuum+CO2 emission model through calibrated local spectral response, not simply central wavelengths. A peak at4.293 vs4.256 micron needs forward modeling.
- Use actual target WCS and ephemeris at each exposure; inspect spatial source and contamination visually.
- Retain OUTLIER sensitivity explicitly: moving uncatalogued sources may be flagged by per-image clipping, but unrelated artifacts may too.
- Compare larger/smaller apertures and annuli, ePSF/extended-coma models, detector-row controls and co-moving repeated detections.
- Propagate background, aperture, relative calibration and temporal variability errors, not just independent pixel variance.
- Do not infer a CO2/H2O abundance ratio from noncontemporaneous optical OH without physical modeling.

A confirmed, previously unreported quantitative spectrum could add real observational information, but ordinary CO2 detection alone is not evidence for extraordinary composition or a high-impact discovery. Stronger scientific question: whether post-perihelion gas-to-dust or CO2/H2O evolution differs robustly from physically matched observations; pre-register and test this only if data span permits.
