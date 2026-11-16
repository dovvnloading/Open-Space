# Extreme radio spectra at 23 MHz

## Investigation plan

Use the LoDeSS DR1 15–30 MHz continuum images (nominal 23 MHz) to test whether a source's integrated low-frequency emission departs from its optically thin spectrum. A single DR1 wideband image supplies one spectral datum. Start with the good fields and bright, compact, isolated sources. Preserve all selection counts and explicitly distinguish directions outside another survey's footprint from upper limits.

1. Record the original source and Gaussian catalogues, image identifiers, field grades, angular resolution, calibration and observation provenance.
2. Match in sky coordinates to LoTSS DR3 and to independently observed surveys. Give priority to 38–74 MHz data and GLEAM where its northern footprint overlaps. Count multiple components inside the lower-resolution beam.
3. Fit the reference spectrum using higher-frequency detections, with a shared flux-scale term for each survey. Treat upper limits through their likelihood or use image measurements with signed flux. Compare with a curved or absorbed spectrum only when the data support the extra parameters.
4. Test the 23 MHz residual against comparable sources in the same field. Check field position, local noise, ionospheric coherence losses and bright-source artefacts. Use image measurements and alternative resolutions to test blends. Compare overlapping pointings and independent observing epochs where available.
5. For each survivor look up aliases, earlier spectra, variability and source classification. Independent evidence at frequencies constraining the claimed feature is necessary for a physical interpretation. Document exclusions and preserve any null or inconclusive result with the sample denominator.

## Overlap check

Checked the full, untruncated default-branch tree (4,029 entries at a62807165dd86b1c405751042dbd0a98be00a777), the project results index, code searches for LoDeSS and decametre, and the four open PRs (#102, #104, #105 and #106). The nearest repository work concerns LoTSS images of WISE J0458 and multi-frequency pulse timing of GLEAM-X J1627. Neither tackles this question.

Cordun et al. already compared the 23 MHz flux scale with extrapolations based on NVSS, LoTSS, VLSSr and 8C. A repeat of the global comparison cannot establish novelty. The goal here is a named or well-defined population-level result with independent low-frequency evidence and a documented test of alternatives.

[Release](https://lofar-surveys.org/lodess.html) · [Cordun et al., especially §§5.4 and 6](https://arxiv.org/html/2610.04698v1)
