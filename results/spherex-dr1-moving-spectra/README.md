# SPHEREx DR1 moving-object spectral search

Status: investigation checkpoint, **no discovery claim**. Created 2026-11-16.

## Question and observing unit
Can newly released, single-exposure SPHEREx DR1 observations expose repeatable
3-micron absorption or 4.26-micron gas emission in solar-system minor bodies
that is absent from the object's published record?
The observing unit is an independently identified exposure plus an ephemeris
prediction at its midpoint, not a static all-sky-cube pixel.

## Novelty boundary
Before selection we inspected the default branch's `results/` directory names
and open pull requests. Existing work includes many asteroid dynamical,
comet geometric, and solar-wind studies; no SPHEREx moving-object spectral
campaign was present. This is a new instrument/data-analysis lane, not an
extension of an existing candidate. Scientific novelty remains unestablished.

## Primary data and fixed gates
SPHEREx DR1 was released November 12, 2026. Use the NASA/IPAC IRSA spectral
images (DOI 10.26131/IRSA678), not the static coadd, for moving sources.
Archive: https://irsa.ipac.caltech.edu/data/SPHEREx/docs/overview_dr1.html
Moving-object access: https://irsa.ipac.caltech.edu/Missions/spherex.html

1. Establish public access, image coverage and an immutable input manifest.
2. Require motion-aware position and exposure time, detector/channel wavelength,
   uncertainty, masks, and a local background measurement.
3. Freeze target/feature selection before testing independent epochs.
4. Compare off-track positions, stars crossed by the track, and featureless
   continuum controls. Separate rotational variability from spectral structure.
5. Test aperture/background sensitivity, saturation, PSF and thermal continuum.
6. Search literature by object aliases and features; require an independent
   measurement before escalating any possible scientific discovery.

An unavailable product or missing wavelength coverage is an inconclusive
result, never evidence of absent absorption/emission. This checkpoint records
the question before looking at target spectra. Subsequent work selected 29P as a known-gas positive control and 174P/Echeclus for an independent sensitivity pilot; see SCREEN_STATUS.md and CONTROL_RESULT.md. Neither constitutes a discovery.
