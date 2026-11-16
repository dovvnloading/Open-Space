# First real-image control: 29P spectral excess, not a discovery

2026-11-16 checkpoint. We selected a known active comet as a positive control
before applying this new lane to less-studied targets. Carbon dioxide and
variable activity in 29P are already in the literature (e.g.
https://doi.org/10.3847/PSJ/adce07). No novelty claim is made.

## Measurement

The frozen MOST query has 178 image intersections. We selected its 90 D4–D6
rows for acquisition; the archive clamps the end date to October 16, 2026,
so this is NOT coverage of the November 2026 outburst. First 51 downloaded
cutouts are frozen in `control-diagnostics.json`, with original URLs and hashes.

At the moving ephemeris position, 2.5-pixel-radius aperture sums after a local
8–12 pixel median background subtraction are:

| UTC date/time | Central wavelength (µm) | Sum (MJy/sr summed over pixels) | Formal error |
|---|---:|---:|---:|
| 2025-12-16 05:19:33 | 4.168 | 4.355 | 0.224 |
| 2025-12-16 05:26:00 | 4.276 | 13.061 | 0.297 |
| 2025-12-18 14:23:00 | 4.225 | 6.406 | 0.244 |
| 2025-12-18 14:25:09 | 4.261 | 15.647 | 0.325 |
| 2025-12-18 14:27:18 | 4.297 | 5.663 | 0.245 |

The localized excess repeats in different exposures two days apart. The
second sequence brackets it within about two minutes on each side. This is
consistent with the *expected known* gas-band region, not yet a calibrated
line identification, gas production rate, or a discovery. Channel response
integration, pixel solid angles, aperture corrections and systematics are
not included in these sums or their formal errors. No spectral inference
should join observations days apart as if simultaneous.

## Quality trap preserved rather than hidden

Many moving-source core pixels carry OUTLIER bit 19. The official explanatory
supplement §3.2.10 describes single-image clipping after masking known sources;
thus this is not automatically a cosmic ray, but also not permission to
ignore a flag. We preserve strict nominal-pixel apertures AND explicitly
labelled OUTLIER-inclusive diagnostics. Other bad bits remain excluded.
The example 3.074-µm image also has a nonfunctional pixel inside the aperture
and remains rejected. An independent ePSF fit of that image has reduced
chi-square about 5.1, so a point-source model is inadequate there.

The 4.261-µm sequence's core is SOURCE-only flagged and survives the strict
mask: strict sum 15.665 ± 0.324 versus diagnostic 15.647 ± 0.325. Valid
same-image off-target apertures give formal scores -1.49, -2.06 and -1.85;
a fourth is rejected. Four nearby positions do not establish a survey false
alarm probability, and coma/background correlations remain possible.

## Reproduce and limitations

Requires Python, NumPy and Astropy. Run from repository root:

```
python -m unittest tests/test_spherex_moving_screen.py
python examples/spherex_moving_screen.py results/spherex-dr1-moving-spectra/most-29p.tbl --cache /tmp/spherex-cache --output /tmp/29p-diagnostics.json
```

Three focused tests pass (constant-background subtraction, explicit OUTLIER
gate and edge rejection). This is not a full-repository test result. Downloads
sometimes return HTTP 503; failures are retained. The extraction gives
reproducible screening measurements, not publication-grade comet spectroscopy.

Primary product documentation:
https://irsa.ipac.caltech.edu/data/SPHEREx/docs/SPHEREx_Expsupp.pdf
Data DOI: 10.26131/IRSA678. SPHEREx is a JPL/Caltech project funded by NASA.

## New-target gate

A separate 174P/Echeclus MOST query returns 257 intersections, but the predicted
visual magnitude is around 21 at ~13.8 au in one season; this is a feasibility
candidate, not a promised detection. Before promoting any new target, establish
sensitivity, close-in-time feature/continuum coverage, known work and independent
repeats. A failed sensitivity gate must prevent spending a full campaign on it.
