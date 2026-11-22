# Two-aperture photometry in ObsID 22929

During the registered trial, the equal disjoint cores contained 5 photons at the extraction center and 0 at 2CXO J114044.9-263004. Of the 8 original trial photons, 0 landed in the CSC core and 3 fell outside both cores.

The CSC center is 10.794 arcsec from the extraction center and outside the original aperture. Its separation from the reported X-ray localization is 10.326 arcsec. The 2.780-arcsec X-ray error circle does not intersect the CSC positional ellipse (0.2955632455783342 × 0.2815176193991466 arcsec). The CSC ellipse misses the original extraction aperture in the local tangent plane. The named observation target was SPIDERWEB. The CSC same-ObsID stack match type and broad-band variability index are u / unreported. Code u denotes a unique CSC stack association.

The core radius is 9.873 event SKY pixels. Each core uses the same local background annulus around the midpoint. Scaling divides nominal sky areas. Edge or dither effects can change the effective areas.

| Phase | CCD live s | Original aperture | Trial / CSC core | Shared annulus | Expected background / core | Trial share (exact 95%) | Common supported s |
|---|---:|---:|---:|---:|---:|---:|---:|
| before | 40909.25 | 11 | 5 / 10 | 95 | 5.78 | 5/15 (0.12–0.62) | 2368.12 |
| trigger | 1776.48 | 8 | 5 / 0 | 4 | 0.24 | 5/5 (0.48–1.00) | 0.00 |
| after | 8544.95 | 2 | 1 / 0 | 23 | 1.40 | 1/1 (0.03–1.00) | 235.72 |

The two-sided Fisher comparison of the trial with the pooled flanks gives p = 0.03509.
The table conditions on photons in the two cores. Both centers were selected after inspection, and the flare center was determined using the discovery data. Sparse counts, diffuse emission, PSF wings, astrometric offsets and incomplete geometry limit source attribution. The geometry column gives times supported by independent detector witnesses around the entire construction. The analysis uses no PSF model or exposure map. Its count test has no discovery false-alarm interpretation.

Blue marks the blind extraction core, orange the CSC core, the dashed line the original aperture, and dotted lines the background annulus. Points are individual selected event rows.

![Three phases in matched SKY coordinates](positions.png)

[Counts, associations, detector support and input receipts](measurement.json) · [Tagged EVT2 rows and both offsets](photons.json.gz) · [CSC position and error ellipse definitions](https://cxc.cfa.harvard.edu/csc2.2/columns/positions.html)
