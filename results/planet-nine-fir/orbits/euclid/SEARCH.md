> **Timing provenance:** All orbit-dependent searches, figures, coverage counts and recovery fractions below are historical calculations conditional on the unverified AKARI epoch 2006-12-26 21:03:07.5 UTC. The published Figure 6 shows different December scan labels; the catalog centroid epoch is unverified. These images were already inspected. Scenario propagation does not validate the historical masks or recovery fractions. See [timing provenance](../TIMING.md).

## Euclid VIS moving-source search

No source passed the motion, quality and independent-counterpart tests. Seven 560.52-second calibrated VIS exposures in pointings 5509 and 5510 cover the target on 2024-12-27 from 06:34 to 08:41 UTC (midpoints). These are exposures of one observing sequence. They constitute a 2.12-hour motion test rather than independent nights. ESA released the DR1 images in November 2026. [Euclid DR1 VIS processing](https://arxiv.org/abs/2611.15890)

Using published central positions and 500–700 au, **3,597 of 3,917** independently drawn bound orbits recovered an artificial moving source of VIS_AB 24.0. This yields 91.83% conditional completeness, with 90.93–92.65% Wilson 95% interval for sampling error. With a Gaussian 18-by-7 arcsec IRAS centroid ellipse at PA 79° and a uniform illustrative 40-arcsec AKARI error disk, the measured sensitivity is **2,976/3,928 (75.76%)**. The AKARI error disk has no measured confidence level. Fixing the IRAS effective time to the July 9–10 group found in SCANPI gives **4,042/4,385 (92.18%)** with central astrometry. Unsampled farther orbits, uncertain FIR identifications and planets with different spectra are outside these conditional statements.

| VIS AB | Nominal family recovered | Fraction | After astrometric stress |
| --- | ---: | ---: | ---: |
| 24.00 | 3,597 / 3,917 | 91.83% | 75.76% |
| 26.00 | 3,499 / 3,917 | 89.33% | 73.80% |
| 26.25 | 3,245 / 3,917 | 82.84% | 68.64% |
| 26.50 | 1,843 / 3,917 | 47.05% | 39.31% |

For orientation, an independent photon-counting convolution of the 1-au CALSPEC solar spectrum through the SVO Euclid VIS full system throughput gives the Sun VIS_AB = −26.986. Under the explicit grey geometric-albedo cases in the parent orbit analysis, the selected planet models span VIS_AB = 18.64–23.65. Coloured or dimmer reflectance spectra require a different band convolution. The faintest injection at which 90% of the nominal full family was recovered lies between 25.75 and 26.00 AB. [CALSPEC](https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/calspec) · [SVO throughput](https://svo2.cab.inta-csic.es/theory/fps/fps.php?ID=Euclid/VIS.vis)

### What went into the measurement

1. Intersect the 720 dated ESA DR1 image products in the acquisition cone against 3,617 bound tracks, using geocentric directions replaced by the actual Euclid barycentric observer. JPL Horizons vectors are sampled daily and interpolated to the exposure midpoint. Exactly 39 product footprints overlap, all on 2024-12-27. At least one footprint contains every track. Products include repeated exposure, stacked, VIS and NIR rows and are not counts of independent detections.
2. Fetch science, RMS and flag rows from the original calibrated VIS quadrant FITS files. Foreign sources are hashed in the range manifest. The bad-bit mask uses bits 0, 1, 2, 3, 4, 5, 7, 8, 9, 12, 15, 17, 21, 22. The STARSIGNAL and OBJECTS flags describe real sources and are retained. It takes four valid bilinear interpolation neighbours to produce an analysis pixel. The native photometric conversion is VIS_AB = MAGZEROP − 2.5 log10(ADU / EXPTIME). Local sky is removed on each quadrant and the header's TPV WCS is checked against the reprojected grid. [VIS data product](https://euclid.esac.esa.int/dr/q1/dpdd/visdpd/dpcards/vis_calibratedquadframe.html) · [VIS flag definitions](https://euclid.esac.esa.int/msp/dpdd/v1.1/visdpd/dpcards/vis_flagmap.html)
3. Rasterize the convex hull of the nominal tracks on a 0.1-arcsecond grid at 07:30 UTC. Fit, at each sky cell, a static brightness for every pixel plus a moving point source. Inverse-variance weighting accounts for the actual RMS and masks. Three two-hour motion templates approximate every nominal trajectory to within 0.023 arcsec at each exposure, with a one-pixel Gaussian sigma in the fitting kernel. Of 5,912,008 cells, 5,652,302 meet the coverage test of at least three contributing exposures spanning at least 40 minutes and local coverage greater than 70%. The 0.1-arcsec rasterization and convex hull remain approximations. Individual inaccessible cells and the outside area are counted as unrecovered.
4. Apply a raw threshold of 8 nominal noise units and retain a local maximum in each 25-by-25-pixel window. Of 58 raw peaks, 26 were within 1.5 arcsec of independent DES DR2 sources. The other 32 lose the detection when a single exposure is removed. The minimum allowable score across the seven leave-one-out fits is 3. None survive both checks. The 26 masked positions remain untested for a coincident planet. The thresholds define an empirical search rule, not a calibrated global false-alarm probability. The product dense_vetting.csv gives every peak and its position, reference-source separation and least robust exposure. The brightest static structures generate large positive and negative residuals when dithered frames differ in PSF or astrometry.
5. Sample 8,192 independent uniform triples in effective IRAS epoch, heliocentric distance and signed normalized endpoint-radius change. Retain bound orbits. Inject off-grid moving sources using the true spacecraft-observed coordinates of each trial, with per-exposure empirical stellar PSFs. The images supply 22 qualifying Gaia-star measurements. Exposures with fewer than three usable PSFs use the pooled profile. The on-image 13-pixel diameter normalization uses the median recorded aperture correction of −0.06561 mag. Shot noise uses the 3.48 e-/ADU gain. Add each injection directly to the registered, calibrated pixel data, before static-scene subtraction, and re-run the full and leave-one-out fits. Source-dependent changes in upstream detector reduction, masking or nonlinear response are outside the injection. Each recovered signal must exceed both 8 and any pre-existing competing local maximum. Testing only the nearest spatial trial gives a conservative allowance for subpixel mismatch. Measurements use one independently noisy injection per draw and brightness. Fractions over all accepted draws include masked, inaccessible, and unsearched locations. The quoted Wilson intervals describe the finite sampling and do not include systematic model uncertainties.
6. The geometric intersection, quality-selected area, conditional population completeness and exact-grid recovery are different outputs. For joint survey work, fixed_grid_trials.csv.gz lists shared grid_id and explicit recovered_VIS* decisions after all image, counterpart, persistence and crowding tests. It is a regular grid, and its equally weighted counts need not match conditional Monte Carlo completeness. The field-level evidence uses offgrid_recovery.csv and its stress variants.

### Reproduction

From this directory, with NumPy, SciPy, Astropy, pyERFA and Matplotlib, run the following commands. Disk usage is approximately 2 GB. Queries use public endpoints without credentials. The CSV query receipts include full query URLs, and the original Euclid metadata response is included in compressed form. Products written to the cache can be regenerated.

```bash
python intersections.py euclid_dated.csv.gz
python vis_pixels.py
python mosaic.py
python dense.py
# Packaged Gaia and DES response files are sufficient; catalog.py refreshes them.
python empirical_psf.py
python dense_vet.py
python offgrid.py
python offgrid.py --scenario centroid_stress
python offgrid.py --scenario july_scans
python offgrid.py --scenario fixed_grid
python solar_vis.py
python known_objects.py
```

The original exposure files are multi-gigabyte products. Each requested public byte range, HTTP response length and SHA-256 is listed in range_manifest.json. Saved cutouts are intermediate products. The fixed FITS layout is checked against each downloaded extension's FITS ID, data type and dimensions before use.

### Matched DES and Euclid planet injections

The two-survey comparison uses the same 400 bound states in each DES experiment. It counts recovery in either survey directly. Each orbit has its own Sun and observer distances at each survey epoch. The solar VIS calibration above and a wavelength-independent geometric albedo place the artificial planets on the native Euclid flux scale. We compare the exact target cell with a fresh empirical point source injection, competing neighboring peaks and every leave-one-exposure-out fit. The DES input draw table and its SHA-256 are identified in paired_summary.json.

| Model (mass, radius, grey albedo) | Nominal DES | Nominal Euclid | Nominal either | Stress DES | Stress Euclid | Stress either |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 5 Earth masses, 2.71 Earth radii, 0.10 | 379 | 362 | 397 / 400 | 371 | 271 | 388 / 400 |
| 5 Earth masses, 2.94 Earth radii, 0.30 | 389 | 362 | 400 / 400 | 379 | 271 | 392 / 400 |
| 10 Earth masses, 3.70 Earth radii, 0.30 | 389 | 362 | 400 / 400 | 379 | 271 | 392 / 400 |

For the darkest compact model, the nominal conditional recovery is 99.25% (Wilson 95% 97.82–99.74%) and the DES centroid-stress recovery is 97.0% (94.83–98.28%). The DES experiments draw from a separate 13 × 32 × 56 proposal grid. Their stress test samples combinations of IRAS ellipse-axis extrema and AKARI cardinal displacements of 40 arcsec. The standalone Euclid stress analysis uses continuous offsets. These experiments answer different questions about historical astrometry. An exactly 400/400 result still has a finite sample (Wilson 95% lower bound 99.05%).

The JPL small-body identification service returned no matching known bodies within the field rectangle at 07:30 UTC on 2024-12-27, using the Euclid observer state and the service's two-pass procedure. This is a known-body check at one epoch. [JPL SB Identification](https://ssd-api.jpl.nasa.gov/doc/sb_ident.html)

Reproduce the matched injections with `python paired.py /path/to/planet_draws.csv.gz`. The paired_trials.csv.gz table gives the nominal VIS_AB and both independent outcomes for every draw and model.
