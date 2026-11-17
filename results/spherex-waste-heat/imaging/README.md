# Where the WISE emission falls

**In Hephaistos B, C and G the 12 µm centroids fall 3.06, 3.28 and 5.54 arcsec from the stars. At 3.4 µm the offsets are 0.13, 0.13 and 0.12 arcsec.** These measurements corroborate wavelength-dependent contamination reported from archival imaging and radio work. The centroids alone cannot apportion the 12 or 22 µm flux between the stars and nearby sources. [Native images](atlas.png) · [Measurements including signed components](centroids.csv) · [Archival diagnostics](https://arxiv.org/abs/2607.03619) · [Radio interferometry at G](https://arxiv.org/abs/2501.05152)

| Field | W1 offset (arcsec) | W3 offset (arcsec) | W4 offset (arcsec) | W3 local placement p90 | W4 local placement p90 |
| --- | ---: | ---: | ---: | ---: | ---: |
| A | 0.07 | 0.84 | 1.10 | 0.20 | 0.51 |
| B | 0.13 | 3.06 | 4.00 | 0.21 | 1.67 |
| C | 0.13 | 3.28 | 3.72 | 0.41 | 1.97 |
| D | 0.06 | 1.07 | 2.42 | 0.40 | 1.86 |
| E | 0.12 | 0.89 | 1.71 | 0.39 | 2.29 |
| F | 0.16 | 0.79 | 2.12 | 0.61 | 3.09 |
| G | 0.12 | 5.54 | 3.13 | 0.68 | 2.15 |
| H | 0.17 | 2.78 | 0.75 | 1.79 | 2.75 |
| I | 0.10 | 2.03 | 4.16 | 6.92 | 3.82 |
| J | 0.10 | 2.16 | 5.32 | 1.24 | 4.13 |

A's W3/W4 displacements exceed the local placement scales, but remain below 1.1 arcsec. F's W4 fit changes by 4.48 arcsec when its adopted PSF width increases 20%, larger than its fitted 2.12 arcsec separation. At I the corresponding W3 width sensitivity is 2.37 arcsec and the local W3 placement p90 is 6.92 arcsec. In D and E, background galaxies about an arcsecond from the stars dominate the published longer-wave photometry even though the WISE offsets are modest. [JWST/MIRI separation](https://arxiv.org/abs/2607.09460)

## Images, controls and limitations

Each wavelength panel spans 32 arcsec. The Gaia position is propagated from J2016 to 2010.5. We fit a circular Gaussian and a sloping background simultaneously on the native Atlas pixels, with adopted FWHM 6.1, 6.5 and 12.0 arcsec in W1, W3 and W4. Twenty-four local placements per band and target test recovery of a Gaussian at the fitted brightness. The table shows the 90th percentile of radial displacement in those placements. PSF shape error, inter-band astrometric registration and blends at the target remain outside this statistic. A small fitted offset can be significant relative to local noise without establishing a source identity. [PSF-width sensitivity, offsets and control scales](centroids.csv) · [AllWISE Atlas](https://irsa.ipac.caltech.edu/ibe/docs/wise/allwise/p3am_cdd/)

The separate 48 arcsec Legacy Survey DR10 and unWISE-NEO6 views cover all ten targets. C is conspicuously crowded in the optical image. Fourteen Gaia entries lie within the ten 15 arcsec cones. Three extra sources are around C and one is around I. The target-only Gaia results elsewhere do not exclude faint, obscured, or unresolved neighbors. SIMBAD's initial combined query timed out. All ten individual cones subsequently succeeded through the separate warm-target follow-up. They recover the known radio associations in A, B and G. [Optical and unWISE contact sheet](optical_unwise.jpg) · [Gaia entries](gaia_neighborhood.csv) · [Completed SIMBAD cones](../warm_fgk/simbad_receipts.json) · [Identifications and offsets](../warm_fgk/simbad_neighbors.csv)

## Replay

Use `python neighborhoods.py` for the viewing products and Gaia queries. To repeat the WISE fits, unpack `../wise/archive_inputs.zip` under `../wise/archive/`, then run `python atlas.py` with NumPy, SciPy, Astropy and Matplotlib. To replay the archived pixels, first extract [atlas_inputs.zip](atlas_inputs.zip) here to restore `stamp_cache/`. The script verifies the SHA-256 of every cached cutout. Eighteen 135-pixel cutouts were returned by the IBE service. For twelve others it read IRSA's publicly mirrored native Atlas files and computed the same-size cutout using the FITS WCS. Each of the 30 target-band receipts identifies its bytes, acquisition route and FITS hash. [Image query and hash receipts](queries.json) · [Optical display receipts](viewer_queries.json) · [IRSA cloud archive](https://irsa.ipac.caltech.edu/cloud_access/)
