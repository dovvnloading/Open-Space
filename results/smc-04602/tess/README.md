# OGLE-SMC-ECL-04602 | TESS follow-up

NASA TESS / STScI MAST TESScut returned seven 11 × 11-pixel full-frame time series at 15.7324167°, −72.8732778°. For the 2025 images the OGLE ephemeris was fixed at HJD 2460001.4679 + 1.7295874 E, with primary phase 0.00084 and secondary 0.61597. Gaussian width parameters were 0.021 and 0.023 cycles. The two 2025 sectors were the latest available.

## Results

Across 17,501 quality-filtered images, both fixed-template difference signals had their largest positive coefficient in the nearest pixel in each 2025 sector. Ten orbital cycles in sector 94 and eleven in sector 95 passed the free-center folding coverage test. The 3 × 3-pixel free fit places primary and secondary minima at 0.00384 and 0.61434 cycles; the difference is 0.61050, with a percentile 95% interval of 0.59898–0.62024 resampling accepted cycles. This supports the OGLE-IV late eclipse spacing. Every earlier TESS interval overlaps the 2025 spacing, and some earlier difference images peak well away from the source. The TESS images alone do not establish a long-term trend.

The source identification remains aperture-limited. The Gaia DR3 positional counterpart is 0.129 arcsec away (G = 18.345); another Gaia source is brighter (G = 18.067) and 2.115 arcsec away. There are 52 Gaia DR3 sources with quoted magnitudes within a 21 arcsec radius. A TESS pixel is about 21 arcsec wide. The exact OGLE period and eclipse phases provide additional constraints.

## Replay

Python 3.12, NumPy, SciPy and Astropy are required. The two compact 2025 files hold 144-bin, per-orbit detrended flux sums in 3 × 3, four- and one-pixel apertures plus native cadence counts. To rerun the combined 3 × 3 fit and bootstrap with these files:

    python -c "import numpy as np; from test_tess import summary; f=[np.load('folds-%s.npz'%s) for s in (94,95)]; print(summary([v for x in f for v in zip(x['flux'],x['count'])],1200))"

For all 7 original FITS stamps (around 124 MB uncompressed), run:

    python acquire.py 1 27 28 67 68 94 95
    OPENBLAS_NUM_THREADS=1 python test_tess.py 94 95 > tess_result.json
    OPENBLAS_NUM_THREADS=1 python controls.py 94 95 > controls.json
    OPENBLAS_NUM_THREADS=1 python early_test.py > historical2.json

The original FITS SHA-256 receipts appear in the JSON files; TESScut may change packaging headers. Source URLs follow the input parameters in acquire.py. The Gaia DR3 query URL is recorded in the TSV header (Gaia Collaboration, CDS/VizieR I/355). All extractions and stellar timing estimates here are new reductions; the original images are NASA/TESS and STScI/MAST products.

## Extraction

The reduction adapts Open-Space's published 06090 TESS method. It removes detector trends from quality-zero finite pixel time series using 6 off-target common modes plus a block-specific intercept and slope every 0.65 days. The phase calculation corrects TIME (BTJD) by an approximate 70 seconds for HJD_UTC. A phase-binned sinusoidal baseline plus two Gaussians gives each aperture's independent primary and secondary phase and depth, with 0.055-cycle center half-windows. Full-cutout fixed-template fits assess whether the target pixel contains the source of maximum variability. A separate reduction uses a faint-outer-pixel median sky model and a moving median. In 2025 the two secondary centers from that extraction are 0.61643 and 0.61210. The 3 × 3 aperture's undisplaced two-eclipse search statistic ranks first of 144 common phase displacements in both sectors (correlated controls). Absolute and relative depths are aperture- and detrending-dependent.

Original method: https://github.com/dovvnloading/Open-Space/tree/main/results/smc-06090/tess
Archive: https://mast.stsci.edu/tesscut/
Gaia catalogue: https://doi.org/10.26093/cds/vizier.1355
