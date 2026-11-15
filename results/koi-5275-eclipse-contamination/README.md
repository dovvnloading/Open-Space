# KIC 6367628 contaminates the photometry of KOI-5275.01

Kepler observed this known eclipsing binary and the candidate target KIC 6367663 simultaneously in quarters 14 and 15. The stars are 29.97″ apart. Gaia DR3 source 2077647484070832512, at the position of KIC 6367628, has a fitted binary period of 3.77990376 days. KOI-5275.01's reported period is 3.781123 ± 0.001038 days, and its catalog transit depth is 90.4 ppm. The EB's Kepler detection is also cataloged as KOI-6696.01 at approximately half the binary period.

| Quarter | Matched cadences | Relative target response to unit fractional EB flux change | 95% orbit-bootstrap interval | Highest-significance stamp pixel | Its transfer amplitude (electrons/s per unit fractional EB flux change) |
|---|---:|---:|---:|---|---:|
| 14 | 3,918 | 251.6 ppm | 108.4–366.8 ppm | (0, 3) | 3.57 ± 0.20 |
| 15 | 3,971 | 239.3 ppm | 156.6–328.3 ppm | (0, 3) | 4.25 ± 0.23 |

A 40% drop in the binary's flux would produce about 101 ppm and 96 ppm in the two target light curves at the measured transfer coefficients. This is comparable to the candidate's 90.4 ppm depth. These observations establish a contaminating component large enough to explain the shallow event. A separate signal on KIC 6367663 could still contribute. The source lies beyond the downloaded pixel stamp, so these measurements do not yield a PRF position at the binary.

## Method and checks

The script joins quarter-specific calibrated PDCSAP light curves by exact long-cadence number. It normalizes each to its own median and fits the target flux against the simultaneously measured EB flux, using one intercept per binary orbit. No relative time shift is fit. The response has the same sign on odd and even orbits in both quarters. The percentile intervals resample complete orbits 1,000 times.

For the pixel check, the script joins the EB light curve to the target's calibrated pixel file and removes an independent quadratic in time per orbit from the EB predictor and each pixel. The largest transfer-to-bootstrap-error ratio in both quarters occurs at corner (0, 3), toward the off-stamp variable. The KOI is near (2.76, 1.64) and (2.38, 1.10). The EB projects to (−0.83, 8.20) and (−1.22, 7.66). The aperture-independent pixel pattern and source light curve both repeat after the quarterly spacecraft roll.

The original DR25 automated KOI table classified KOI-5275.01 as a false positive. The DR25 supplemental table, vetted in August 2018, and the cumulative table list it as a candidate. A November 30, 2016 DR25 positional-probability report assigned 0.9959 to KIC 6367663 and 0.004137 to generic background. Its table contains no separate entry for KIC 6367628. Here the simultaneous EB flux supplies information absent from that positional model. The published positional probability was conditional on its catalog of possible source positions. [NASA positional-probability report archive](https://exoplanetarchive.ipac.caltech.edu/data/KeplerData/docs/appPdfReports_20170420.tar), member `kplr006367663_K05275_01_dr25_app.pdf`. [NASA cumulative KOI table documentation](https://exoplanetarchive.ipac.caltech.edu/docs/API_kepcandidate_columns.html). [Gaia DR3 eclipsing-binary catalog](https://doi.org/10.1051/0004-6361/202245330).

## Reproduction

Run `python koi5275.py` with Python 3.10+, NumPy, pandas and Astropy. The script downloads six SHA-256-verified MAST products from `input_manifest.csv` and regenerates both CSVs. The tabulated random uncertainty is conditional on the observed Kepler orbits and pipeline photometry.

A [cross-quarter depth test](CROSS_QUARTER.md) trained the dilution in one quarter and predicted the other. The binary predicts 77.0 and 77.7 ppm, compared with measured target dips of 77.6 and 79.0 ppm. The uncertainty on each residual is about 45–50 ppm.
