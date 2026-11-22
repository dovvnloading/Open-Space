# Optical, infrared, ultraviolet and radio measurements

Positions refer to the 2MASX nucleus. WISE uses a 3 arcsec cone. PS1 uses the closest mean object within 2 arcsec and ZTF requires a median centroid within 2 arcsec for each light-curve ID. The host contributes unresolved light. Different filters, zero points, apertures and angular resolutions preclude direct cross-survey magnitude differences.

## Same-survey temporal comparison

Medians before the 2014 slew (MJD < 56658, through 2013) and near eRASS1 (MJD 58484–59214, calendar 2019–2020). Positive differences mean fainter. Dates are MJD. Values pool individual detections, so uneven cadence, correlated calibration errors and non-detections limit interpretation. First and last refer to medians in the first and last occupied 180-day bins, including sparse bins. The seasonal table gives counts and scatter.

| Archive | Band | N | MJD range | First bin | Last bin | Before 2014 (N) | 2019–20 (N) | Difference |
|---|---|---:|---|---:|---:|---:|---:|---:|
| AllWISE | W2 | 41 | 55365.15–55558.30 | 12.121 (27) | 12.272 (14) | 12.152 (41) | — | — |
| NEOWISE | W1 | 10 | 56654.74–59575.37 | 12.630 (1) | 12.588 (1) | 12.630 (1) | 12.697 (2) | +0.066 |
| NEOWISE | W2 | 282 | 56651.12–60303.72 | 12.361 (14) | 12.061 (14) | 12.361 (14) | 12.401 (59) | +0.040 |
| Pan-STARRS1 DR2 | g | 7 | 55649.52–56697.66 | 17.440 (2) | 17.580 (2) | 17.539 (5) | — | — |
| Pan-STARRS1 DR2 | i | 15 | 55645.53–56874.26 | 16.911 (4) | 16.853 (7) | 16.911 (8) | — | — |
| Pan-STARRS1 DR2 | r | 6 | 55649.54–56369.41 | 16.950 (2) | 16.820 (2) | 16.950 (6) | — | — |
| Pan-STARRS1 DR2 | y | 7 | 55267.64–56333.66 | 15.893 (2) | 16.391 (1) | 16.159 (7) | — | — |
| Pan-STARRS1 DR2 | z | 11 | 55586.55–56709.66 | 16.605 (2) | 16.489 (2) | 16.533 (9) | — | — |
| WISE (both releases) | W2 | 323 | 55365.15–60303.72 | 12.121 (27) | 12.061 (14) | 12.207 (55) | 12.401 (59) | +0.194 |
| ZTF | zg | 701 | 58202.31–60857.23 | 17.582 (16) | 17.685 (13) | — | 17.659 (250) | — |
| ZTF | zi | 216 | 58232.32–60455.32 | 16.249 (1) | 16.255 (2) | — | 16.342 (129) | — |
| ZTF | zr | 821 | 58216.24–60857.24 | 16.697 (22) | 16.751 (10) | — | 16.737 (310) | — |

## Archive results

- **AllWISE** (measured). 41 rows within 3 arcsec, retained 0 W1 and 41 W2. allwise_p3as_mep
- **CRTS** (unavailable). HTTPSConnectionPool(host='nesssi.cacr.caltech.edu', port=443): Max retries exceeded with url: /cgi-bin/getcssconedb_release_img.cgi?RA=194.7 / HTTPSConnectionPool(host='nunuku.caltech.edu', port=443): Max retries exceeded with url: /cgi-bin/getcsscon
- **FIRST** (unavailable). VizieR query error
- **GALEX MCAT** (unavailable). Request failed.  If this problem persists, please notify archive[at]stsci.edu.
- **GALEX footprint** (coverage_only). 2 footprint records, 4 distinct MCAT products, 0 matched band measurements.
- **Gaia DR3** (no_epoch_product). 3956880840692699392 at 0.49 arcsec, has_epoch_photometry=false
- **NEOWISE** (measured). 293 rows within 3 arcsec, retained 10 W1 and 282 W2. neowiser_p1bs_psd
- **Pan-STARRS1 DR2** (measured). Nearest mean object 136701947144589613 at 0.48 arcsec; 46/53 PSF measurements retained.
- **Swift UVOT fields** (unavailable). expected a CSV table: <?xml version="1.0" encoding="UTF-8"?> <VOTABLE xmlns="http://www.ivoa.net/xml/VOTable/v1.3" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.iv
- **Swift UVOTSSC** (no_match). 0 catalog rows within 3 arcsec, 0 retained bands. Serendipitous catalog has limited mission dates.
- **VLASS** (no_match). 0 peak fluxes; VizieR J/ApJS/255/30 within 4 arcsec. Empty cone gives no flux limit. VLASS component catalog is epoch 1.
- **VLASS CADC fields** (coverage_only). 5 distinct VLASS observation/time records; footprint intersection with CADC ObsCore.
- **ZTF** (measured). 2089 public DR rows, 1738 unique quality-filtered epochs in band; coordinate query.

The WISE release boundary changes processing. W1/W2 magnitudes use Vega profile fits. The 2010–2011 AllWISE multiepoch catalog and 2013–2024 NEOWISE single-exposure catalog are reported separately. Residual blending, scan angle, host light and calibration changes can mimic small offsets. A shared slow change across bands warrants inspection of seasonal sampling. Absence of a catalog row alone gives no upper limit.
