## Sector 109 transit search

The three white dwarfs with usable apertures had zero periodic or isolated triggers in Sector 109. Their 33,519 quality-zero exposures span 2026 October 5–31. The brightest star rejected by the crossmatch, 40 Eri B, is part of a close multiple. Five of the remaining eight failed the dilution cut. The three retained TIC IDs are 237397195, 43529091 and 55997020.

| Sample step | Count |
|---|---:|
| All-sky Gentile Fusillo EDR3 main catalogue with Pwd > 0.75 and Gaia G < 14 | 157 |
| On the predicted Sector 109 silicon, each with one verified 9 × 9 pixel cutout | 9 |
| Unique unflagged TIC 8.2 entry with matching Gaia identifier and <2″ J2000 separation | 8 |
| At least 80% model target light across the adopted PRF and colour range | 3 |
| Quality-zero exposures per accepted target | 11,173 |

The roll and detector footprint are from tess-point 0.9.5.post1. Gaia proper motions propagate source coordinates from 2016.0 to 2026.8 for these two sectors and separately to J2000 for the TIC match. The 3′ Gaia DR3 cone is centered at the catalogue position. The detector WCS transform uses the earlier asymmetric-dip search's `pixels.target_pixel`. A 1.5-pixel-radius aperture is fixed on the propagated WD position. Each cadence has a local sky estimate from the 9 × 9 perimeter median. The code also computes the positive-flux centroid within three pixels, the aperture error, and the pipeline quality flag. The search uses the total measured aperture flux. Its injection model multiplies theoretical flux losses by the target fraction, which avoids amplifying detector noise. BTJD values from TESScut use the FFI timing correction. Cross-cadence comparisons remove the respective TIMECORR values before comparing spacecraft times.

The nominal fraction of target light uses the TIC target magnitude, Stassun et al. (2019) equation 1 for Gaia neighbours and a pixel-integrated Gaussian with sigma 0.85 pixels. The crowding cut requires at least 80% in all models with sigma 0.65, 0.85 and 1.05 pixels. The adverse-colour model makes the target 0.3 mag fainter and assigns BP−RP = 3.5 to neighbours with missing colours, compared to 1.0 nominally. All neighbour colours are bounded to the calibration range [−0.2, 3.5]. These are approximate PRFs and broad-band colour conversions. Unresolved companions, extreme colours, scattered light and inter-pixel variations can alter the actual dilution. The TIC match deliberately excludes ambiguous or renamed Gaia IDs. The record includes each target's largest predicted contaminant contributions.

For physical injections the pure-H mass and gravity columns give the radius via R = sqrt(GM/g). Stellar model errors matter especially for TIC 43529091, the previously published 9.5-hour ELM WD binary. The H-atmosphere catalogue values are the assumptions used for injection recovery, not spectroscopic measurements of this binary. Kawka et al. report orbital beaming without eclipses. [Kawka et al. 2020](https://doi.org/10.1093/mnrasl/slaa068)

### Auxiliary checks on the same three stars

TESScut returned five earlier 200-s visits from Cycle 5 onward through Sector 108, all passing the aperture cut. TIC 237397195 had Sectors 87, 98 and 108. TIC 43529091 and TIC 55997020 each had Sector 98. The searches used 93,884 additional quality-zero exposures with zero periodic triggers. Sector 98 yielded 24 individual box triggers clustered within BTJD 3991.087–3991.710 across all three stars. The strongest individual thresholds in S87, S108 and S109 were below the detection cut. The image records give central-pixel, three-by-three-pixel and science-aperture depths, background, centroid changes, and nine-by-nine difference images for each trigger. Each trigger's 2.4-hour rolling stamp-edge noise was 1.68–2.39 times the median for that visit.

All 572 contributing 20-s in-window exposures in the 24 trigger boxes carried SPOC bit 4096, “Scattered Light Exclude”. The 20-s record lists counts per event and the independent target-pixel file hashes. This target-dependent quality flag is unavailable in the FFI quality word. Multiple targets and detectors share the affected time interval. These flagged windows cannot support a debris or planetary claim. We do not assign a unique physical cause to the excess noise. No same-period ground fold is defined for these isolated, excluded events. The auxiliary observations have no occurrence denominator. [MAST quality flags](https://outerspace.stsci.edu/spaces/TESS/pages/14563420/2.0+-+Data+Product+Overview)

Nine Sector 110 geometric targets were also requested via TESScut. At the recorded retrieval time each returned the JSON status “No FFIs found for input target. Cutout not performed.” Sector 110 was still underway. Individual calibrated FFI availability is documented separately in the main results.

### Replay the FFI measurements

Use Python 3.12 with numpy, scipy, astropy and tess-point==0.9.5.post1. Run from the package containing these scripts and `wd_search.py`. Supply the path to the existing `results/tess-asymmetric-dips/pixels.py` if this package sits below a nested scripts directory.

```bash
python ffi_acquire.py cache --sectors 109 110 --workers 4
python ffi_photometry.py cache lc109 --roster 's109.csv' --pixel-tools ../tess-asymmetric-dips/pixels.py
python ffi_images.py cache lc109 wd_search.py s109_images.json
python ffi_history.py cache lc109/ffi_visits.csv
python ffi_photometry.py cache historic --roster history.csv --pixel-tools ../tess-asymmetric-dips/pixels.py
python ffi_images.py cache historic wd_search.py historical_images.json
python ffi_quality.py cache historic historical_images.json historical_quality.json
```

For `python ffi_fastcheck.py cache historical_images.json historical_fastcheck.json`, first download the three 20-s TPFs named in the frozen `historical_fastcheck.json` to `cache/` from MAST. Each record gives a filename, archival URI and full file SHA-256.

The receipts record the query URLs, response sizes, times, and SHA-256 values for input FITS and catalogue responses. The catalogue data-body hashes exclude changing query timestamps in response headers. FITS products generated by TESScut may change checksums when service-written metadata changes. Accepted/rejected measurements with a TIC match contain separate TIME, FLUX and QUALITY array hashes. The frozen all-sky source selection and per-sector rosters are included. Re-running S110 queries as the archive ingests will legitimately return new data. Compare the S109 scope to its frozen receipt and keep any new sectors separate.

[Gentile Fusillo et al. EDR3 catalogue](https://doi.org/10.1093/mnras/stab2672) · [TIC v8.2](https://cdsarc.cds.unistra.fr/viz-bin/cat/IV/39) · [Gaia DR3](https://doi.org/10.26093/cds/vizier.1355) · [TESScut](https://mast.stsci.edu/tesscut/) · [TIC colour transform](https://doi.org/10.3847/1538-3881/ab3467)
