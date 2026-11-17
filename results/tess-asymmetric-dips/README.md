# No validated deep aperiodic occultation in a 384-star TESS Sector 105 search

The strongest apparent isolated stellar eclipse recurs. TIC 316937197 showed a **17.2% fitted QLP depth**, but a separate 3 × 3 pixel extraction gave **0.71%** at the same time and template. Archival TESS visits show eclipses about **18.76 days** apart. The blended source of the eclipse remains unresolved. No event survived as a validated new member of the requested class.

These results cover one June–July 2026 sector, not the population of all stars observed in Cycles 6–8. Sector 105 was the newest verifiable QLP target roster among 107, 106 and 105. We selected 384 distinct TICs, ranks 1–384 by SHA-256 of each unpadded decimal TIC ID in the 399,217-star roster. There was no magnitude cut. Ranks 1–128 were a raw-photometry pilot. Ranks 129–384 formed a separate confirmation sample whose detector was frozen before opening those light curves. [Selection](SELECTION.md) · [Pilot protocol](PROTOCOL.md) · [Confirmation protocol](CONFIRMATION.md) · [Availability receipts](availability.json) · [Mission sector dates](https://tess.mit.edu/observations/sector-105/)

| Sample | Stars measured | Wall-clock star-days | Screen triggers | Stars triggered | Individual injections |
| --- | ---: | ---: | ---: | ---: | ---: |
| Raw-aperture pilot | 128 | 3,554.1 | 34 | 14 | 9,216 |
| Corrected-flux confirmation | 256 | 7,108.1 | 65 | 39 | 18,432 |

## How many events could the search detect?

In the confirmation sample a 2%-deep, 12-hour triangular event was recovered **1,149 times in 2,048 trials (56.1%)**. Injections sampled the entire observing span, including gaps. Effective exposure was **3,987.9 star-days**. Counting all 65 triggers as possible cases gives a **pointwise one-sided 95% upper limit of 0.02195 events per star-day** (8.02 per star-year). This is a conservative ceiling for the specified apparent signal in published DET_FLUX. It allows common rates and offers little constraint on a rare Boyajian-like population.

| Apparent peak depth | 6 h recovered | 12 h recovered | 24 h recovered | 12 h upper limit per star-day |
| --- | ---: | ---: | ---: | ---: |
| 1% | 16.9% | 22.8% | 22.6% | 0.05901 |
| 2% | 55.6% | 56.1% | 56.2% | 0.02195 |
| 5% | 58.8% | 59.6% | 60.6% | 0.02059 |

Each confirmation cell has 2,048 trials. Ingress occupies 20% or 80% of the event with equal probability. Event-plus-shoulder windows must be fully covered on a 30-minute grid. Single-event injections act on delivered flux, so they omit any loss before QLP produced that flux. We used a lower confidence bound for the effective exposure and an upper Poisson count bound. Both bounds are one-sided at 97.5%, with joint coverage at least 95% under independent, homogeneous rare event arrivals for each template population. The count includes recognized binaries and artifacts, so the calculation does not need an assumed efficiency for their human adjudication. No zero-count upper bound was inferred after vetting. The 1% cells sit directly on the detector amplitude threshold. Multi-day, irregular-depth, flat-bottomed, heavily diluted or pre-detrending dips require different calibrations. [Pilot grid](report.json) · [Confirmation grid](confirmation/report.json)

## What produced the strongest signals?

- **TIC 316937197, CD−54 9621.** The Sector 105 peak near BTJD 4217.63 is real in the pixels but much shallower there. The 3 × 3 fit gives 0.71%, with a positive-deficit centroid 0.27 pixel from the TIC position. These are a fixed-template amplitude and a descriptive position. A shallower eclipse also appears in earlier TESS visits, with a clock near 18.76 days spanning multiple years. Gaia resolves another star about 1.56 arcsec from the TIC position, well inside one TESS pixel. Catalog identity and the photometric period do not identify which star eclipses. ASAS-SN recorded no single 5% dip in 2,694 selected archival measurements; its coverage ended about a year before the Sector 105 event. AllWISE W1−W2 is −0.039 ± 0.031 mag, and the W4 entry is an upper limit. The absence of a conspicuous warm excess does not constrain cold debris.
- **TIC 116658378** provides a positive control. Six detections include deep primary eclipses and shallower secondaries. Four primary epochs have a spacing near 7.71 days. Independent pixel fits recover 7–13% for the four primaries, with their deficit centroids near the stellar position. SIMBAD and VSX identify an eclipsing binary there. The measured clock differs from the retrieved VSX period, which we keep as a catalog discrepancy.
- **TIC 2027133756, π² Gru B.** The largest confirmation template in corrected flux gave 18.4%. The pixel fit gives 0.18%, and 13 of 84 off-event control windows were as deep. A Gaia source about 5 arcsec away is over five G magnitudes brighter. The corrected QLP light curve alone cannot establish an intrinsic 18% eclipse.
- **TIC 257556370.** A positive pulse in the raw pilot produced two apparent dips. Capping the 46 high cadences in the pulse removed both triggers. The QLP and small-aperture products show a brightening near BTJD 4209.30. The cataloged rotational variability of this star provides context; the pulse's physical origin was not needed to reject the putative dips.
- **TIC 441246218.** A 2.2% confirmation fit at BTJD 4208.34 had a 2.5% pixel-aperture counterpart. The small and large QLP fits differed (0.9% and 6.5%), and nearby stars show structure at the same time. The pixel result alone cannot assign the deficit to the target. The event remains included in the rate numerator.

## Completed measurements and remaining limits

Quality-zero data only entered either detector. We fitted local linear baselines, measured the response in different QLP apertures, checked recurrence, examined the background and detector peers, and retrieved separate 9 × 9 full-frame cutouts for each of the 53 triggered stars. Each cutout supplied a target-centered 3 × 3 extraction, a per-cadence perimeter sky estimate, a difference image, a descriptive centroid and off-event fits. The independent extraction uses the same telescope and photons. The six nearest QLP targets outside 2 arcmin were requested for each triggered target. The sector data release describes momentum dumps at downlinks and scattered light in the camera fields; quality cuts alone did not remove every background-related signal. [Sector 105 data release](https://archive.stsci.edu/missions/tess/doc/tess_drn/tess_sector_105_drn143_v01.pdf)

All 53 positions returned zero public ZTF light-curve measurements in the queried sector window. ASAS-SN supplied 45 historical light curves, all ending before the TESS events. Other queries returned service errors, unmatched positions or an ambiguous pair. Gaia science photometry ended before this sector; the anonymous ATLAS forced-photometry queue required authentication. These results leave the ground-based epochs unconfirmed. The follow-up also queried Gaia DR3 neighbors, variability and NSS flags, VSX, SIMBAD and AllWISE. Catalog nonmatches provide limited evidence about prior identification. High-resolution eclipse localization and contemporary independent photometry are still missing for the blended cases. Numerical rates retain every trigger.

## Reproduce the two searches

Use Python 3.12, NumPy, SciPy and Astropy. Run from this directory. The fetcher verifies the roster, ranked membership and original FITS byte hashes, then writes a manifest. The pilot originally used lossless NumPy projections of those FITS columns. The replay reads original FITS; container checksums differ while the measured columns are identical. The confirmation runner streams the original FITS and emits their hashes.

```sh
python acquire.py work/pilot --sample sample.json --products products.csv
python scan.py --manifest work/pilot/scan_manifest.csv --out work/pilot-scan
python confirm.py --pilot-sample sample.json --out work/confirmation
```

All raw screening hits and injection outcomes are retained in the per-sample CSVs, compressed injection tables and JSON reports. The published QLP inputs are provided by the MIT TESS Science Office via MAST. TESS is funded by NASA. [QLP archive and data citation](https://archive.stsci.edu/hlsp/qlp) · [TESSCut](https://mast.stsci.edu/tesscut/)
