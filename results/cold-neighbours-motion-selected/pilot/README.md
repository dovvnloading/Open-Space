# Blind CatWISE search near the Euclid Q1 fields

Two previously published cold dwarfs reproduce their expected motion in the independent unTimely W2 extraction. A third bright catalogue lead has zero eligible unTimely centroids within 2″ of its predicted position across the sampled epochs. These three exhaust the specified bright tier in 84.8036245 deg² of CatWISE2020. Search depth here refers to catalogue thresholds. No effective volume or survey completeness has been measured.

| Cone (ICRS center, degrees) | All CatWISE rows | W2 S/N ≥8 | Loose colour and motion | Primary selection | Bright tier |
|---|---:|---:|---:|---:|---:|
| North (269.733, +66.018) | 1,589,282 | 1,083,093 | 8,665 | 617 | 2 |
| South (61.241, −48.423) | 1,122,018 | 487,950 | 1,509 | 24 | 1 |
| Fornax (52.932, −28.088) | 890,260 | 334,388 | 1,125 | 17 | 0 |
| **Total** | **3,601,560** | **1,905,431** | **11,299** | **658** | **3** |

Each cone has radius 3°. The sum of their disjoint spherical areas is 3 × 2π(1 − cos 3°) sr. Euclid Q1 imaging covers irregular subsets of these circles. Instantaneous sky position for zone selection is CatWISE's catalogue position. Individual trajectories use the motion-fit position at MJD 57170 (2015.405).

## Selection and limits

The first pass requires motion-fit W2 S/N ≥8 and raw total motion ≥0.20″/yr. It retains W1 S/N <3, missing W1 S/N, or measured motion-fit W1−W2 ≥1.5 Vega mag. Three separate count queries check the source, S/N and prefilter totals in each cone. Primary selection requires colour ≥2.0, tile-corrected motion ≥0.30″/yr and at least one motion component exceeding five times its catalogue uncertainty in quadrature with 0.035″/yr. The latter is a screening threshold, without a calibrated false-alarm probability. The two bands share a magnitude zero point of 22.5 for the motion-fit fluxes. For W1 S/N <3, the retained colour bound uses 2.5 log₁₀[(F₂−1.645σ₂)/(max(0,F₁)+1.645σ₁)]. The separate one-sided bounds on the fluxes have nominal 95% Gaussian coverage individually. For measured W1 (S/N ≥3), the cut uses the point colour. The published tile corrections are added before applying the final motion and position cuts.

The bright tier additionally requires W2 ≤16.5 Vega, W2 motion-fit reduced χ² ≤3, W2 unWISE artifact mask zero, AllWISE W2 contamination character zero or absent, and both raw motion components below 99.99″/yr in absolute value. A blank AllWISE flag supplies no positive quality information. The full primary table retains one row with its motion pegged at the fit boundary. Fit quality, confusion, W1 limits, high motion and low tangential velocity can all reduce sensitivity. At 5, 10 and 20 pc the bright cut alone corresponds to absolute W2 ≤18.01, 16.50 and 14.99 respectively (ignoring extinction). These conversions do not establish a temperature or distance limit.

## Epochal cross-check

The unTimely catalogues provide W2 detections extracted independently from individual unWISE epoch coadds spanning 2010–2020. CatWISE and unTimely share WISE exposures through 2018. Every file used has a SHA-256 digest. For each tile the script chooses the catalogue with the most detections per 182.625-day interval, anchored at MJD 55197. At the ecliptic poles a chosen catalogue can have partial or zero local coverage. Within an epoch each tested track receives credit for at most one centroid with flux/error ≥4 and quality factor >0.85. The summary compares a 1.5″ window on the CatWISE trajectory to a stationary trajectory and three trajectories rotated about the catalogue reference position. This association diagnostic carries no astrometric distance information. A failed match can reflect blending, an artifact, low signal or imperfect local coverage.

| Bright source | W2 | W1−W2 | Total CatWISE motion (″/yr) | Epochs within 1.5″ moving / fixed | Interpretation |
|---|---:|---:|---:|---:|---|
| J181125.34+665806.4 | 16.024 | 3.037 | 0.411 | 16 / 9 | Previously published object |
| J041358.20−475038.1 | 15.600 | 3.555 | 0.321 | 17 / 11 | Previously published T9 |
| J182532.09+671906.2 | 16.477 | 2.834 | 2.686 | 0 / 0 | Catalogue motion lacks epochal support, closest eligible centroid 4.925″ |

## Replay

Run `python fetch_catwise.py`, then `python select_catwise.py` from this directory for the standard-library catalogue stages. Each loose-result CSV is validated against an independent count and source-ID uniqueness. Raw TAP responses and their query texts are retained with hashes in `catwise_manifest.json`. The offset input is the [published CatWISE tile correction table](https://irsa.ipac.caltech.edu/data/WISE/CatWISE/gator_docs/CatWISE2020_Table1_20201012.tbl), SHA-256 `c3a03ee558ee31ac81c4b905615f7f2f027094fea9df3d49609def7f33931247`.

For the trajectory audit install NumPy and Astropy, then run `python check_untimely.py` and `python score_tracks.py`. Set `UNTIMELY_CACHE` to a download cache with several GB of free space. Setting `UNTIMELY_EXTRA_TILES` to an empty string limits the replay to the two bright-tier tiles. The epoch FITS inputs and archive index are recoverable from the exact URLs and hashes in `untimely/manifest.json`. The extraction preserves nearby competing detections. Consult the main result for external catalogue matches, spectra, image-level evidence and control recoveries.

[CatWISE2020, Marocco et al.](https://doi.org/10.3847/1538-4365/abd805) · [unTimely, Meisner et al.](https://doi.org/10.3847/1538-3881/aca2ab) · [unTimely data](https://doi.org/10.26131/IRSA580) · [Euclid Q1 field geometry](https://www.euclid-ec.org/science/q1/)

The catalogue of detected neighbours also supplies an empirical per-epoch depth proxy. Within 25″ of the reference coordinate the diagnostic takes the median of `22.5 − 2.5 log10(4 × dflux)` over detections with quality factor >0.85 and at least three contributing frames, then takes a median over epochs. For J182532 this is W2=16.19 and four of twelve usable epochs have proxy limits fainter than the source's W2=16.477. The two known sources are brighter than their proxies in all sixteen and seventeen usable epochs respectively. Noise structure, blank positions and crowding can differ within each aperture, so this estimate cannot measure completeness.
