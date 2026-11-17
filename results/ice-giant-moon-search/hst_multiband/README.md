# Hubble WFC3/UVIS search of Uranus and Neptune

**No unidentified moon passed independent-image checks in 104 searchable filter epochs.** We retrieved and hashed all 1,427 WFC3/UVIS exposures in the archive inventory. The remaining five Neptune filter epochs had only three visits with a registered planetary disk. The search tested 75,199,973 supported circular, equatorial orbit cells. All ten cells surviving the initial artifact checks came from Uranus FQ619N images. Their frozen trajectories failed follow-up in independent exposures. [Epoch index](epoch_index.csv) · [FITS hashes and MAST URIs](exposure_hashes.csv) · [Per-exposure and per-template tables](epoch_details.tar.gz)

| Planet | Retrieved exposures | Registered | Searched filter epochs | Supported orbit cells | FQ619N cells followed up |
|---|---:|---:|---:|---:|---:|
| Uranus | 563 | 563 | 43 | 18,115,418 | 10 |
| Neptune | 864 | 749 | 61 | 57,084,555 | 0 |

The census contains 109 filter epochs in nine passbands. An epoch groups observations with the same planet, program and filter, splitting at gaps longer than three days. Registration accepted 1,312 exposures and rejected 115 Neptune images whose measured bright region failed the expected planetary disk-area test. Of these registered exposures, 1,297 contributed to searchable epochs and 15 belonged to the five three-visit groups. Each supported orbit cell is one radius/phase template counted in one epoch. Neighboring cells share pixels. [Search and registration code](screen.py) · [Retrieval and selection](../ARCHIVE.md)

## Searches and measured depth

Uranus radii span 74,500–84,000 km and 110,000–155,000 km. Neptune radii span 68,000–190,000 km. Detector masks remove the disk and its immediate halo, the adopted rings and arcs, flagged pixels, sharp events and known moons. Each frame uses a limb registration, a distortion-aware local projection and its own light-time correction. The radial grid spacing is 100 km and the azimuthal spacing is approximately 550 km. Sampling with alternating disjoint visits, exposure dropout and a common orbit tests time coherence. The tested orbits follow each planet's equatorial plane. Orbits with inclination or eccentricity can have different recovery fractions. [Detector masks, projection and score](../hst16995/reduce.py) · [Horizons query identities and hashes](position_queries.csv)

In the 79 epochs with repeated exposures per visit, 19,389 initial cells triggered and 19,246 retained sufficient contrast after the strongest exposure was dropped. Subtracting the median planet-aligned image from *other* visits eliminated all 19,246. The separate sparse procedure admitted epochs with at least four visits and as few as one usable exposure per visit. Each alternating half required a normalized stack of five, at least three visits had to exceed 2.5 individually, and the stack excluding the strongest visit had to exceed five. In the 25 searchable sparse epochs, 112 cells crossed this screen and ten survived the planet-aligned image subtraction. Residual scores use a radial median and median absolute deviation. Correlated image residuals and the template count preclude interpreting these scores as Gaussian significances. [Leave-visit-out subtraction](stationary.py) · [Sparse selection](sparse.py)

Validated depth comes from 1,792 blind image-level plants in the two program 16995 F845M epochs. The measured 50% transitions are total AB 22.39–22.89 in the inner radial interval of each planet and 22.89–23.39 in the outer interval. Illustrative diameter brackets range from 26–71 km at Uranus and 61–166 km at Neptune, depending on radius, bracket edge and albedo (0.04–0.12). Other epochs and filters have different exposure times, halos and masks. Their search coverage is measured. Equivalent image-plant completeness curves remain unmeasured. [Per-planet recovery fractions, uncertainties and size assumptions](../hst16995/README.md)

## The two Uranus FQ619N epochs

**Program 16995, November 2022.** Five adjacent cells at 138,900–139,200 km follow the same three excursions. The highest-ranked cell at 139,000 km has a normalized FQ619N stack of 8.55 across six usable visits. The frozen ephemeris gives −0.92 in 16 independent F845M exposures and −0.65 in eight F763M exposures. None of the five cells reaches 2.5 in any individual follow-up visit once exposures are combined. In the original FQ619N frames, the strongest cell gives twice roughly 500 electrons in a 2.5-pixel aperture (AB 20.59 and 20.63), then −39 ± 27 electrons at another supported visit with identical exposure time. One of the bright excursions has 71% of its local 5 × 5-pixel positive signal in one pixel, about four pixels off the frozen prediction. The nominal known moons are more than 87 pixels from the three strongest contributions. These measurements provide no repeatable point-source photometry on this orbit.

**Program 18456, October 2026.** Five cells span 74,500–119,700 km. The closest to independent corroboration was the 110,800-km cell. Its F845M score was 3.63 in 20 usable exposures across eight visits, including two visits above 2.5. F763M gave −0.79, F467M 0.04, F547M −0.33, F657N 0.09 and FQ727N 0.32. In the eight original FQ619N exposures, a valid 2.5-pixel aperture gives 2,186 ± 50 electrons at visit 08 (AB 18.99). Four other supported exposures with identical integration times give −7 to 32 electrons and empirical 5-MAD limits of AB 22.21–22.68. The apparent source in the highest-scoring visit has its intensity centroid 4.2 pixels from the nominal trajectory, with bad pixels in the aperture. The independent F845M excess and the bright FQ619N measurements do not describe one repeatable moon. Across all six independent filters, the other four cells produce normalized scores between −3.16 and 2.18. [Frozen independent-filter scores for all ten cells](sparse_followup.csv)

For every survivor, the follow-up freezes phase at the mean *emission* time of the original accepted exposures, propagates at its own mean motion and verifies numerical replay of the source score before changing filters. Localizations at distinct visits were also tested against fixed detector pixels and fixed-sky positions. The planet-tracking pointings usually carried the sidereal test position outside the supported field in different visits. Within-visit comparisons in other filters provide the useful fixed-sky tests. The data distinguish very weak repeatable emission from the reported bright excursions only to the sensitivity of those exposures. The FQ619N aperture values are native detector aperture magnitudes, without an infinite-aperture correction. [Forced tracks, morphology and motion tests](sparse_vet.py) · [Native aperture photometry and background controls](sparse_photometry.py)

## Reproduce

From this directory, unpack the diagnostic tables and run the frozen-source replays. The cache accepts a directory of previously retrieved native FLT files. The reducer validates each complete FITS hash. NumPy, SciPy, Astropy and Requests are required.

```sh
tar -xzf epoch_details.tar.gz
python sparse_vet.py uranus-16995-fq619n-59892 --cache /path/to/hst-cache
python sparse_vet.py uranus-18456-fq619n-61332 --cache /path/to/hst-cache --filters f845m f763m f467m f547m f657n fq727n
python sparse_photometry.py /path/to/hst-cache uranus-16995-fq619n-59892 uranus-18456-fq619n-61332
```

Inside the archive, each epoch directory contains the input checksums and individual registration decisions. The two FQ619N survivor directories also contain `sparse_vet_filters.csv` (forced scores and same-radius comparison tracks), `sparse_vet_exposures.csv` (positions, masks, nearest moons and morphology), `sparse_vet_motion.csv` (other-visit anchored tests), `sparse_vet_anchors.csv` (same-visit cross-filter motion tests), and `sparse_vet_apertures.csv` (native electrons, empirical noise and aperture AB). The aperture errors combine the robust dispersion of nearby apertures with source Poisson noise. Templates were selected using the original filter, so apparent tail fractions in that filter are selection-biased.

For a fresh full census, run the following in a copy of this directory without precomputed coverage, vet or sparse summary CSVs. Those files are normally resume checkpoints. Keep the frozen compressed ephemerides and the catalog at their relative paths. The F845M cache is shared by the two planets. `pack.py` verifies that the output represents every HST observation ID in the reconciled archive inventory before it writes the deterministic table bundle.

```sh
for body in Uranus Neptune; do
  python screen.py "$body" /path/to/rolling-cache /path/to/f845m-cache
  python stationary.py "$body" /path/to/rolling-cache /path/to/f845m-cache
  python sparse.py "$body" /path/to/rolling-cache
done
python pack.py
```
