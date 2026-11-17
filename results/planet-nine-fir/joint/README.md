# Joint optical recovery on the same synthetic planets

DES or Euclid detected 397 of 400 simulated dark, compact planets on trajectories drawn with the published centroids, and 388 of 400 after perturbing the two historic centroids. Three and twelve of the injections respectively escaped both searches. Each artificial planet has 5 Earth masses, radius 2.71 Earth radii, and geometric albedo 0.10 across the optical passbands. Distance, geocentric motion, satellite parallax and apparent brightness enter at each observation epoch. The two instruments use independent pixels, local backgrounds, source masks and detection rules. A trajectory missed by one survey is counted as recovered only if it passed the other survey's full search. [Complete matched results](../orbits/euclid/paired_summary.json) · [Per-trajectory outcomes](../orbits/euclid/paired_trials.csv.gz)

| Historic positions and model | DES | Euclid VIS | Both | Either | Conditional sensitivity, 95% Wilson interval |
|---|---:|---:|---:|---:|---|
| Published, 5 M⊕, 2.71 R⊕, albedo 0.10 | 379 | 362 | 344 | 397/400 | 99.25% (97.82–99.74%) |
| Perturbed, same dark compact model | 371 | 271 | 254 | 388/400 | 97.00% (94.83–98.28%) |
| Published, 5 M⊕, 2.94 R⊕, albedo 0.30 | 389 | 362 | 351 | 400/400 | 100% observed (99.05–100%) |
| Perturbed, same brighter model | 379 | 271 | 258 | 392/400 | 98.00% (96.10–98.98%) |
| Published, 10 M⊕, 3.70 R⊕, albedo 0.30 | 389 | 362 | 351 | 400/400 | 100% observed (99.05–100%) |
| Perturbed, same 10 M⊕ model | 379 | 271 | 258 | 392/400 | 98.00% (96.10–98.98%) |

The draw is uniform with replacement from a 13 × 32 × 56 grid of date, range and radial change, retaining bound solutions. For the stress experiment it also samples uniformly from 25 combinations of five discrete displacements per historic position. The data supply neither an orbital posterior nor a calibrated error distribution for AKARI. The intervals describe finite injection sampling. They exclude systematic errors in historical astrometry, planet photometry and image models. The separate continuously sampled Euclid experiment uses a different distribution of centroid shifts. [DES draw definitions](../archives/OPTICAL.md) · [Euclid injections](../orbits/euclid/SEARCH.md)

For the coarser 3,617-node central-position grid, at least one of DES or Euclid covered 3,601 nodes after masks. Sixteen failed both masks. Fixed-magnitude injections on 320 shared grid trials give a separate check without assuming a planet colour. Geometry and the measured magnitude pairs are reproducible via `python joint/join.py` from the campaign directory. [Geometry](matched_geometry.json) · [Magnitude table](matched_recovery.csv)

The physical experiment is reproducible by running `archives/planet_injections.py` after the DES pixel preparations, then `orbits/euclid/paired.py archives/planet_draws.csv.gz` after the Euclid pixel and empirical-PSF preparations. Retain each set of emitted files in its own survey directory. File paths and replay details appear in the linked survey reports.

The 2025 DECam narrowband exposure supplies a third test of the dark-model misses. Of three nominal and twelve perturbed misses, one and six were usable in N395. The highest real-source signal-to-noise among them was +1.02. No seeded synthetic planet passed. Sixty-four PSF/noise repeats at each usable location gave three recoveries on one perturbed draw and none elsewhere. These tests used the actual saved states and each draw's blue-band brightness. [Measurements and flags for all 15 draws](../n395_cal/paired_misses.csv) · [Calibration](../n395_cal/README.md)
