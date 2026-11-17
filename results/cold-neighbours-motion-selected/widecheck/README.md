# Six-pass unTimely W2 trajectory screen

**All 1,017 selected CatWISE rows in the 24 bright-search caps were tested using 792 epoch W2 catalogues across 132 tiles. Twenty-nine met the trajectory threshold. Twenty-four of those have positional alerts in the comparison literature products. The remaining five received native-image tests; none produced a repeatable motion residual in both scan seasons.** The six-visit counts are a triage statistic. Individual cold-object, temperature and distance claims require independent astrometric and spectroscopic evidence.

All 1,017 had six trial catalogue files. In 1,012 cases the local four-sigma depth proxy reached the CatWISE W2 magnitude in at least three chosen epochs. There were no nominal-track matches for 686 sources. Rotated 90°, 180° and 270° trajectories produced zero, one and three alerts, respectively. The [five-target native-image audit](../wide_deep/README.md) contains the held-out evidence and injection tests.

## Sample and measurements

For each CatWISE tile, six disjoint MJD windows target 2010–2011, early reactivation, mid-mission, late CatWISE, and two post-CatWISE periods. The exact boundaries are [55197,55750), [56580,57024), [57024,57754), [57754,58466), [58466,58850) and [58850,59216). Within each window the code prefers the alternating scan direction if a catalogue from that direction covers at least 95% of the tile, then prefers any catalogue meeting the coverage threshold. Among eligible alternatives it chooses the file with the most detections. At most one file per window enters the test. The count of local epochs with eligible neighbours is recorded separately. A tile-level footprint does not guarantee usable data on a particular moving source.

The primary measure counts distinct epochs with flux/error ≥4, quality factor >0.85 and a detection within 1.5″ of the predicted position. A stationary source and the same motions rotated 90°, 180° and 270° share the starting coordinate at MJD 57170. An alert requires at least three associated epochs over five years, at least one after CatWISE, and two more matches than the strongest alternative. The same rule is applied to each rotated trial for a negative control. A catalogue entry can satisfy several trial paths within a visit. For every source with at least three nominal moving matches, the offset test optimizes **both** the moving and stationary reference coordinates over the same 0.375″ grid inside a 3″ radius. Their maximum epoch counts are descriptive, without a calibrated false-alarm probability.

The depth proxy is the per-epoch median of 22.5 − 2.5 log₁₀(4 × dflux) from detections within 25″ of the reference coordinate with quality factor >0.85, positive uncertainty and at least three input frames. It can differ from the sensitivity at the target and supplies no completeness estimate.

## Replay and data

Run `python check.py` followed by `python counts.py` with NumPy and Astropy installed. The defaults read `../wide/selected.csv` and `../catalogues/wide_matches.csv`. Flags `--catalog`, `--literature`, `--cache`, `--state`, and `--workers` allow explicit input and temporary directories. Files are downloaded from NERSC, hashed, examined in memory and discarded unless they already exist in the optional cache. A compact per-tile checkpoint allows an interrupted run to resume. The full files can be retrieved from the exact URLs in `manifest.json`, which also records the input sample and catalogue-alert table digests. `summary.csv` lists nominal, rotated and free-intercept results. `detections.csv.gz` retains the individual local W2 extraction records. Spectra and native-image tests are reported in the main result and dedicated follow-up.

CatWISE and unTimely share WISE exposures through 2018. The last two windows use data collected after the final CatWISE exposures. The independent extraction algorithms therefore supply processing crosschecks at early epochs as well as added time baseline at late epochs.

[unTimely data, DOI 10.26131/IRSA580](https://www.ipac.caltech.edu/doi/irsa/10.26131/IRSA580) · [Meisner et al.](https://doi.org/10.3847/1538-3881/aca2ab)

The committed W2 neighborhood table is split into independently readable gzipped CSVs. Run `python unpack.py` to reconstruct the byte-identical merged table. [Chunk digests and merged SHA-256](observations.json)
