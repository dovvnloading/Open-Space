# OGLE-LMC-ECL-13803 depth prediction

The two-day orbit has catalog ephemeris HJD 2460001.0230 + 1.9964371 × cycle. The OGLE-IV catalog coordinates are 05:22:02.52, −69:08:09.7. The parent [science note](../../../docs/ARCHIVAL_FOLLOWUP_CONTROLS.md) gives the scope, coverage, result and interpretation. The analysis started from Open-Space revision 5f8bc79e188ff1940d87a59e43fe2e96a9083059. [OGLE-IV archive](https://www.astrouw.edu.pl/ogle/ogle4/OCVS/lmc/ecl/)

## Replay the observation

Run these commands in this directory with NumPy and Python 3.11 or later. The first regenerates holdout.json using the target's released OGLE-IV I magnitudes and the original catalog row. The second regenerates lmc-early.json with the catalog ephemeris, earlier OGLE-III I-band data, and EROS-2 source lm0010m13681 at a 0.420-arcsecond offset.

```sh
python holdout.py
python lmc_early.py
```

OGLE-IV times use full HJD. OGLE-III and EROS-2 report HJD minus 2,450,000. EROS R and B are measured in their original passbands. The historical test rejects uncertainties above 0.10 mag for OGLE-III and 0.5 mag for EROS. The compressed EROS input retains all 654 original target rows from the coordinate query. The EROS catalog receipt records nearby sources. Source URLs and SHA-256 values for the original downloads are in [sources.json](sources.json). A fresh VizieR response may have a different timestamp header.

## Replay the population screen

The signed prediction is preserved in [prediction.txt](prediction.txt). The 1,280 training ranks are in [screen.csv.gz](screen.csv.gz). For a full population rerun, place the four OGLE-IV release files ecl.dat, ident.dat, remarks.txt and phot.tar.gz in a separate directory, e.g. a scratch directory named lmc-screen/raw/. The photometry archive is approximately 1 GiB. SHA-256 values and download URLs are in sources.json. From this code directory, run the following with the scratch directory path adjusted as needed.

```sh
python screen.py --cache /absolute/path/to/lmc-screen
```

The screen writes selected light curves and a fresh screen.csv to the specified cache. Its ranking requires at least eight different observing nights in each selected phase window in both training eras, and at least twenty off-eclipse samples per era. The fixed phase half-width of each eclipse core is 0.03 cycles. Two quadrature windows collect the local comparison measurements. Nightly medians give each night equal weight.

Original photometry and catalog credit go to OGLE and EROS-2. [OGLE-III LMC catalog](https://www.astrouw.edu.pl/ogle/ogle3/OIII-CVS/lmc/ecl/) · [EROS-2 release](https://doi.org/10.1051/0004-6361/202557321) · [EROS-2 VizieR catalog](https://cdsarc.cds.unistra.fr/viz-bin/cat/II/390)
