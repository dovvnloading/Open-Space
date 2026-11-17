# Replaying the measurements

The observing units are a unique Gaia DR3 source for the kinematic census, a survey observation for spectral fitting, and an explicitly simulated draw for injection recovery. Gaia positions use ICRS epoch J2016.0. Parallaxes are in mas, proper motions in mas/yr, source IDs are exact decimal integers, and speeds are in km/s. The code and archived observations were assembled by a research agent. Only the injection grid and uncertainty draws are simulated.

### Frozen inputs and software

The contribution builds on repository revision `f8ceda0c69a4a5a48b35faee0d5b68efc5c7982f`. Run from `results/gaia-dr3-supernova-survivors` in a checkout that also contains `results/dr1-dr19-unbound-stars`. The Git revision of this directory identifies the analysis code. `SHA256SUMS` covers the checked-in study files apart from the manifest itself. Recorded queries and source file receipts identify the external bytes. Re-fetching can change timestamps, header text and row order. Reproduce the frozen tables before refreshing them.

The numerical run used Python 3.12, NumPy 2.3.5, pandas 2.2.3, SciPy 1.17.0, Astropy 8.0.2, Gala 1.12.0, galpy 1.12.1, Matplotlib 3.10.8 and gaiadr3-zeropoint 0.1.0. The literature extraction uses Poppler's pdftotext. Archive access requires no account. Keep several GB of local scratch for both compressed and expanded SDSS FITS.

```sh
sha256sum -c SHA256SUMS
python control/recovery.py
python search/analysis.py
python search/inject.py
python plot.py
```

### Full catalogue, spectra and orbit runs

Set `SURVIVOR_CACHE` to a writable absolute path outside this result folder. The DESI, SDSS and union analyses share the cached Gaia rows. The release SHA-256s and the deterministic ordered batches are verified before the common cuts are applied.

```sh
python desi/fetch.py "$SURVIVOR_CACHE"
python desi/screen.py "$SURVIVOR_CACHE"
python sdss/screen.py "$SURVIVOR_CACHE"
python spectra/lamost_parent.py --cache "$SURVIVOR_CACHE/lamost" --output spectra
python union/census.py "$SURVIVOR_CACHE" "$SURVIVOR_CACHE/lamost/lamost_wd.csv.gz"
python origin/run.py
python origin/remnants.py
python literature/crossmatch.py
python literature/recent.py "$SURVIVOR_CACHE"
```

The detailed instructions include replay of the published [control and blend labels](control/README.md), [Gaia candidate and injection queries](search/README.md), [232 survey coverage queries, XP spectra and Mg/Balmer fits](spectra/README.md), [original Shen spectra](shen/README.md), [DESI metal and carbon selections](metals/README.md), [SDSS classifications and line fits](sdss_followup/README.md), and [hot-subdwarf RV selections](hot/README.md). Air/vacuum wavelengths, observatory frames, masks, per-observation IDs and error models are specified with each fitter. Formal line-fit errors do not account for all pressure shifts, gravitational redshifts or binary motion. Atmospheric fits to XP inherit assumptions about ordinary Galactic stars.

### Archives and citations

This analysis uses data from the European Space Agency's Gaia mission, processed by the Gaia Data Processing and Analysis Consortium, including the DPAC BP/RP products, Lindegren parallax correction and Bailer-Jones distances. It uses DESI DR1 and the work of Swan et al., the Sloan Digital Sky Survey and its Astra pipeline, the Guo Shou Jing Telescope (LAMOST), the SIMBAD and VizieR services at CDS, the GAVO data centre, the NOIRLab Astro Data Lab, NASA/IPAC IRSA, and the original observing teams' released spectra. Each directory links its original data, publication and release. [Gaia DR3](https://doi.org/10.1051/0004-6361/202243940) · [DESI collaboration and acknowledgments](https://data.desi.lbl.gov/doc/acknowledgments/) · [SDSS](https://www.sdss.org/collaboration/citing-sdss/) · [LAMOST](https://www.lamost.org/lmusers/) · [SIMBAD](https://simbad.cds.unistra.fr/simbad/) · [VizieR](https://vizier.cds.unistra.fr/)
