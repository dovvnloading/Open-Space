# Search for an unreported transiting exoplanet — research in progress

This investigation searches newly public TESS Sector 108 two-minute light
curves and checks surviving signals against independent observing sectors,
instrumental diagnostics, contaminating stars and prior candidate reports.
**No discovery claim is made by this initial checkpoint.** A periodogram trigger
does not by itself establish an exoplanet.

The frozen selection has 10,798 archive targets, complete TIC metadata for
those targets, and 1,062 eligible stars after the documented stellar and known
host cuts. `selection.csv` preserves every exclusion. `plan.json` records the
ordered target list and input snapshot identifiers. The broader literature
comparison remains a separate mandatory novelty check.

Read [PROTOCOL.md](PROTOCOL.md) for the question, selection, thresholds,
falsification tests and reporting standard. `search.py` downloads original
photometry, records full-product and scientific-column checksums, and performs
the exploratory BLS search. `vet.py` measures predicted events directly in
native samples with local baselines and reports coverage, odd/even differences,
secondary events, SAP/PDCSAP consistency, descriptive red-noise inflation,
background and centroid changes. Its phase ranks are conditional diagnostics,
not search-calibrated false-alarm probabilities. `catalogs.py` preserves exact
catalog bytes and query receipts.

## Reproduce the search

Use Python 3.12 with NumPy, SciPy, Astropy, pandas and requests. Keep the cache
outside the checkout. The initial plan was frozen before new photometry was
inspected; regenerating a plan with newer catalogs defines a different cohort.

```sh
python results/exoplanet-search-20261125/search.py inventory --cache /path/to/cache
python results/exoplanet-search-20261125/search.py metadata --cache /path/to/cache
# Retrieve the catalog snapshots using catalogs.py before making a new plan.
python results/exoplanet-search-20261125/search.py plan --cache /path/to/cache
python results/exoplanet-search-20261125/search.py fetch --cache /path/to/cache --limit 1062
python results/exoplanet-search-20261125/search.py screen --cache /path/to/cache --limit 1062
python results/exoplanet-search-20261125/search.py summarize --cache /path/to/cache
```

Final measured outcomes, rejection ledger, controls and retained candidate
evidence will replace this in-progress summary as the investigation completes.

Data credit: NASA/TESS, MAST/STScI, SPOC, TIC, NASA Exoplanet Archive, ExoFOP,
and the credited authors of the comparison catalogs. Public observations and
catalog misses do not establish priority; all candidate interpretations remain
subject to independent astronomical review.
