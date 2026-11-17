# Bright infrared-motion search over 4,235.10 square degrees

Across 24 disjoint caps, 2,829 CatWISE2020 rows meet the wide query constraints. The full bright selection retains 1,017 entries with W2 ≤15.5. The literature match search flags 29 and leaves 988 for independent follow-up. These numbers count catalogue entries. A nearby published source can produce several fits, and a compact patch of structured emission can produce many entries.

Each 7.5°-radius field was selected by a reproducible equal-area placement in latitude with successive golden-angle longitudes. No candidate positions enter the choice of fields. All 24 queries completed. Each response contains fewer than 50,000 requested records, with unique source IDs. The exact ADQL, raw byte hashes and row counts for every cap are in inputs_manifest.json. The caps overlap other parts of this investigation, so their area is reported separately.

The selection uses the pilot motion, colour, mask and fit thresholds, with a brighter W2 cutoff appropriate for checking single exposures. The thresholds are catalogue-level gates. Background structure, correlated scans, negative W1 measurements, reference-frame shifts and centroid switching limit interpretation. At 20 pc, even a source with WISE 0855's W2 luminosity would be much fainter than this cap search.

~~~bash
python wide/fetch.py
python wide/select_wide.py
python catalogues/crosscheck.py wide/selected.csv catalogues/wide_matches.csv
python counterparts/gaia_wide.py
~~~
