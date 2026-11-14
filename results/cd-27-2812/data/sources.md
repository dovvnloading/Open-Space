# Data provenance

Measurements were retrieved on 2026-11-14. Gaia source identifiers are decimal strings throughout. The astrometric replay uses catalog parallaxes directly without zero-point or covariance corrections between sources. The transverse speed is a projected quantity.

- `hosts.csv.gz` contains Gaia DR3 fields for all 1,312 `GaiaDR3` IDs from VizieR `J/A+A/699/A153/table2`. `wdms-row.csv` preserves the selected row using the VizieR column names. The 2025 catalog uses low-resolution spectra, broadband constraints and source-quality checks. The published component parameters are available for only a subset. [VizieR table](https://vizier.cfa.harvard.edu/viz-bin/VizieR?-source=J/A%2BA/699/A153/table2)
- `cone-returns.csv.gz` contains the 2,320 returned rows, including hosts, tagged with their zero-based query group. Sort the host IDs lexicographically as decimal strings and group them by nines. For each host, place a Gaia DR3 cone of radius `min(300, max(20, 25 * host_parallax_mas))` arcsec about its DR3 coordinates, then apply the ADQL prefilters below. Combine the nine predicates by OR. Cone centers and radii in degrees were rounded to seven decimal places, and astrometric filter centers to three. Final velocity and volume tests are in the offline analysis. The escape scale used for all screening is computed for a nominal two solar masses. [Gaia DR3](https://doi.org/10.1051/0004-6361/202243940)
- `gaia-pair.csv` contains the astrometry, photometry and transit quality indicators for both Gaia DR3 IDs from `gaiadr3.gaia_source`, and `nss-orbit.csv` contains selected columns for the bright source in `gaiadr3.nss_two_body_orbit`. All Gaia TAP queries used `https://gea.esac.esa.int/tap-server/tap/sync` with `REQUEST=doQuery`, `LANG=ADQL`, `FORMAT=csv`.
- `control-field.csv.gz` keeps the IDs, positions, parallaxes and proper motions for every `gaiadr3.gaia_source` returned by the field query below. `comparison-colors.csv.gz` keeps source IDs, red magnitudes, colors and parallaxes from the complete nearby comparison returned by `comparison-query.adql`. Its absolute RP range is 9.94–10.34 mag. The comparison is empirical and includes any objects satisfying the cuts, without a spectral-type label or extinction correction.
- `2mass.tsv`, `catalog-wide-pair.tsv`, `msc.tsv`, and `wds.tsv` are the archived VizieR ASCII query responses. Each header includes its full request URL. The two empty catalog queries describe their tables and contain no source rows. The confidence statistic `R` in the wide-pair table is a density ratio calibrated against chance alignments. [El-Badry et al. 2021](https://doi.org/10.1093/mnras/stab323)

The following cone predicates were applied in addition to the positional constraint. The proper-motion differences are per coordinate at this stage, in mas/yr.

```sql
ABS(parallax - host_parallax) < 0.5 + 3 * parallax_error
AND ABS(pmra - host_pmra) < 3 + 3 * pmra_error
AND ABS(pmdec - host_pmdec) < 3 + 3 * pmdec_error
```

The half-degree field query selected `source_id,ra,dec,parallax,pmra,pmdec` and extra quality columns from `gaiadr3.gaia_source` under this predicate.

```sql
1 = CONTAINS(POINT('ICRS',ra,dec), CIRCLE('ICRS',93.24857865,-27.88041694,0.5))
AND parallax > 0 AND parallax_over_error > 5
```

Gaia measurements come from the European Space Agency mission Gaia and the Gaia Data Processing and Analysis Consortium (DPAC). Funding for DPAC has been provided by national institutions participating in the Gaia Multilateral Agreement. This investigation used the CDS VizieR and SIMBAD services and data products from the Two Micron All Sky Survey, a joint project of the University of Massachusetts and IPAC/Caltech funded by NASA and the NSF.
