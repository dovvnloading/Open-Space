# Literature and counterpart screening

Ten catalogue products were read in full. The comparison uses rows with parseable coordinates. The manifest gives download URLs, byte counts, and SHA-256 digests of the inputs. VizieR includes the query time in its headers, so replay also checks its header-independent payload digest.

| Input | Positions read | Source |
|---|---:|---|
| UltracoolSheet v2.1.0 main table | 3,890 | Best J2000 positions in [Zenodo 15802304](https://zenodo.org/records/15802304) |
| Kirkpatrick et al. 20 pc census | 3,058 | Positions with coordinates in [Table 4](https://exoplanetarchive.ipac.caltech.edu/docs/20pcCensus.html) |
| Backyard Worlds motion discoveries | 3,086 | CatWISE designations in [Schneider et al. 2026, Table 1](https://doi.org/10.3847/1538-3881/ae5b60) |
| SPHEREx SPIFF | 16,917 | All four groups, including published candidates, in [Gagné et al. survey tables](https://zenodo.org/records/19051216) |
| unTimely motion discoveries | 258 | Newly reported high-motion objects in [Karpov et al. 2025, Table 2](https://doi.org/10.1051/0004-6361/202452011) |
| Euclid Q1 spectroscopy | 201 | [Mohandasan et al. 2026](https://doi.org/10.1051/0004-6361/202660157) |
| Euclid Q1 photometry | 5,366 | Benchmark and candidate tables in [Zerjal et al. 2026](https://doi.org/10.1051/0004-6361/202556980) |
| Backyard Worlds companions | 187 | Components in [Rothermich et al. 2024](https://doi.org/10.3847/1538-3881/ad324e) |
| Backyard Worlds Spitzer follow-up | 93 | [Meisner et al. 2020, Table 6](https://vizier.cfa.harvard.edu/viz-bin/VizieR?-source=J/ApJ/899/123) |
| CatWISE cold-dwarf Spitzer follow-up | 169 | [Meisner et al. 2020, candidate table](https://vizier.cfa.harvard.edu/viz-bin/VizieR?-source=J/ApJ/889/74) |

Each infrared trial position is propagated along its measured tangent-plane motion to the published astrometric epoch when available. Euclid Q1 uses approximate epoch 2024.25. For entries with an uncertain reference epoch, a closest approach over 1985–2027 flags possible aliases. A 4″ match triggers review. A chance crossing or a nearby binary component can produce an alert. Conversely, an inaccurate motion, missing table entry or unsearched publication can hide a prior identification. An unmatched row has no automatic novelty classification. The search uses the W2 time series and SPHEREx observations to test it.

The Gaia search propagates the infrared position to J2016.0 and returns every Gaia DR3 source within 8″. The tables retain G magnitude, proper motion and parallax where available. A static neighbour can contaminate a WISE measurement while the cold dwarf itself is below Gaia's detection limit. The bright Euclid-field tier also has all Q1 MER sources within 12″ of its approximate 2024 position, with aperture photometry and morphology. The light of a nearby Euclid object can be blended within the WISE PSF. An empty Q1 cone can also fall outside the released mosaics.
