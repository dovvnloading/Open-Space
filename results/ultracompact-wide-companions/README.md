# No high-quality wide companion passes the Gaia screen

We searched the 153 mass-transferring ultracompact white-dwarf binaries and candidates in the 2 March 2026 edition of the Green, van Roestel & Wong catalogue. Of 132 Gaia identifications, 95 had a proper motion and positive parallax with signal-to-noise at least two. The search finds no high-quality pair whose **nominal** relative tangential speed is below the escape speed for a total mass of two solar masses. The survey does not constrain unresolved companions, companions below Gaia's detection limit, or orbits hidden by velocity errors. [Catalogue](https://zenodo.org/records/18842123) · [Survey paper](https://doi.org/10.1051/0004-6361/202554925)

| Stage | Neighbor rows |
| --- | ---: |
| DR3 parallax entries in 95 cones, excluding each host | 44,259 |
| Projected separation at most 100,000 au | 44,201 |
| Both parallaxes at least 5σ, both RUWE below 1.4, neither flagged as duplicated | 5,673 |
| Parallaxes agree within 3 combined standard errors with a 0.025 mas floor in quadrature | 420 |
| Nominal projected speed within the escape envelope for two solar masses | 0 |

The closest case by speed ratio is ASASSN-21fhx with Gaia DR3 6917172068614964608. Its measured transverse difference is 0.475 km/s against a projected escape speed of 0.271 km/s at 48,239 au. Measurement uncertainties and projection of the system's space motion affect this comparison. We use the host distance for both stars and do not interpret this diagnostic as a probability of binding or a population limit. No positive triple-system claim follows. We selected the numerical thresholds during exploratory analysis.

[All cones with counts and checksums of the returned source-ID sets](fields.csv) · [All 420 parallax matches ordered by speed ratio](parallax_matches.csv) · [Pre-screen checkpoint](SEARCH.md)

With Python, NumPy, pandas and Astropy installed, run the following from this directory. Each query returns the same Gaia release from ESA or Heidelberg. Query radii allow 100,000 au with a minimum of 60 arcsec and a maximum of 0.5 degree. The source catalogue is pinned by its published checksum. The script explicitly reports and removes the host row before counting neighbors.

```sh
python search.py /tmp/ultracompact-gaia --output /tmp/ultracompact-replay
```
