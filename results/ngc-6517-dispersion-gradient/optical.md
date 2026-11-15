# Stellar colors across NGC 6517

The radio fit predicts falling dispersion measure toward the east. In the same field three stellar colors increase eastward. All use proper motions to reduce contamination by Galactic field stars. This is a test of the sign expected if dust reddening and the column of free electrons trace the same monotonic spatial variation.

| Photometry | Color | Stars | Eastward slope (mag/arcmin) | Range omitting one of eight sky sectors |
|---|---|---:|---:|---:|
| Hubble Source Catalog v3, WFC3 | F555W − F814W | 230 | +0.086 | +0.065 to +0.092 |
| Pan-STARRS1 DR2 | r − i | 259 | +0.092 | +0.084 to +0.096 |
| Gaia DR3 | BP − RP | 331 | +0.107 | +0.090 to +0.119 |

The sectors are 45° wedges centered on the cluster, with boundaries at multiples of 45° east of north. The ranges describe sensitivity to omitting sky areas. The catalogs have different instruments and passbands but share stars, Gaia selection and the underlying sky structure.

## Observations and controls

The archive queries selected 27,480 Hubble matches, 2,294 Pan-STARRS objects, and 2,187 Gaia sources. The analysis uses an annulus of 10–85 arcsec about J2000 18h01m50.52s, −08°57′31.6″. Select Gaia sources within 1 mas/yr of (−1.554, −4.452) mas/yr, with RUWE <1.4 and both proper-motion uncertainties <0.6 mas/yr. Match each catalog to the nearest Gaia source within 0.35 arcsec, keeping the closest qualifying photometric match per Gaia ID. These are membership candidates. The proper-motion vector is from Vasiliev & Baumgardt's cluster measurements. [Cluster parameters](https://people.smp.uq.edu.au/HolgerBaumgardt/globular/orbits.html) · [Gaia cluster analysis](https://doi.org/10.1093/mnras/stab1475)

Fit color to east and north position, magnitude, magnitude squared, and cluster-centric radius using a robust loss. WFC3 and PS1 use the redder filter as the magnitude coordinate and Gaia uses G. The function in optical.py specifies the catalog flags, color/magnitude cuts and loss. The eastward coefficients remain positive with proper-motion radii of 0.5 and 2 mas/yr, with a tighter 0.15 arcsec match, inside 60 arcsec, outside 30 arcsec, and for every omitted sector. Color has units of magnitudes. Filter responses, crowding, sky coverage and sensitivity to the stellar sequence differ across catalogs. In particular, this measurement omits the central 10 arcsec where several of the pulsars lie. It compares directions across the field and does not assign reddenings to individual pulsars.

An earlier CFHT BVI analysis already reported differential reddening across a larger 2.35 × 3.76 arcmin field. It described the northeast as bluer than the southwest and a best gradient axis 20° west of north. Our measurements are local, depend on aperture, and use different passbands and samples. Neither a uniform relationship across the full historical field nor a dust-to-electron conversion follows from these data. The relative signs disfavor a single component in which increasing dust and free-electron columns dominate both gradients. Multiple foreground components, changes in ionization, or cluster material remain possible. [Kavelaars, Hanes & Bridges](https://arxiv.org/abs/astro-ph/9411090)

## Reproduction

Run `python optical.py` here with NumPy, pandas and SciPy. To re-query each archive, supply `--fetch /tmp/ngc6517-optical` with an empty cache. The local replay uses optical-matches.csv, including the looser proper-motion selection so the alternatives can be reproduced. The original query URLs, returned row counts, and exact receipt SHA-256 values for this run are in optical-inputs.json. Live services can change row ordering. The analysis orders the selected rows by catalog and Gaia ID. Retrieved 2026-11-15. [Hubble Source Catalog](https://archive.stsci.edu/hst/hsc/) · [MAST Pan-STARRS](https://catalogs.mast.stsci.edu/panstarrs/) · [Gaia DR3](https://www.cosmos.esa.int/web/gaia/dr3)
