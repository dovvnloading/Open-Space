# Testing planetary-nebula membership in nearby star clusters

An apparent excess of planetary nebulae near compact clusters in the bright PHANGS–MUSE sample fails to reproduce in 544 additional nebulae from the same galaxies. The closest globular-cluster candidate with agreement between the two nebular classifications, NGC 1433 cluster 883 and PN 89, remains consistent with a chance projection. The analysis establishes neither membership nor a progenitor mass.

## Population test

We compared the 899 Scheuermann planetary nebulae (PNe) in 19 galaxies with the union of the PHANGS-HST DR5 human and machine class 1+2 catalogs. We deduplicated identical pointing and cluster IDs, retaining the human catalog row where available. The combined source file contains 39,209 clusters before selecting the host galaxy. For each PN we counted at most one neighbor at a given radius. The table uses clusters with DR5 ages of 30–3,000 Myr and a projected angular radius of 0.45 arcsec.

| Nebular sample | Sources | Nebulae with a close cluster | Mean under shifts | Fraction of simulations with at least the observed count |
|---|---:|---:|---:|---:|
| Scheuermann PNe | 899 | 18 | 8.071 | 0.0018 |
| Scheuermann PNe within 0.5″ of a Congiu PN | 236 | 4 | 1.519 | 0.0612 |
| Additional Congiu PNe farther than 1″ from every Scheuermann PN or remnant | 544 | 1 | 5.534 | 0.9966 |

In each of 4,999 draws we translated the whole nebular pattern of each galaxy together, in a random direction 3–10 arcsec from the observed position, uniform in annular area. This preserves grouping of nebulae within each galaxy and samples local cluster density. Sources near edges can leave coverage. We examined four radii, three cluster selections and multiple catalogs. The fractions in the table are descriptive single-cut upper tails. Even multiplying the smallest quoted fraction by the 12 radius/cluster combinations gives 0.0216 for the primary sample alone. The second reduction uses the same MUSE exposures and a different extraction algorithm. Its additional objects test transfer across selections and have different detectability. A difference between the two catalogs does not identify the physical cause.

## Individual tests

| Host, nebula, HST cluster | Catalog separation | Independent HST narrow-band check |
|---|---:|---|
| NGC 1433, PN 89, cluster 883 | 0.172″ (15.8 pc at 18.94 Mpc) | A F657N residual peaks 0.161″ from the cluster center. The aperture sum is 3.42 local robust standard deviations above background. 1.75% of neighboring control positions yield an equal or larger aperture sum after a local peak search. |
| NGC 4303, PN 11, cluster 3040 | 0.137″ | The F658N check yields 1.00 local robust standard deviation. 90.4% of local controls equal or exceed the residual. |
| NGC 628, PN 99, cluster 767 | 0.097″ | The F658N check with a single continuum band yields 2.74 local robust standard deviations. 20.9% of local controls equal or exceed the residual. |

Congiu assigns PN model probabilities of 99.9% and 100.0% to the neighbors of PN 89 and PN 11 (as rounded in the published table). These probabilities apply to the emitting gas. In NGC 1433 the HST red continuum source is a class 1 machine selection flagged as a globular candidate by Floyd et al. The DR5 photometric age is missing. Its older DR3 age fit gives 6.69 Gyr. Among all 899 PNe, one lies within 0.2″ of a globular candidate. The shifts give at least one such match 38.2% of the time. HST independently records photons from the candidate's direction, but there is no cluster/nebula velocity comparison and the projected offset is large relative to a typical globular-cluster half-light radius.

The bundled 73×73 pixel images are native count-rate cutouts with embedded source URL, WCS, floating-point location and pixel origin. The broad-band mosaics come from PHANGS-HST and the narrow-band mosaics from HST skycells. The code aligns by WCS, robustly fits local continuum plus a constant in a 9–32-pixel annulus, searches within five pixels of the catalogued nebula and sums a 2.5-pixel aperture. Residuals can include [N II], registration or PSF mismatch, or diffuse emission. Control apertures overlap. Their empirical fractions do not represent independently calibrated false-alarm probabilities.

## Earlier M31 search

We also tested all 97,020 pairs in the 441-object Caldwell & Romanowsky M31 catalog. At an adopted 780 kpc, 35 lie within 150 projected pc. Of the 25 with both velocities and masses, none has a nominal absolute line-of-sight velocity difference below the upper bound on the two-body escape speed from the projected distance. The closest complete pair is B127-G185/NB89 at 15.712 pc, with a velocity difference of 190.4 km/s and escape speed of 42.64 km/s. This is a nominal-mass, nominal-velocity screen. It does not exclude all 35 pairs or include the host tidal field. An exploratory PHAT match also places the known M31 PN M06 3073 1.05″ from AP 2050 (2.5 half-light radii); the single projection and lack of cluster velocity did not justify promotion.

## Reproduction

Run with Python, NumPy, SciPy and Astropy. Catalog results were reproduced with Astropy 8.0.2. The external data services must be accessible for the first command.

```sh
python prepare.py /tmp/phangs-nebula-inputs
python analyze.py /tmp/phangs-nebula-inputs > /tmp/catalog-results.json
python pixels.py > /tmp/pixel-results.json
```

Compare these outputs with `results.json` and `pixels.json`. The input manifest gives URLs and SHA-256 hashes for the exact catalog content. VizieR hashes exclude changing comment and blank lines. The image bundles store the original archive URLs. The narrow-band FITS products are hst_skycell-p0365x13y06_wfc3_uvis_f657n_all (NGC 1433), hst_skycell-p1458x13y05_wfc3_uvis_f658n_all (NGC 4303), and hst_skycell-p1684x15y02_wfc3_uvis_f658n_all (NGC 628).

[PHANGS-HST compact-cluster catalog and release documentation](https://archive.stsci.edu/hlsp/phangs/phangs-cat) · [Scheuermann et al. 2022](https://academic.oup.com/mnras/article/511/4/6087/6510829) · [Congiu et al. 2023](https://doi.org/10.1051/0004-6361/202245153) · [PHANGS-HST imaging](https://archive.stsci.edu/hlsp/phangs/phangs-hst) · [Caldwell & Romanowsky 2016](https://doi.org/10.3847/0004-637X/824/1/42) · [Veyette et al. 2014](https://doi.org/10.1088/0004-637X/792/2/121) · [Johnson et al. 2015](https://doi.org/10.1088/0004-637X/802/2/127).
