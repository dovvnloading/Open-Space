# All four apparently fast SOHO comets admit ordinary near-Sun paths

A heliocentric parabolic orbit fits every recorded position of each of the four rapid batch outliers. The unweighted per-coordinate RMS ranges from 22.47″ to 37.61″. Each value is lower than the corresponding unrestricted three-second orbit in the catalogue replay. This comparison uses all 123 directions and the 123 accompanying spacecraft positions. The observations do not require a strongly hyperbolic trajectory.

| Object | Positions | Arc (days) | Quick-fit nominal excess speed (km/s) | Quick-fit RMS (″) | Rounded published orbit RMS (″) | Independently fitted e=1 RMS (″) | Fitted q (au) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| C/2004 D4 | 49 | 1.3383 | 47.3 | 98.4 | 48.40 | 22.47 | 0.00498 |
| C/2005 D3 | 21 | 0.6622 | 133.0 | 68.18 | 43.13 | 28.17 | 0.00698 |
| C/2006 S8 | 34 | 0.7167 | 72.6 | 40.10 | 90.60 | 35.49 | 0.00823 |
| C/2006 S9 | 19 | 0.4333 | 325.4 | 43.9 | 74.53 | 37.61 | 0.00502 |

The [independent orbit fitter](../soho_models.py) parses every optical direction and its paired MPC code 249 spacecraft vector. Adding that vector to ERFA's heliocentric Earth position gives the observer coordinate at each exposure. SOHO's geocentric range was 0.01084–0.01123 au across the sample. The prediction uses a two-body parabolic solar orbit at emission time with iterated light time. We vary perihelion distance, perihelion time, inclination, ascending node and argument of perihelion at fixed eccentricity, with equal weight per coordinate. The five parameters can vary substantially on a short arc. We use this calculation to test whether the extreme speed is necessary for a fit. It does not estimate an original barycentric energy or the probability of a particular energy.

The angular RMS is the square root of the sum of squared east and north residuals divided by twice the number of positions, matching [the pinned Find_Orb routine](https://github.com/Bill-Gray/find_orb/blob/1de55cca12e7a585c9c7384acbb80538ea2afc3c/orb_func.cpp#L950). An unconstrained global least-squares optimum could always fit at least as well as a nested restricted model with the same residual definition and data. The three-second search's larger residual shows why its local formal uncertainty cannot decide the energy. The machine-readable [residuals, elements, input hashes and MPC receipts](results.json) accompany the four complete, unedited designation-specific observation excerpts in this directory.

## Image and orbit history

Each object already has identified coronagraph images, a published short observed arc and a known Kreutz-group classification. Planetary perturbations from a long solar-system passage cannot be inferred from these hours of near-Sun imaging. The comet's concentration of light can move as the nucleus loses material. An orbit formal error from a poorly fitting short-arc solution does not incorporate that systematic.

| Object | Mission and source images | Contemporary photometric evidence | Previous independent orbit fit |
| --- | --- | --- | --- |
| C/2004 D4 | SOHO-744, J. Sachs, LASCO C3 and C2 | Knight and colleagues measured peak C3 magnitude 4.5 at a heliocentric distance of 12.9 solar radii | 49/49, mean residual 22.248″, fitted e = 0.999120 ± 0.000836 |
| C/2005 D3 | SOHO-910, K. Černis, LASCO C3 | The original project pixel list spans 06:42:05–22:35:41 UTC on February 22 | 20/21, mean residual 25.745″, e fixed to one at the constraint epoch |
| C/2006 S8 | SOHO-1192, R. Matson, LASCO C3 and C2 | In C3 the object was elongated and reached approximately magnitude 4.5. It appeared diffuse in C2, showed a faint tail and faded rapidly. | 33/34, mean residual 33.924″, fixed e = 1 |
| C/2006 S9 | SOHO-1189, R. Matson, LASCO C3 and C2 | Approximately magnitude 6 in C3. It appeared diffuse and faint in C2 with a possible tail. | 18/19, mean residual 34.613″, fixed e = 1 |

The mission's year tables give group membership and discovery metadata. The ongoing SOHO project lists perihelion distance 0.0048 au for S8, whereas the frozen SBDB/MPC SAO_2008 orbit gives 0.0093 au. The two catalogs used different solutions. The fitted epoch eccentricity printed for D3 by Project Pluto is 1.0044642 with a displayed zero error after imposing its constraint. We draw no original-energy inference from this number. [NRL Kreutz tables for 2004](https://sungrazer.nrl.navy.mil/comets_table_2004), [2005](https://sungrazer.nrl.navy.mil/comets_table_2005), [2006](https://sungrazer.nrl.navy.mil/comets_table_2006). [D3 contemporary image centroids](https://sungrazer.nrl.navy.mil/sites/sungrazer/files/comet_txts/soho0910_xy.txt). [Knight et al., Table 2](https://pages.astro.umd.edu/~dphamil/research/preprints/KniAHeBie09.pdf). [IAUC 8772 observer reports archived by the BAA](https://people.ast.cam.ac.uk/~jds/kreutz06.htm).

The prior independent numerical fits and their observation selections appear in the [Project Pluto C/2004 D4](https://www.projectpluto.com/soho/2004D4.htm), [C/2005 D3](https://www.projectpluto.com/soho/2005D3.htm), [C/2006 S8](https://www.projectpluto.com/soho/2006S8.htm) and [C/2006 S9](https://www.projectpluto.com/soho/2006S9.htm) pseudocirculars. The original discovery measurements were published in [MPEC 2004-M14](https://www.minorplanetcenter.net/iau/mpec/K04/K04M14.html), [MPEC 2005-G40](https://www.minorplanetcenter.net/iau/mpec/K05/K05G40.html) and [MPEC 2006-U63](https://www.minorplanetcenter.net/iau/mpec/K06/K06U63.html).
