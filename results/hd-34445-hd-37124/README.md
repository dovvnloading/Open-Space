# A predicted flyby between the planetary systems HD 34445 and HD 37124

HD 34445 and HD 37124 are currently 16.99 pc apart. Their Gaia space motions give a nominal closest approach of **0.0477 pc in 205,350 years**, at **80.9 km/s** relative to one another. The median across 200,000 draws from the Gaia measurement uncertainties is **9,968 AU**, with a central 90% interval of **7,542–16,800 AU**. Substituting Keck velocities for the Gaia velocities gives a nominal separation of **0.0474 pc**.

The NASA Exoplanet Archive lists six planets around HD 34445 and three around HD 37124. Their largest listed semimajor axis is 6.36 AU. Five of the HD 34445 planets are disputed. Rosenthal et al. linked signals to stellar activity and instrumental effects, and Howe et al. retained only planet b. The stellar encounter calculation uses the host stars’ motions and does not depend on those planetary fits. These data provide no evidence of planet disruption, exchange, or a comet cloud around either star. [Frozen planetary parameters](planets.csv) · [Rosenthal et al., Appendix G.4](https://doi.org/10.3847/1538-4365/abe23c) · [Howe et al., §2.1](https://ntrs.nasa.gov/citations/20250006033) · [HD 34445](https://exoplanetarchive.ipac.caltech.edu/overview/HD%2034445) · [HD 37124](https://exoplanetarchive.ipac.caltech.edu/overview/HD%2037124)

## How the pair was found

An exploratory 50 pc search compared known planet hosts with field stars. The two stars appeared as each other's perturber. We then compared every pair of qualifying planet hosts within 100 pc to establish their ranking. The second comparison includes the discovery pair and therefore measures its rank in the larger sample. It supplies no independent significance test.

The Gaia DR3 query required parallax ≥10 mas, parallax divided by its uncertainty >20, RUWE <1.4 and a radial velocity. This returned 133,468 stars. An exact Gaia DR3 source-ID join to the NASA Exoplanet Archive composite table on 2026-11-14 left 948 distinct host sources. The frozen table includes every one of those hosts. An exact identifier match excludes hosts lacking a populated DR3 ID. Both the magnitude limits on Gaia velocities and the planet survey selection affect the sample. [All 948 hosts](hosts.csv) · [Input queries and SHA-256 hashes](provenance.json)

All 448,878 unordered pairs were tested. We retained flyby candidates with relative speed above 5 km/s and approach epoch within ±2 Myr of J2016.0. Of the 416,117 pairs passing those cuts, 154 had a nominal miss under 1 pc and one fell inside 0.1 pc. The runner-up, GJ 1061 and ε Indi A, had a nominal distance of 0.1050 pc. This ordering uses catalog central values. Sample ranks and thresholds do not give the probability that a particular pair will encounter. [Computed rankings and measurements](results.json)

## Measurements and checks

For each star, we converted inverse parallax, sky position, tangential motion and radial velocity into Cartesian position and velocity. For relative position **r** and velocity **v** the time to closest approach is −(**r**·**v**)/|**v**|². At that time the distance is |**r** + **v**t|. Statistical draws use each source's full five-parameter Gaia covariance and a separate Gaussian draw for radial velocity. Stars' errors are sampled independently. The central 90% intervals are 204,229–206,473 years and 80.66–81.13 km/s. In 99.253% of draws the distance was under 0.1 pc. This is conditional on the adopted measurement model.

| Inputs or physical model | Nominal distance (pc) | Time ahead (kyr) |
| --- | ---: | ---: |
| Gaia DR3 astrometry and velocities | 0.0477 | 205.35 |
| Earlier Gaia DR2 astrometry and velocities | 0.0540 | 208.36 |
| Gaia DR3 astrometry and Chubak et al. Table 3 velocities | 0.0458 | 205.47 |
| Gaia DR3 astrometry and Keck velocities using the Table 1 present-work value for HD 37124 | 0.0474 | 205.43 |
| Hipparcos–Gaia long-baseline proper motions, Gaia DR3 distances and Keck velocities | 0.0497 | 205.42 |
| Hipparcos-only proper motions, Gaia DR3 distances and Keck velocities | 0.0353 | 205.49 |
| Hipparcos-only astrometry and Keck velocities | 0.5644 | 188.67 |
| Gaia DR3 propagated in a Galactic potential | 0.04685 | 205.35 |
| Gaia DR3 with differential light travel time | 0.0489 | 205.22 |

Times refer to each catalog's reference epoch (J2016.0 for DR3, J2015.5 for DR2 and J1991.25 for Hipparcos). The older Hipparcos-only parallaxes were 21.52 ± 0.66 and 29.70 ± 0.70 mas. DR3 measured 21.837 ± 0.021 and 31.609 ± 0.024 mas respectively. Hipparcos alone does not resolve a 0.05 pc miss. Earlier Gaia releases and long-baseline Hipparcos–Gaia proper motions share observations with DR3. [Hipparcos catalog](https://doi.org/10.1051/0004-6361:20078357) · [Brandt, long-baseline proper motions](https://doi.org/10.3847/1538-4365/abf93c)

Chubak et al. list −78.867 km/s for HD 34445 and −23.054 km/s for HD 37124 in Table 3, with 57 and 27 spectra over 1,973 and 1,990 days respectively. For HD 37124, Table 1 distinguishes the present Keck/HIRES measurement (−23.032 km/s) from the earlier Nidever measurement (−23.076 km/s). Their average equals both the Table 1 final value and the Table 3 entry. The Keck and Hipparcos checks above use the present-work value for HD 37124 and the Table 3 value for HD 34445. The Keck and Gaia spectra were taken with different instruments. Keck adopted the Nidever velocity zero point, and Gaia set its zero point using ground-based standard stars. Both spectroscopic estimates can share stellar line shifts. [Chubak et al., Tables 1 and 3 and §4](https://arxiv.org/abs/1207.6212) · [Gaia DR3 spectroscopic calibration standards](https://gea.esac.esa.int/archive/documentation/GDR3/Data_processing/chap_cu6spe/sec_cu6spe_input/ssec_cu6spe_auxdata.html)

The Gaia DR3 radial velocities are −78.940 and −23.078 km/s, from 24 and 32 transits. RUWE is 0.91 and 1.03 respectively. The Gaia radial-velocity variability p-values are 0.81 and 0.56. The cataloged speeds are spectroscopic estimates of the stellar center-of-mass motion. Convective line shifts, gravitational redshift and orbital reflex motion can introduce offsets. Independently moving the two Gaia velocities by −0.5, 0, or +0.5 km/s gives 0.0398–0.1479 pc across the nine trials. Independently offsetting their parallaxes by −0.03, 0, or +0.03 mas gives 0.0349–0.0730 pc. These are discrete systematic sensitivity tests. [Gaia follow-up rows](gaia_pair.csv)

We also integrated the nominal orbits in a Miyamoto–Nagai disk plus spherical bulge and halo using the masses and scale lengths of Fernandez-Puig et al., Table 2. The local initialization places the Sun 8.2 kpc from the Galactic center and 20.8 pc above the plane, with Galactic velocity (11.1, 245.0, 7.25) km/s. The tidal correction changes the miss by 0.00087 pc. For the separate light-time check, we advanced each star from the emission epoch to a common epoch by its distance divided by the speed of light and converted apparent transverse velocities to emission-time velocities using the factor 1 + vᵣ/c. For two solar masses, the focusing scale GM/v² is approximately 0.27 AU at the nominal relative speed. [Galactic model](https://arxiv.org/abs/2605.16496) · [Light-time equations](https://doi.org/10.1051/0004-6361/201424483)

## Prior work and falsification

Searches on 2026-11-14 for both HD names, their Hipparcos designations and Gaia IDs alongside encounter, flyby and closest-approach terms found the planetary-system papers and related encounter surveys. We located no earlier prediction for this pair. The recent all-star neighborhood study centered its star–star search on objects within 25 pc of the Sun. These stars are 45.80 and 31.64 pc away. The debris-disk and Solar encounter searches use different target samples. A literature search cannot establish absolute priority. [Fernandez-Puig et al.](https://arxiv.org/abs/2605.16496) · [Bertini et al.](https://doi.org/10.1051/0004-6361/202245415) · [Bailer-Jones](https://arxiv.org/abs/2207.06258)

The prediction can be tested with independently improved absolute parallaxes, proper motions and systemic velocities. The frozen Gaia model predicts the two-star miss will remain under 0.1 pc, and gives its expected time and distance ranges above. A later solution outside that distance or those ranges would revise the geometry. The wider systematic sensitivity range describes how such a revision could arise.

## Reproduce

From this directory, with Python 3.10+, NumPy and SciPy installed, run the following command. It verifies the frozen input hashes, screens all pairs, reruns the uncertainty and physical-model checks and prints the results. The code reads only local inputs. The upstream archive receipts include the full queries and the original retrieved byte hashes. The VizieR receipts include a request timestamp, so a later identical request can have a different byte hash. To reconstruct the 948-host table from the two upstream survey responses, strip the literal prefix `Gaia DR3 ` from the NASA ID, retain the first hostname for each source ID, inner-join to the Gaia rows by the exact integer ID and sort by that integer. [Replay code](replay.py) · [Frozen output](results.json)

```sh
python replay.py
```

This work used data from ESA's Gaia mission, processed by the Gaia Data Processing and Analysis Consortium, the ESA Hipparcos mission, the NASA Exoplanet Archive operated by Caltech under contract with NASA, the Keck Observatory and the CDS VizieR catalog access tool. [Gaia acknowledgements](https://www.cosmos.esa.int/web/gaia/dpac/consortium) · [NASA Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu/)
