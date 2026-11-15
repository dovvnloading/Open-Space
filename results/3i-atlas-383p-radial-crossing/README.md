# 383P/Christensen crossed the heliocentric direction of interstellar comet 3I/ATLAS

On 4 October 2025, 3I/ATLAS crossed the heliocentric J2000 ecliptic direction **230.860° longitude, +4.886° latitude**. Periodic comet 383P/Christensen crossed that same direction farther from the Sun about a week later. JPL Horizons places the bodies at 1.63547 and 2.87254 au at their respective crossings. A parcel moving at constant radial speed would need **301.17 km/s** to cover the intervening **185.06 million km** in **7.1121 days**. This identifies the specific receiving comet and the transport condition. Neither a nucleus flyby nor a measurement of interstellar material follows from the geometry. The simultaneous nucleus separations at the two JPL epochs were 192 and 210 million km.

| Orbit calculation | 3I source crossing | 383P receiver crossing | Radial gap (au) | Required speed (km/s) |
|---|---|---|---:|---:|
| JPL Horizons, DE441 | Oct 4 10:44 | Oct 11 13:26 | 1.237063 | 301.17 |
| MPC elements, eight-planet integration | Oct 4 10:43 | Oct 11 13:30 | 1.236970 | 300.99 |
| LTE/IMCCE Miriade, INPOP19A | Oct 4 11:23 | Oct 10 01:35 | 1.236283 | 382.79 |

Dates are in 2025, using TDB for Horizons, TT for Miriade and ephemeris days for the MPC integration. TT and TDB differ by under two milliseconds here. The MPC receiver epoch agrees with JPL to four minutes. Miriade's 383P position differs by 2.53 million km on Oct 11 at 00:00, predominantly along the orbit. That discrepancy moves the receiver epoch by 35 h 51 min while preserving the same directional crossing and nearly the same radial gap. This is a service-to-service discrepancy, not a confidence interval. The calculations can share underlying astrometry. JPL's 383P solution #21 includes 349 observations through July 2026. Its 3I solution #55 includes 117 observations with an Oct 2025–Oct 2026 arc. LTE reports a 2019 perihelion epoch for its 383P elements. Remaining differences between orbital solutions were not attributed to a particular cause.

## How the pair was selected

The exploratory parent search examined comet–comet close approaches and radial connections across many bodies. It found numerous candidate alignments. We selected this pair after seeing the identity of the source. The following focused screen reports its denominator and a rotation control. It is not a blind interstellar-object test.

The frozen MPC CometEls product contained 973 rows. We selected rows with perihelion distance at most 4.5 au and removed designations with a fragment suffix or a hyphen in the numbered prefix. This leaves 743 distinct catalogue labels, including 3I. Alternative provisional and permanent labels may still refer to the same physical object. For each of 742 other rows, both mutual nodal directions of its fixed osculating orbital plane and 3I's were examined. A true candidate has an outbound radial parcel with the receiving comet farther from the Sun, receiver radius below 4 au, source radius above 0.1 au, positive gap below 1.5 au and inferred speed from 250 through 850 km/s. Epochs were bounded by 2010-01-01 and 2046-01-01. The 3I source dates are fixed by its single unbound passage. We kept 383P alone, with a first-pass speed of 301.60 km/s.

Rotating only 3I's node by each integer angle from 1° to 359°, with all other elements, time intervals and filters fixed, yielded 284 orientations with zero matches, 67 with one, seven with two and one with three. Thus **75/359 (20.9%)** had at least one match. Adjacent rotations are correlated. The control measures how often the search geometry admits a result under longitude rotation. It gives neither a chance-coincidence probability for this selected pair nor evidence for a population excess. Screening used two-body orbits and can miss candidates requiring a perturbed orbit calculation. Only the selected pair received the three-solution check.

## A same-day wind measurement challenges the transport interpretation

MAVEN sampled the solar wind at Mars on 4 October. Mars was **0.20346 au and 6.46°** from 3I as seen from the Sun at the JPL source epoch. We retained valid three-component SWIA velocity samples on the dayside (MSO X > 0), moving away from the Sun, at least 0.15 Mars radii beyond the Vignes nominal bow shock along the ray from its focus (X₀ = 0.64, L = 2.04, eccentricity = 1.03). The quality flags for all three components had to be one. The 3,679 selected eight-second samples span fourteen UTC hours. Their outward MSO component ranged from **407 to 651 km/s**, with median **569** and 10th–90th percentiles **527–603**. Both required speeds fall below this entire sample. The hours bracketing the source passage with selected data, 09 and 12 UTC, had medians of 601 and 599 km/s. These are observations at Mars, on another radial path. No wind measurement at the 3I source point or in-situ compositional measurement at 383P was recovered.

The principal uncertainty in an ion-tail interpretation is transport. A constant 1° transverse deviation over the JPL radial gap shifts an endpoint by **3.23 million km**; 0.1° shifts it by 323,000 km. No local tail width, ion survival, field configuration or actual parcel trajectory was measured. The root is obtained by solving for two epochs, so the near-zero directional residual supplies no independent evidence for plasma transfer. Mars imagery establishes that 3I had a coma during this period. The data here do not establish how far that coma's ions propagated along this direction. [NASA observations](https://science.nasa.gov/solar-system/comets/3i-atlas/) · [Bow shock geometry used in MAVEN analysis](https://lasp.colorado.edu/maven/files/2017/02/Martian-electron-foreshock-from-MAVEN-observations.pdf)

## Prior work and reproduction

Directional alignment with a retarded launch time is a published ion-tail search technique. MSSL's Tailcatcher searches spacecraft trajectories. Earlier work on 3I explicitly predicted spacecraft tail encounters. Searches on 15 November 2026 for the paired identifiers “3I/ATLAS” and “383P/Christensen”, including “October 2025”, “tail”, “solar wind” and “Christensen”, located no report of this particular receiving comet. Coverage of search indexes is incomplete. [MSSL Tailcatcher](https://www.ucl.ac.uk/mathematical-physical-sciences/mssl/research/solar-system/planetary-science/tailcatcher-service) · [Prior spacecraft analysis of 3I](https://arxiv.org/abs/2508.15768)

Run these commands from this directory with Python, NumPy and SciPy. The first two operate offline. The third streams the 44 MiB daily MAVEN file, checks its SHA-256 and reproduces the committed hourly summary.

```sh
python screen.py
python measure.py
python wind.py
```

The two-time fit equates the directions of the normalized heliocentric positions along two tangent-plane axes and verifies that the unit vectors share the same hemisphere. Hermite splines interpolate the service positions and velocities. The MPC check starts from its two elements at their printed epochs, uses DOP853 with direct and indirect point-mass gravity of all eight planetary system barycenters and interpolates their JPL vectors at four-day spacing. It omits nongravitational forces and relativity. See [source receipts](sources.json) for queries, file roles and hashes. [JPL Horizons API](https://ssd-api.jpl.nasa.gov/doc/horizons.html) · [LTE Miriade API](https://ssp.imcce.fr/webservices/miriade/api/ephemcc/) · [MPC catalogue](https://minorplanetcenter.net/iau/MPCORB/CometEls.txt)
