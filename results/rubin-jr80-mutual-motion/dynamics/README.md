# A bound on system mass from two dated relative positions

The four-exposure Hubble mean on June 13 and the first Rubin *i* exposure on June 27 imply a **conditional system-mass lower bound of 1.71 × 10¹⁹ kg**. It assumes a gravitationally bound Kepler pair, that the measured relative positions trace the two bodies, and radial error envelopes of 0.02″ for HST and 0.06″ for Rubin. Allowing 0.10″ at both epochs reduces the bound to 8.33 × 10¹⁸ kg. These radii are explicit sensitivity choices. They have no calibrated confidence level. The geometry alone does not measure an orbital period, density, individual component mass or three-dimensional separation.

## Calculation

Let **r₁**, **r₂** be secondary-minus-primary vectors, R the sum of their lengths, c the endpoint chord and Δt the emission-time difference. Among bound transfers for the same endpoints and mass, the short-way parabolic limit gives the minimum time. Consequently the total mass must exceed

    [((R + c)^(3/2) − (R − c)^(3/2)) / (6 Δt)]² / G.

This is the finite-time Lambert relation. It does not assume a circular orbit or identify the mean projected chord speed with an instantaneous speed. The expression increases with both R and c on R ≥ c ≥ 0. [Lancaster & Blanchard (1969)](https://ntrs.nasa.gov/citations/19690027552) · [Izzo (2015), parabolic limit](https://www.esa.int/gsp/ACT/doc/MAD/pub/ACT-RPR-MAD-2014-RevisitingLambertProblem.pdf)

The two measurements are on different tangent planes. We convert each relative vector into physical metres in its own local east/north basis at the target range. Each projected separation lower-bounds the corresponding true radius. The unit normal to both sightlines measures one component of the chord unaffected by either unknown line-of-sight depth. We take the smaller absolute difference for the two possible assignments of the nearly equal components. We subtract the two radial error allowances from this chord bound and the corresponding allowances from the radii. The calculation also allows independent fractional distance errors of 0.5%, far larger than needed to describe even a million kilometres of binary depth at 39 AU. This includes first-order depth/perspective uncertainty at this scale. We enlarge the available transfer time for the same distance uncertainty in the light-time correction and add 2 s. The HST visit mean is assigned its mean timestamp; treating its first or final exposure individually gives nearly the same baseline.

For the adopted radii, the HST–first Rubin *i* pair has a minimum endpoint chord of 12,744 km, a minimum sum of endpoint radii of 39,975 km and a maximum elapsed time of 13.7351 days. The *z* exposure on June 27, observed 41 minutes after *i*, gives 1.74 × 10¹⁹ kg with the same assumptions. These two bounds share the HST measurement and instrument systematics. We report the first *i* result; all cross-date pairs are in the CSV.

The HST direction and range come from JPL Horizons solution 21 queried from HST (NAIF −48) near the visit midpoint. Rubin directions/ranges come from the ALeRCE ephemeris columns. Rubin's MJD is the exposure midpoint in TAI; HST's MJD is UTC. Optical paths, nearly equal component brightness, and detector orientation limit this inference. A larger unmodeled bias than the adopted radii weakens or invalidates the numerical bound.

Run `python bound.py` to rebuild `bounds.csv` from the included HST and Rubin measurements. Run `python hst_geometry.py` to request the ephemeris again; the live orbit solution and the service timestamp can change.
