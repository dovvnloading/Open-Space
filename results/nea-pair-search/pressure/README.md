# An outward-acceleration candidate in the 19-year linkage

The joint arc for 2007 VE138 and 2026 VV15 favors an additional outward force varying with inverse-squared heliocentric distance. The test is conditional on the two designations belonging to the same object. This is a measurement from a sparse archive with a 19-year gap. The physical origin remains open.

| Six-state fit or extension | A1 in AU/day² at 1 AU | A2 in AU/day² at 1 AU | All-position radial RMS |
|---|---:|---:|---:|
| Gravity | Fixed zero | Fixed zero | 1.142″ |
| Radial only | +1.362 × 10⁻⁸ | Fixed zero | 0.831″ |
| Transverse only | Fixed zero | +5.85 × 10⁻¹³ | 1.104″ |
| Radial and transverse | (+1.424 ± 0.094) × 10⁻⁸ | (−1.61 ± 1.14) × 10⁻¹³ | 0.828″ |

The uncertainties in the last row are one estimated standard error, clustered by observatory and UTC date. The corresponding formal radial uncertainty treating position errors as independent is 0.110 × 10⁻⁸ AU/day². The same two-parameter fit to 2017 RB15 and 2026 RG39 gives a radial term consistent with zero, (+0.17 ± 0.32) × 10⁻⁸, and a transverse term (−1.11 ± 0.28) × 10⁻¹². On dropping F51 and I52 its transverse uncertainty rises to 1.26 × 10⁻¹². The latter result warrants further astrometry before a physical interpretation.

## Withheld astrometry and weighting

| VE138 observations withheld | Positions | Gravity, radial RMS | Radial + transverse, radial RMS |
|---|---:|---:|---:|
| All J95 and 595 observations | 14 | 1.146″ | 0.707″ |
| Entire final 2026 UTC date (Nov 14) | 2 | 5.882″ | 1.104″ |
| Entire final 2007 UTC date (Nov 13) | 3 | 2.905″ | 0.205″ |
| All 2026 observations after the first 2 | 65 | 5.850″ | Underconstrained |
| All 2026 observations after the first 8 | 59 | 2.506″ | 3.103″ |

These subsets overlap. Each refit starts near the selected joint solution, and both models see the historical apparition. The comparison measures how withheld astrometry behaves within that solution family. The first-two fit is poorly conditioned and depends on the stopping rule. At 20 versus 30 iterations the extended-model holdout RMS moved from 1.75″ to 2.15″. The first-eight test favors gravity. A new observing night or a third apparition would provide a stronger independent check.

The simultaneous A1 remains positive when the site or date is varied. Linearized fits deleting each complete observatory give +1.373 to +1.458 × 10⁻⁸ AU/day². Deleting a whole UTC date gives +1.312 to +1.500 × 10⁻⁸. Equal total weight per station-date gives +1.382 × 10⁻⁸. Using the larger of each measurement's two quoted component sigmas, a 0.5″ floor, and a default of 0.8″ for missing entries gives +1.411 × 10⁻⁸. Omitting the original site's fourteen 2007 positions (703) gives +1.385 × 10⁻⁸. Allowing separate constant east and north offsets for each 2007 star catalog with 1″ Gaussian priors yields +1.442 × 10⁻⁸. The replay also reports a common 2007 coordinate offset, derivative-step sensitivity and all deletion estimates. None of these models represents arbitrary time-dependent astrometric or ephemeris systematics.

## Physical scale and interpretation

One AU/day² is approximately 20.04 m/s². The inferred outward acceleration at 1 AU is about 2.85 × 10⁻⁷ m/s². If solar radiation pressure supplies the entire acceleration, then (C_R A/m=a_1 c/S_0) is about 0.063 m²/kg for solar irradiance (S_0=1361) W/m². Here (C_R) is the momentum-transfer coefficient. Taking (C_R=1), (H=26), geometric albedo 0.14, and a compact spherical shape gives a nominal diameter near 22 m and a bulk density near 1 kg/m³. This illustration strongly disfavors that particular compact-sphere interpretation. An object of unusual shape or low mass, non-radiative forces or systematic errors could imitate the signal. There is no supporting spacecraft identification, spectrum, radar cross-section or measured mass in this packet. The JPL 2026 VV15 nominal close-approach speed of 14.98 km/s on 2026-11-14 and an orbital inclination near 19° would make ordinary Earth-escaping debris an unusual explanation. Check the trajectory of a specific launch before claiming an association. The catalog designations and provisional tracklet LVg0046 were absent from Project Pluto’s public artificial-satellite identification list on 2026-11-15. That list is incomplete. [JPL CAD](https://ssd-api.jpl.nasa.gov/cad.api?des=2026%20VV15&date-min=2026-11-14&date-max=2026-11-15) · [Project Pluto identifications](https://www.projectpluto.com/sat_id2.htm)

## Model and replay

Find_Orb's existing two-component inverse-square force adds ((A_1\hat r + A_2\hat t)/r^2) in heliocentric coordinates, where (r) is in AU, (hat r) points outward and (hat t) follows transverse motion. Find_Orb also accounts for Earth shadow in the illumination factor. The acceleration patch sets the parameter count and streams eight values to the existing numerical engine. Independent finite-difference scales, equal weighting, reweighting and station-date cluster covariance use the frozen source astrometry in the neighboring `orbits` directory. The full Cartesian states, cluster covariance and residuals are in `fits/`.

```sh
bash build_acceleration.sh /tmp/nea-force-build
python replay_acceleration.py /tmp/nea-force-build /tmp/nea-force-results
```

Both commands run from this directory. Build inputs and DE440 are pinned in the neighboring orbit package. [Find_Orb gravitational and nongravitational acceleration implementation](https://github.com/Bill-Gray/find_orb/blob/1de55cca12e7a585c9c7384acbb80538ea2afc3c/runge.cpp) · [Solar radiation-pressure constants](https://github.com/Bill-Gray/find_orb/blob/1de55cca12e7a585c9c7384acbb80538ea2afc3c/constant.h)
