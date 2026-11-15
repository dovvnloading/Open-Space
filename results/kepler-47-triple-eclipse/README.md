# Kepler-47 d crossing both stars during the primary eclipse of 2029 July 4

**A testable five-body forecast.** At the minimum projected separation, Kepler-47 d lies in front of star B, and star B lies in front of A. For **98.5 minutes** the entire planetary disk falls inside both stellar disks. The planet therefore blocks the faint secondary while passing over a patch of the bright primary that B already hides. Relative to adding the separate transits of B and d over A, the modeled light curve recovers approximately **0.53% of the system light** near closest approach. The system will still be deep in primary eclipse.

| Nominal phase | Barycentric JD (TDB) | Approximate Earth-center UTC, July 4 |
|---|---:|---:|
| First common area of all three disks | 2462322.157977 | 15:44 |
| Entire planet first inside both stars | 2462322.186667 | 16:25 |
| Smallest d–B separation | 2462322.220862 | 17:14 |
| Entire planet last inside both stars | 2462322.255096 | 18:03 |
| Last common area | 2462322.284046 | 18:45 |

Time precision here identifies the numerical solution. The uncertainty of the real orbit is larger. At midpoint the planet center is 0.457 A radii from the primary center and 0.488 times the sum of the B and d radii from B. The target lies **109.46° from the Sun**. In northern East Asia, the nominal event occurs at night with the target high in the sky. The geocentric conversion subtracts 169.21 light-seconds projected along the target direction and approximately 69.184 seconds for TDB–UTC, assuming the present leap-second offset persists. Use barycentric times for fitting.

## What supports the forecast

- We reintegrated Orosz et al.'s five barycentric states and masses with pairwise Newtonian gravity and the dominant binary 1PN correction. Light-time-adjusted sky positions give disk contacts and foreground ordering. Only the nominal state is available in the inputs. A separate fixed-step integration was used while scouting the geometry.
- All **16** released d transit forecasts from 2019–2026 in their Table 10 reproduce within **1.96 minutes**, with median absolute discrepancy **0.62 minute**. This checks the calculation against the authors' implementation. For six Kepler observations from the same fit, median absolute timing residual is **4.21 minutes**. The worst is **53.26 minutes** at BJD 2455418, a grazing event with one cadence missing. It differs from the source's own fitted prediction by about the same amount. These six observations are in-sample.
- A separate test used genuine, later **TESS SPOC** light curves of TIC 271548206 from ten sectors between 2019 and 2024. Of 34 nominal primary eclipses lying at least 0.4 day from the edges of their sector, 18 passed coverage and half-depth-crossing cuts. Their median absolute observed-minus-calculated timing difference is **1.55 minutes**, with a maximum of **6.11 minutes**. Starspots, changing eclipse shape, and noisy light curves contribute to the timing estimator. None of the Table 10 d transit midpoints falls within these ten sector spans. TESS therefore tests the binary phase independently, leaving d's post-Kepler phase unmeasured here.
- The forecast uses the source's Kepler-band quadratic limb darkening for A, the measured radii, a luminous B at approximately 0.5% of A's flux, a uniform B surface and an opaque planet. At closest approach the recovered primary contribution is **5,436 ppm**. The planet still occults approximately **181 ppm** of B light. The resulting **5,255 ppm** difference refers to binary flux outside eclipse and to a sum of the two isolated primary-star transits. Passband, dilution, radius and surface-brightness uncertainties affect this value. Median formal two-minute uncertainties in the TESS light curves span roughly 2.5–3.9%, so a single TESS-like exposure would not resolve this signal.

## Robustness and falsification

A local stress test shifts d along its trajectory while keeping the binary fixed and allows a displacement perpendicular to the binary's projected motion, in primary radii. Positive phase shift evaluates the planet at a later trajectory time at each fixed binary time. This is a sensitivity calculation, not a draw from the fitted posterior.

| Displacement | Planet-phase shift retaining positive-area overlap while d lies entirely on A | Shift retaining complete d containment in A and B |
|---:|---:|---:|
| −0.10 A radii | −2.15 to +3.13 h | −1.76 to +2.73 h |
| 0 | −2.00 to +2.95 h | −1.54 to +2.46 h |
| +0.10 A radii | −1.72 to +2.63 h | none in ±12 h |

The first interval is a sufficient geometric criterion for some restored light. Near its edges the common area can be arbitrarily small. The published December 2026 forecast for d has a **0.95 h** primary-transit timing error and **0.048 A radii** impact-parameter error. Those are uncertainties for a different encounter, with a different relative planet/star speed. The 2019 supplement contains no joint state samples from which to calculate a calibrated probability in 2029. A new d transit time and chord would directly constrain the uncertain trajectory. Continuous observations across the binary eclipse and its shoulders can test the predicted disappearance and return of d's deficit against neighboring binary eclipses. A shifted d transit can miss B entirely.

The relevant prior art includes the known 2009 grazing transit of d blended with a Kepler-47 primary eclipse, ordinary Kepler-47 transit forecasts through December 2026, and syzygies observed in Kepler-1647. The 2009 event does not place d on B in the nominal solution. The specific 2029 triple-disk containment and its photometric signature were not identified in the sources checked. The NASA Exoplanet Archive lists Orosz et al. 2019 as its sole parameter solution for d. No claim of a new planet or an observed 2029 event is made.

## Reproduction

From this directory, using Python, NumPy and SciPy.

```sh
python reproduce.py > fresh-evidence.json
python tess_check.py --cache /tmp/kepler47-tess > fresh-tess-evidence.json
```

The first command reads the frozen initial state and the published d transit table. It recomputes all contacts, geometry, flux contrasts, cross-checks, and stress tests. The second retrieves ten individual two-minute SPOC light curves, checks their SHA-256 hashes, and repeats the independent eclipse timing measurements. `tess-products.csv` preserves exact MAST product identifiers. `sky-geometry.json` contains the Horizon vectors, directions and exact queries. The coordinate and radius columns use AU, days, and solar masses as indicated by their headers. All event times throughout the numerical evidence are barycentric.

## Sources

- Orosz et al. (2019), Tables 2, 5, 6, 7 and 10, including initial Cartesian states at BJD 2454965.0. [Article](https://doi.org/10.3847/1538-3881/ab0ca0) · [Open full text](https://arxiv.org/pdf/1904.07255) · [Machine-readable release](https://cdsarc.cds.unistra.fr/viz-bin/ReadMe/J/AJ/157/174)
- TESS SPOC mission light curves. [MAST archive](https://archive.stsci.edu/missions-and-data/tess)
- Adopted system coordinates and fit inventory. [NASA Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu/overview/Kepler-47%20d)
- Solar and terrestrial barycentric position vectors. [JPL Horizons](https://ssd.jpl.nasa.gov/horizons/)
- Kostov et al. (2016), Kepler-1647b during a stellar eclipse. [Article](https://doi.org/10.3847/0004-637X/827/1/86)
