# A predicted encounter between JUICE and (41013) 1999 UH22

ESA’s published cruise trajectory brings JUICE **49,751 km from (41013) 1999 UH22 on 25 November 2027**, at **6.324 km/s**. Closest approach is at 19:24:06 TDB, approximately **19:23 UTC** under the current leap-second convention. The two bodies spend **7 h 37 min inside 100,000 km**. The geometry gives a bright asteroid against the sky seen from the spacecraft. Estimating the number of camera pixels across its disk still requires an assumption about the asteroid’s reflectivity. [Numerical results](measurement.json) · [ESA trajectory archive](https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/spk/)

## Three asteroid orbits and nine versions of the spacecraft trajectory

JPL orbit 34 for the asteroid uses **850 positions from November 1995 to January 2026**. We separately requested heliocentric mean J2000 equatorial states from the Paris Observatory Miriade service using the Lowell ASTORB and Minor Planet Center MPCORB orbit solutions and the INPOP19A planetary ephemeris. The three reductions share much of their underlying telescope astrometry. Combining each asteroid orbit with ESA’s same October 2026 spacecraft solution gives the following predictions. [JPL orbit and covariance](https://ssd-api.jpl.nasa.gov/sbdb.api?sstr=41013&cov=mat) · [Miriade definitions](https://ssp.imcce.fr/webservices/miriade/api/ephemcc/)

| Asteroid orbit | Closest approach, 25 November 2027 TDB | Distance |
|---|---:|---:|
| JPL Horizons | 19:24:06 | 49,751 km |
| Lowell ASTORB via Miriade | 19:24:08 | 49,787 km |
| MPCORB via Miriade | 19:23:59 | 49,748 km |

Holding the JPL asteroid orbit fixed, **eight post-launch ESA trajectory revisions** give 49,603–50,795 km. The revisions are 000100, 000110 and 000112–000117. The last was informed by tracking through 22 October 2026. We also recovered the encounter at **46,520 km** on the pre-launch CReMA 5.0b23.1 trajectory used in earlier flyby work. Direct evaluation of ESA kernel 000117 with CSPICE and separate heliocentric Horizons vectors agrees with the spacecraft-centred Horizons distance to better than 1 km. Each revision represents a different navigation estimate. Their spread measures sensitivity to those revisions. Future manoeuvres can change the encounter.

## What the camera could measure

At closest approach the JANUS camera’s 15-microradian sampling gives **0.746 km per pixel**. JPL reports absolute magnitude H = 16.21 and supplies no measured diameter. Using the standard magnitude–diameter relation yields the conditional sizes below. Angular sampling alone cannot predict resolved surface detail. The point-spread function, tracking, illumination, spacecraft attitude and available observing time also matter. [JANUS camera measurements](https://doi.org/10.5194/angeo-44-1243-2026) · [JPL magnitude–diameter relation](https://cneos.jpl.nasa.gov/tools/ast_size_est.html)

| Assumed geometric albedo | Inferred diameter | JANUS pixels across the full diameter |
|---:|---:|---:|
| 0.04 | 3.81 km | 5.1 |
| 0.14 | 2.03 km | 2.7 |
| 0.40 | 1.20 km | 1.6 |

At the brightest ten-minute sample, **18:40 UTC**, Horizons predicts apparent V = 1.30, solar phase angle 5.7° and solar elongation 174.3° from JUICE. This brightness uses the catalogue H–G model, including its assumed phase law. At 19:20 UTC the asteroid moves **26.2 arcseconds per second** in the inertial sky. Without target tracking, one JANUS pixel would be crossed in 0.12 s. At the same epoch Horizons gives a formal asteroid range uncertainty of 42 km (3σ) and a 3σ sky-plane ellipse semimajor axis of 470 arcsec. Those formal errors hold the spacecraft trajectory fixed. The full observer table includes sky coordinates and rates every ten minutes. No mission observation commitment was found. [Frozen observer table with units](data/obs-fine.json.gz)

## How the target was found

The initial search considered JUICE, Psyche, Europa Clipper and Lucy cruise tracks. Four known Lucy Trojan flybys were recovered. This directory reproduces the JUICE search from 15 November 2026 through 15 July 2029. The JPL numbered and unnumbered asteroid tables supplied **1,467,659** entries with common epoch MJD 61200 TDB and semimajor axes strictly between 1.3 and 5.5 au. For each orbit we propagated osculating two-body elements, divided the spacecraft trajectory into four-day intervals, and minimized the linearized relative distance inside each interval. We used daily heliocentric spacecraft vectors in the J2000 ecliptic frame. The code records the plane and radial pruning cuts. [JPL element table definitions](https://ssd.jpl.nasa.gov/sb/elem_tables.html) · [Full selection code](screen.py)

The 0.012 au cut retained **85 JUICE candidates**. We requested spacecraft-centred, geometric, three-hour Horizons vectors for all 85 over thirteen days around their individual predicted times, then minimized distance using cubic Hermite interpolation of position and velocity. Every refined minimum was inside its window. **(41013) was the closest of these 85**. The runner-up, (274165) 2008 FU128, passed at 125,721 km. The two-body search initially placed (41013) about 365,000 km away, illustrating why a tight cut at that stage could miss a real encounter. The coarse dynamics, element epoch, single retained window per asteroid and semimajor-axis range limit the survey’s completeness. [All 85 objects, including source-response hashes](data/encounters.csv)

## Earlier flyby searches and the next test

Agostini et al. (2022) searched about 140,000 asteroids and reported 71 potential JUICE targets, highlighting (1650) Heckmann and (223) Rosa. Their companion study describes using asteroids with known sizes. The published summary identifies Rosa as the sole object exceeding one JANUS pixel in that analysis. Our search included asteroids with unknown diameters. On 14 November 2026 we checked the authors’ abstracts, the companion open-access paper, presentation, ESA mission material, and scholarly web searches for “41013 JUICE”, “1999 UH22 JUICE” and “25 November 2027 JUICE asteroid”. We found no explicit published identification of this encounter. We checked all 71 target IDs in Table 2.2 of Agostini’s open doctoral thesis, which presents the earlier survey. It omits 2041013 (the numbered-asteroid SPICE ID convention used in that table). As a positive control, it includes (52073) at 324,000 km. Our current trajectory gives 325,220 km. This verifies an addition to that specific published list. Worldwide priority remains unresolved. The pre-launch reconstruction shows that the encounter was already present in the older trajectory. [Earlier survey](https://doi.org/10.1016/j.pss.2022.105476) · [Companion study](https://doi.org/10.1051/0004-6361/202142600) · [Thesis and complete 71-target table, p. 20](https://hdl.handle.net/11577/3472986) · [Author’s trajectory specification](https://cisas.unipd.it/sites/cisas.unipd.it/files/35%20AGOSTINI.pdf)

The next ESA spacecraft solution and updated asteroid astrometry provide a direct test of the predicted range and time. A thermal diameter or occultation chord would turn the angular-sampling estimates into a measured target size. A mission attitude and instrument timeline would be needed to establish observing feasibility.

## Reproduction

Python 3, NumPy and SciPy suffice for the offline replay. It verifies the fine-grid source hashes, interpolates heliocentric positions and velocities, minimizes geometric distance, and reads the apparent conditions from Horizons’ separate light-time-corrected observer product. Vector epochs are TDB. Miriade supplies TT, whose periodic difference from TDB is below 2 ms here. Apparent brightness and sky position use the UTC observer table. Source URLs, byte counts and SHA-256 hashes are in the [manifest](data/sources.json).

```sh
python replay.py > replay.json
diff -u measurement.json replay.json
```

For the population replay, supply the two JPL gzip archives named in the manifest. Live catalogue tables change. The stored candidate list and frozen ephemerides allow review after those catalogue versions roll over. The refinement command caches its service responses and issues one request at a time.

```sh
OPENBLAS_NUM_THREADS=1 python screen.py /path/to/element-tables candidates.csv
python refine.py candidates.csv /path/to/horizons-cache refined.csv
```

The ESA state table contains half-hour samples in heliocentric J2000 equatorial coordinates. To regenerate those rows from the original SPKs, obtain the nine exact kernels in the manifest and run the command below with SpiceyPy installed. Inputs are attributed to ESA Flight Dynamics and the ESA SPICE Service, NASA/JPL Solar System Dynamics, IMCCE/LTE Miriade, Lowell Observatory and the Minor Planet Center.

```sh
python sample_esa.py /path/to/kernels > esa-states.csv
```
