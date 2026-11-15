# Independent integration and direction-dependent catalog corrections

ASSIST and REBOUND integrate an independent orbit using the same DE440 ephemeris and MPC positions as the Find_Orb analysis. Astropy computes the UTC-to-TDB conversion and the station location in celestial coordinates using its bundled Earth-orientation tables. The fitter projects the modeled light-time-corrected direction onto the observed east and north axes, then solves for six initial-state components and, optionally, an outward radial term proportional to inverse squared heliocentric distance. Each fit starts from the six-parameter gravity solution with zero radial acceleration. All 110 positions receive equal coordinate weight. The independent implementation checks software and reference-frame choices; both reductions share observational systematics.

| Astrometry | Gravity radial RMS | Gravity + radial acceleration RMS | Acceleration at 1 au (au/day²) |
|---|---:|---:|---:|
| Submitted positions | 1.1420″ | 0.8307″ | +1.36207 × 10⁻⁸ |
| With EFCC18 corrections | 1.1810″ | 0.8404″ | +1.44222 × 10⁻⁸ |

The second run subtracts the published tile-specific catalog position and proper-motion biases. The table applies to all 43 positions from 2007 (20 UCAC2, 13 USNO-B1, 10 USNO-A2), plus 10 positions from 2026 reduced using UCAC4. The other 57 positions use Gaia or ATLAS2. Gaia1 is among the catalogs JPL leaves uncorrected. ATLAS2 has no entry in the table. The output records each applied correction, its HEALPix index, and the observation identifier. Residual RMS values are directly comparable between force models within each row. [Eggl et al., Icarus](https://doi.org/10.1016/j.icarus.2019.113596) · [JPL debias table](https://ssd.jpl.nasa.gov/ftp/ssd/debias/debias_2018.tgz)

## Replay

Run in this directory after building the gravitational fit. In a fresh Python 3.12 environment, install the dependencies below. ASSIST may need a C compiler at installation. The archived run used NumPy 2.5.4 and SciPy 1.17.0.

```sh
CC=gcc python -m pip install 'rebound==4.6.0' 'assist==1.2.3' 'astropy==8.0.2' 'astropy-healpix==2.0.1' 'astropy-iers-data==0.2026.11.9.0.35.9' numpy scipy
curl -fL https://ssd.jpl.nasa.gov/ftp/ssd/debias/debias_2018.tgz -o /tmp/efcc18.tgz
python independent.py /tmp/nea-find-orb-build/de440.bin ../orbits /tmp/independent.json
python independent.py /tmp/nea-find-orb-build/de440.bin ../orbits /tmp/corrected.json --bias-archive /tmp/efcc18.tgz
```

The source date for the eight optical catalog entries is 2026 November 15. Inputs to the orbit fit are optical directions, observing stations and observation times. No differential aberration or gravitational light-deflection term is added to the reduced astrometric positions. The light-time calculation uses local linear motion over the few seconds of travel time. ASSIST includes planetary gravity, solar relativistic effects and the Earth and Sun harmonics in the chosen implementation. Both integrations omit massive-asteroid perturbers for this check. [ASSIST force model](https://assist.readthedocs.io/en/latest/)
