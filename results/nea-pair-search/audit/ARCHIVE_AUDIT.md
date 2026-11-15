# Designation, prior-publication, and resolved-track checks

On 2026 November 15 the Minor Planet Center Designation Identifier returned each of the eight names below as its own primary designation, with no secondary designations. ESA's NEO Coordination Centre also returned a named orbit element set for every one. The two agencies consume the same underlying optical reports. Separate entries establish catalog status only. These checks cannot exclude a private or pending identification.

| Designation | Discovery circular | Discovery station | Observations in current MPC archive |
|---|---|---|---:|
| 2026 SW11 | [2026-S141](https://www.minorplanetcenter.net/mpec/K26/K26SE1.html) | 703 | 20 |
| 2026 SA12 | [2026-S143](https://www.minorplanetcenter.net/mpec/K26/K26SE3.html) | I41 | 28 |
| 2026 VO | [2026-V25](https://www.minorplanetcenter.net/mpec/K26/K26V25.html) | V00 | 21 |
| 2026 VR1 | [2026-V62](https://www.minorplanetcenter.net/mpec/K26/K26V62.html) | G96 | 19 |
| 2017 RB15 | [2017-S26](https://www.minorplanetcenter.net/mpec/K17/K17S26.html) | F51 | 19 |
| 2026 RG39 | [2026-S03](https://www.minorplanetcenter.net/mpec/K26/K26S03.html) | F52 | 43 |
| 2007 VE138 | [2007-V86](https://www.minorplanetcenter.net/mpec/K07/K07V86.html) | 703 | 43 |
| 2026 VV15 | [2026-V305](https://www.minorplanetcenter.net/mpec/K26/K26VU5.html) | O18 | 67 |

Inspection of the sixteen 2026 circulars referenced by the MPC astrometry showed no explicit pairing of those designations. The September 23 Daily Orbit Update lists SW11 and SA12 separately. The November 6 Daily Orbit Update lists VO and VR1 separately. Search engines found co-listings of each same-apparition pair but no earlier explicit identification in the checked results. [September update](https://www.minorplanetcenter.net/mpec/K26/K26SF1.html) · [November update](https://www.minorplanetcenter.net/mpec/K26/K26VD4.html)

## Same-site measurements

At station H36 on September 23, the two SA12 positions precede the two SW11 positions. A quadratic fit to RA and Dec across the 44.7-minute series gives a maximum tangent-plane residual of 0.375 arcsec. At station 291 on November 6, three VR1 points precede three VO points. A single straight angular track through the 46.2-minute series misses no measurement by more than 0.738 arcsec. A quadratic reduces this to 0.353 arcsec. Neither pair has a same-station exposure pair within 180 seconds. H36 offers only one residual degree of freedom per coordinate for its quadratic fit.

The control 2018 NC15 and 2026 LD1 gave eight simultaneous angular positions separated by 1.250–1.270 degrees from station G96, four on June 5 and four on June 21. These prove that two sources were measured. The first pair of positions was reported at 2026-06-05T05:32:18.505Z, with separation 4,500.858 arcsec. The original June 5 measurements appear in separate June 5 and June 10 circulars, and both bodies are retained in the June 21 daily update. These angular observations alone do not establish a common physical origin. [June 5 circular](https://www.minorplanetcenter.net/mpec/K26/K26L10.html) · [June 10 discovery circular](https://www.minorplanetcenter.net/mpec/K26/K26L56.html) · [June 21 update](https://www.minorplanetcenter.net/mpec/K26/K26M59.html)

The script `verify_archives.py` requires Python 3 and internet access. Run `python verify_archives.py` to recheck the official aliases, ESA orbit names, observation counts and same-station same-timestamp separations. OBS80's lower coordinate precision changes the last decimal relative to ADES. [MPC Designation Identifier API](https://docs.minorplanetcenter.net/mpc-ops-docs/apis/query-identifier/) · [MPC Observations API](https://docs.minorplanetcenter.net/mpc-ops-docs/apis/get-obs/) · [ESA orbit download](https://neo.ssa.esa.int/computer-access)
