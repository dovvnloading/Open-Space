# Two candidate eclipses of Nereid by Triton

For **31 January 2026**, one Horizons orbit reconstruction **predicts** that Nereid would enter Triton’s penumbra for **27 min 21 s**, with a predicted Earth-reception maximum at **00:34:52 UTC**. A uniformly lit projected disk model gives **67% light loss**, about **1.21 mag**. An archived older orbit from a different group misses Nereid entirely. The disagreement makes these timed predictions directly testable with historical imaging. **This shadow candidate remains unconfirmed: no eclipse was measured, and no exposure covering its window was found in the four public archive indices searched.** Those are bounded archive searches, not proof that no observation exists elsewhere. The alternative orbit with no eclipse prevents treating the computed light loss as an observed event.

| Earth-reception date | Triton and Nereid ephemerides | Peak UTC | Shadow axis to Nereid centre | Modeled loss of disk light |
|---|---|---|---:|---:|
| 2026 Jan 31 | Horizons vectors returned for both moons | 00:34:52 | 1,174 km | 67.2% |
| 2026 Jan 31 | Horizons Triton, JPL NEP105 Nereid | 00:34:08 | 1,178 km | 66.7% |
| 2026 Jan 31 | JPL NEP097 for both | 00:16:57 | 1,500 km | 29.0% |
| 2026 Jan 31 | SAI 2015 Triton, SAI 2011 Nereid | 00:01:01 | 2,272 km | 0% |
| 2024 Feb 07 | Horizons vectors returned for both moons | 17:49:04 | 1,381 km | 42.4% |
| 2024 Feb 07 | Horizons Triton, JPL NEP105 Nereid | 17:48:33 | 1,434 km | 35.9% |
| 2024 Feb 07 | JPL NEP097 for both | 17:38:17 | 2,544 km | 0% |
| 2024 Feb 07 | SAI 2015 Triton, SAI 2011 Nereid | 17:29:44 | 4,001 km | 0% |

### What a telescope could have seen

The 2026 direct-Horizons partial phase spans **00:21:16–00:48:37 UTC**. Triton covered about 68% of the Sun as seen from Nereid’s centre. Earth saw Nereid **3.72 arcsec from Triton** and **14.77 arcsec from Neptune**, at a solar elongation of **49.0°**. The Sun-facing shadow-axis calculation places Neptune itself about **289,000 km** off the relevant ray. A spatially resolved light curve could therefore distinguish shadowing by Triton from blending or obscuration on the Earth–Nereid line of sight. Flux in Triton’s wings and Neptune’s halo would require local PSF subtraction. A second, shallower 2024 opportunity has fewer agreeing orbit solutions.

The 2026 first-contact threshold is about **2,042 km** from the axis, using a 1,352.6 km Triton, 170 km Nereid, and a 519 km projected solar radius. The independent SAI paths miss this threshold by **230 km**. The 2024 first-contact threshold is about **2,003 km**. The model spread reflects the choice of orbital solution. These fits use overlapping astrometry. There is no calibrated event probability, and the second-resolution timestamps describe numerical predictions within each fit. Nereid’s adopted radius has a quoted 25 km uncertainty. The light loss assumes uniform brightness across Nereid’s projected disk and across the solar photosphere. Terrain, phase, solar limb darkening and scattered light change the inferred flux. [JPL satellite physical parameters](https://ssd.jpl.nasa.gov/sats/phys_par/)

### Selection and coverage

A discovery screen searched 1900–2160 TDB for projected close approaches of Triton and Nereid along both the Sun and Earth sight lines. A one-day Triton and four-day Nereid/JPL observer grid was Hermite-interpolated with returned velocities and sampled every 90 minutes. Local minima with sampled separation below 35,000 km were refined. This report selects the two recent solar-direction candidates in the 2000–2040 window. The coarse sweep supplies seeds. It is not a completeness or event-rate estimate. `screen.py` reproduces the 2000–2040 discovery subset. The dates of Earth reception differ from local shadow passage by over four hours.

`archive.py` checks image exposure indices within one hour on either side of each Horizons midtime. `archives.json` retains the exact queries, returned candidates, service failures and limits. ESO and CADC returned zero overlapping footprints on both nights. NOIRLab returned zero pointings within four degrees. MAST returned 19 distinct 2024 footprints in a four-degree boresight box and two in 2026. Every footprint missed. The closest 2024 TESS catalogue aperture is centred 0.164° away with radius 0.00139°, and both 2026 HST footprints lie near RA 239.7°, Dec +35.45° despite zero-valued catalogue boresights. The MAST boresight prefilter can omit off-centre targets in wide-field imaging. An index miss leaves independent archives, unavailable exposure footprints and private imaging unconstrained. The **23 November 2024** JWST NIRSpec observation of Nereid occurred months after the first eclipse window. Its spectra do not sample either event. [JWST observation and analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC13189092/)

### Reproduction and provenance

Run from this directory with Python 3, NumPy and SciPy.

```sh
python replay.py > reproduced.json
python screen.py
python acquire.py refreshed-states.csv
python archive.py > refreshed-coverage.json
```

`states.csv.gz` freezes the local geometric ICRF position and velocity series at 30-minute spacing. `acquire.py` retrieves JPL Horizons inertial states and evaluates local Chebyshev polynomial records read from the public SPK files. Native centres are Neptune or its system barycentre, as specified by each SPK segment. Replay keeps the reception event on Nereid, solves light travel from Triton and then to Earth, and computes partial-solar-disk overlap. The solar light-time uses local linear solar barycentric motion. It computes projected Nereid disk flux via 32 radial quadrature nodes and 192 azimuths. Geometric midtime is defined by minimum distance to the shadow axis. UTC conversion uses TT minus UTC of 69.184 s for the two actual dates, neglecting the millisecond TDB–TT term. Geometry and adopted spherical radii dominate timing precision.

The Horizons 802 response in the search labeled its source `nep098_merged`. JPL separately lists NEP105 as its designated Nereid solution. Both are evaluated explicitly here. Direct Horizons and NEP105 share an ephemeris provider. NEP097 is a prior JPL solution. The SAI Nereid and Triton series were published in 2011 and 2015. The 2021 FORCES-8-MAIN-2020 model is another relevant comparison. Its listed download returned HTTP 503 during this investigation and has not been used. [JPL Horizons](https://ssd.jpl.nasa.gov/api/horizons.api) · [JPL ephemeris inventory](https://ssd.jpl.nasa.gov/sats/ephem/sep.html) · [JPL SPK files](https://ssd.jpl.nasa.gov/ftp/eph/satellites/bsp/) · [SAI kernels](https://ftp.imcce.fr/pub/ephem/satel/SAI/NEPTUNE/) · [FORCES paper and data](https://doi.org/10.12149/101078)

A literature check on 15 November 2026 searched both moon names with “mutual eclipse,” “occultation,” and each event date. The 2017 exo-syzygy study includes the geometric possibility of shadowing between these moons. Mallama’s 2007 predictions concern Neptune casting its own shadow on Nereid. The search did not locate a dated report of either Triton-shadow event. Publication priority and actual occurrence remain unverified. [Veras & Breedt 2017](https://academic.oup.com/mnras/article/468/3/2672/3070428) · [Mallama 2007](https://doi.org/10.1016/j.icarus.2007.02.001)
