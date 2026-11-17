## Spacecraft astrometry and blind masks

The October 2025 and October 2026 JWST 8975 cubes use the observatory at 500@-170. The HST 16995 check uses 500@-48. The JPL census supplies Uranus IDs 701–727, 75051 and 75052, and Neptune IDs 801–814, 85051 and 85052. The observer table requests ICRF astrometric RA/Dec, angular diameter and heliocentric and observer distances. The short visit groups use spacings of at most six minutes and cubic interpolation to uniformly sampled shutter times. Horizons responses must match the target identity and every requested epoch. The receipt files give a SHA-256 digest of every consumed response or kernel excerpt. The observables tables preserve all numerical samples. [JPL ephemerides](https://ssd.jpl.nasa.gov/sats/ephem/) · [Horizons observer manual](https://ssd.jpl.nasa.gov/horizons/manual.html)

**S/2025 U1 requires a special case.** During this run, numeric Horizons command 75052 selected asteroid *75052 (1999 UM50)*. JPL's published URA184 kernel contains the correct J2000 segments from the Uranus-system barycenter to S/2025 U1 and Uranus. The code differences those two states at the emission time inferred from the spacecraft-specific Horizons distance, with a differential light-time correction for the moon. For all 24 HST Uranus exposures, the identical calculation for Cordelia agrees with direct Horizons to a median 0.025 and maximum 0.044 milliarcseconds. This checks geometry. Uncertainty in S/2025 U1's predicted orbit away from its discovery epoch requires a separate allowance. For this reason the HST sensitivity report excludes a wider annulus. [URA184 kernel](https://ssd.jpl.nasa.gov/ftp/eph/satellites/bsp/ura184.bsp) · [JPL satellite elements](https://ssd.jpl.nasa.gov/sats/elem/)

JWST calibrated SCI frames are transformed using their native SIP distortion and individual planet-limb fits. The grid has east to the left and array row index increasing toward celestial north. The input surface brightness in MJy/sr is multiplied by the output pixel solid angle and by a million to obtain Jy/pixel. The ERR-plane variance uses the squared bilinear weights. Correlations between resampled pixels remain and motivate the blank-position tests and in-image injections. The permissive calibrated mask accepts DQ values 0 and 4 (jump detected and calibrated) alone. The strict comparison accepts only 0. All variants exclude nonfinite or nonpositive errors, clipped images, planet disk, projected ring envelopes and the sampled trail of every recognized moon. The 8975 search adds five output pixels around the disk, ring envelope and each moon's angular radius. No catalog covariance is implied by that padding. Optical ghosts or PSF wings outside those masks can still generate scores.

The Uranus ring masks use 26,800–51,650, 65,400–69,200 and 89,200–106,200 km. Neptune uses 41,000–43,000, 53,150–57,300, 61,850–62,050 and 62,850–63,020 km. The entire Adams ring is masked together with its arcs. Ring and disk masks reduce the observable fraction. The radial envelope and five-pixel padding are hard exclusions; extended light just outside them is checked separately. [Uranian rings](https://pds-rings.seti.org/uranus/uranus_rings_table.html) · [Neptunian rings](https://pds-rings.seti.org/neptune/neptune_rings_table.html)

## Orbit conventions

The circular search uses a central gravitational field with J2 and J4 and the angular-momentum pole of the regular moons. Uranus's angular-momentum pole is opposite its IAU north pole. The phase origin is the ascending node of the equator on the ICRF equator. Geometry changes with observer direction and range. The epoch is the median TDB emission time of the sampled exposures. A 2-D displacement at each shutter sample is projected onto that observation's tangent plane. Optional initial osculating Kepler elements are integrated with the J2 and J4 force using DOP853. The circular survey samples only coplanar circular orbits. The optional integrations omit inter-moon perturbations, solar perturbations and high-degree gravity; long-baseline orbit fitting requires astrometry, alternate solutions and uncertainty propagation. No fitted orbit or covariance exists for a new body here. [NAIF pole model](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00011.tpc) · [Gravity fields](https://doi.org/10.1093/mnras/stac628)

## Reproduce registered JWST inputs

Run from this result directory with Python, NumPy, SciPy, Astropy and jplephem. In these commands `cache` holds archive products and temporary files. Create it before writing outputs. The source scripts never overwrite the archive images. The center tables identify the eight contributing calibrated exposures per band.

~~~bash
mkdir -p cache/horizons cache/ranges
python moon_spk.py 2460954.42 2460954.75 cache/ura184-2025oct.bsp cache/ranges
python ephemeris.py Uranus JWST 60954.090 60954.205 cache/horizons cache/f162m-eph.npz --ura184 cache/ura184-2025oct.bsp
python prepare.py cache/f162m-eph.npz cache/f162m.npz cache/jw0897500[13]*_cal.fits --centers jwst8975/f162m_centers.csv --scale .031 --size 421 --samples 33 --psf-sigma .90 --allow-corrected-jumps
python moon_spk.py 2461321.32 2461321.70 cache/ura184-2026oct.bsp cache/ranges
python ephemeris.py Uranus JWST 61321.008 61321.172 cache/horizons cache/f140m-eph.npz --ura184 cache/ura184-2026oct.bsp
python prepare.py cache/f140m-eph.npz cache/f140m.npz cache/jw0897500[57]*_cal.fits --centers jwst8975/f140m_centers.csv --scale .031 --size 421 --samples 65 --psf-sigma .85 --allow-corrected-jumps
python moon_spk.py 2459892.70 2459894.40 cache/ura184-2022hst.bsp cache/ranges --check-body 706
python hst16995_spk.py hst16995 cache/ura184-2022hst.bsp
~~~

The shutters for each 8975 CAL span four integrations with intervening resets. Templates approximate uniform exposure weighting across the full shutter span. The source and injections both use normalized Gaussian cores sampled along the trial orbits. The well-known-moon controls and integration-level 6379 check constrain interpretation of this approximation. The search never counts integrations or detector dithers from the same visit as independent visits.

## GO 8975 wider F162M search and glare veto

The wider F162M survey evaluated 66,600 trials from 36,000 to 128,000 km in 500 km × 1° steps, with 33 shutter samples per exposure. Only 8,807 met the two-visit support criterion. Five neighboring phase solutions at 55,000 km passed the first score thresholds; the inverted images produced zero passing trials. The highest-ranked hypothesis (104°) had 5.65 and 6.57 in the individual visits and 8.46 combined. Its modeled track crosses eight detector frames (seven contribute to its score) and spans 14.38 arcseconds in sidereal coordinates. At each exposure, the closest recognized moon along this track is at least 46 registered pixels away. These are properties of the predicted track, not measured point-source astrometry.

Each exposure shows this track close to the outer edge of the ring glare. Sampling ten adjacent trails at the same orbital radius, ±6°, 8°, 10°, 12° and 15° away in phase, measures the diffuse flux at comparable positions. Subtracting their median within each exposure reduces the scores to 1.77 and 1.76 (2.47 pooled) for the 104° track. All five phase solutions fail the second-visit requirement under this control. The calculation retains the original weights without adding the background estimate's variance, so the corrected scores are optimistic. The numerical results and all frame coordinates are in [f162m_broad](jwst8975/f162m_broad/). The radial grid has no standalone completeness magnitude. Consult the separate population-matched injection tables for any depth statement.

~~~bash
python broad_scan.py cache/f162m.npz cache/f162m-broad --grid 36000 128000 500 1.0
python vet_tracks.py cache/f162m.npz cache jwst8975/f162m_centers.csv cache/f162m-broad/passing.csv cache/f162m-broad
~~~

At the predicted S/2025 U1 location in this October 2025 dataset, unmasking that moon alone gives 5.74 in the first visit (three supported exposures) and 1.33 in the second (two exposures), or 5.27 combined. The supplied ephemeris and photometry fail to recover this known object at the two-visit cuts in 8975. These numbers do not give an upper limit on the moon. To reproduce, add `--control-body 75052` to the F162M preparation command with a separate output file, then run `known_control.py CONTROL.npz 75052 OUTPUT.json`. [Forced-control output](jwst8975/s2025u1_f162m_control.json)
