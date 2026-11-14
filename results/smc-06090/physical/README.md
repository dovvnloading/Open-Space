# A geometrical light-curve test for OGLE-SMC-ECL-06090

This directory fits the flux removed when two spherical stars cover one another on an eccentric orbit whose line of apsides rotates uniformly. The original I-band data are 637 OGLE-III observations and 711 OGLE-IV observations over 22.6 years. The V-band check uses another 48 and 47 measurements. The input files retain their published values and are checked by SHA-256 on every run. [OGLE-III photometry](https://www.astrouw.edu.pl/ogle/ogle3/OIII-CVS/smc/ecl/phot/) · [OGLE-IV photometry](https://www.astrouw.edu.pl/ogle/ogle4/OCVS/smc/ecl/phot/)

## Orbit and eclipse models

For each individual exposure, solve Kepler’s equation, rotate the position by the evolving argument of periastron and project it at the fitted inclination. The sign of the line-of-sight coordinate determines which disk is in front. Compute the exact overlap area of two uniform disks, then perform Gauss–Legendre integration over uniform subdisks for linear limb darkening. Both stars use the same fixed band-specific coefficient. Between 32 and 64 radial quadrature nodes, the maximum per-observation difference across the dilution tests was 0.011 millimag. Masses, spectroscopic flux ratios, radiative atmosphere models, Roche distortion and detailed tidal dissipation are outside this fit.

The 13 fitted coefficients in array order are the phase zero point and drift per decade, precession period (years), eccentricity, argument of periapse at HJD 2456000 (radians), sum of radii divided by semimajor axis, radius ratio, cosine of inclination, central surface brightness ratio, two baseline magnitudes, and cosine semiamplitudes in magnitude at twice and once the sidereal frequency. Coordinates use mean longitude, so that the precession term is consistently subtracted from mean anomaly. Phase uses P=1.1696719 days and HJD 2460000.0369. Uniform extra light is specified as a fraction of the total uneclipsed flux. Plus and minus designate the sign of the precession rate relative to the mean orbital motion. The sense is observationally unresolved.

The six OGLE exclusion windows use boundaries HJD 2453500, 2455000, 2457000, 2458300 and 2459700. Refit the full geometry to the other five blocks with fixed initial reference values chosen from existing knowledge of the system. Phase, drift, centers and shapes are all refitted. Retrospective testing conditional on the selected star does not provide a population significance. The score is twice the soft-L1 loss on magnitude residuals divided by the observational error and a 0.012 mag floor added in quadrature. Its scale is 1.8. Scores are comparable within each experiment and are not chi-square values.

| Geometrical model | Full-data robust loss | Pooled withheld robust loss |
| --- | ---: | ---: |
| Forward precession | 1419.5 | 1485.6 |
| Reverse precession | 1448.6 | 1532.0 |
| Static eccentric orbit | 2685.0 | 3478.2 |

At zero dilution and linear limb coefficient 0.4, the forward fit yields U=28.45 yr, e=0.0669, inclination 83.89 degrees, sum of fractional radii 0.384, and radius ratio 0.786. It predicts five of six OGLE blocks better than reverse precession. For 120 paired bootstrap resamples of observing nights within each of six blocks, the 2.5th and 97.5th percentiles of reverse-minus-forward loss are −10.1 and +60.0. Forward won 109 out of 120 paired draws. This is a conditional robustness check and does not establish the direction. The analogous ranges for U and e were 27.79–29.28 yr and 0.0639–0.0696. The range is conditional on this geometric and noise model, the 22.6-year observing window and the reference orbital timing solution.

A possible systematic dominates the inclination. Raising fixed extra light from zero to 30% improves the full forward score from 1419.5 to 1377.0, while the inclination shifts from 83.89 to the edge-on model boundary, the radius ratio from 0.79 to 0.95, and U only from 28.45 to 28.52 yr. Limb coefficients of 0 and 0.7 also favor the forward solution in the zero-dilution fit. Fractional extra light may absorb unmodeled shape or reduction effects. No light-source identification or precise inclination follows from this exercise.

For an external passband test, freeze the I-band geometry, set the V-band limb coefficient to 0.5 and recalibrate only the V brightness ratio and two zero points. The 95 V measurements give forward/reverse robust losses 107.2/123.3 with zero assumed extra light and 97.7/124.6 with 30%. The I and V exposures are distinct. Shared telescopes, sparse sampling, and the three V nuisance parameters limit this test.

## Falsifiable duration prediction

Predict the full width at half the lost flux of each eclipse on the dates used by TESS and 2027 July 1. The mean epoch for each TESS group is from the archived timing analysis, with no TESS brightness input in any geometric fit. For each assumed dilution fit in the range from zero to 0.3, vary the transferred limb coefficient from 0.2 to 0.6. The range refers to those choices and is not a statistical prediction interval.

| Epoch | Forward secondary minus primary FWHM (min) | Reverse (min) |
|---|---:|---:|
| 2018–2019 | −12.3 to −8.5 | +11.9 to +13.5 |
| 2020 | −9.6 to −6.2 | +9.4 to +10.6 |
| 2023 | −1.6 to +1.2 | +2.1 to +2.5 |
| 2025 | +4.1 to +7.1 | −4.1 to −3.5 |
| 2027 July 1 | +8.7 to +12.0 | −9.4 to −8.1 |

In 120 paired resamplings of OGLE observing nights with the geometry refitted (u=0.4), forward width differences for 2027 were +9.2 to +12.3 min and +8.1 to +11.1 min at zero and 30% assumed extra light, respectively (2.5th to 97.5th percentiles). Reverse differences were −10.9 to −7.0 and −11.6 to −5.3 min. The two senses predict closely similar center separations in 2027. The predicted signs of the width differences differ. High-cadence, resolved photometry covering both contacts, both eclipse profiles and an out-of-eclipse baseline over multiple orbits would test this prediction. Predicted physical flux widths and fitted Gaussian widths in a detrended TESS light curve have different observational responses. With independent-cycle resampling as a supplementary check, the 2025 TESS Gaussian width difference was +22 minutes with a 95% interval of −33 to +92 min. A block bootstrap that accounts for adjacent orbits is supplied separately in this investigation.

## Rate and out-of-eclipse shape

With U=28.45 years, the fitted advance is 0.04053 degrees per orbital cycle. The leading point-mass relativistic apsidal rate would require a total mass of approximately 746 solar masses to reproduce all of it. For illustrative total masses of 10 and 20 solar masses the fractions are 5.6% and 9.0%. These numerical examples are conditional estimates of the relativistic share, not measurements of the masses. rate.py also tabulates the effective classical internal-structure constant under additional illustrative assumptions of equal masses, aligned synchronous spins and classical tides plus rotation. None of those conditions is established by this photometry. [Giménez 1985](https://adsabs.harvard.edu/pdf/1985ApJ...297..405G) · [Claret et al., apsidal motion theory](https://doi.org/10.1051/0004-6361/202141484)

For a separate physical consistency check, mask every OGLE I observation at which the inferred disks overlap. Regress the remaining magnitudes against a separate baseline for each OGLE survey and the two sine/cosine pairs at the orbital and twice-orbital frequencies. Increasing the exclusion margin from one to 1.18 or 1.38 times the sum of the radii tests leakage from contacts. Six hundred within-block night bootstraps give uncertainty intervals for the second harmonic. The leading circular equilibrium-tide term for each star is proportional to the companion mass ratio, cubed fractional radius and squared sine of the inclination. Weight each term by the corresponding luminosity and the assumed constant dilution. The minimum over all positive mass ratios of A*q+B/q is 2*sqrt(A*B). Compare empirical and predicted semiamplitudes. For the zero-dilution geometry, 1,017 out-of-eclipse points give +0.45 millimag with a 95% interval of −2.82 to +3.68, compared with a minimum 5.84 millimag from the leading equilibrium formula at limb coefficient 0.2 and zero gravity darkening. At 30% assumed dilution the corresponding result is −0.50 (−3.75 to +2.84) and 5.08 millimag. A quasi-static eccentric template using the instantaneous star separation to the inverse third power gives +0.51 (−2.68 to +3.66) and −0.36 (−3.86 to +3.13) millimag. Additional independent third/fourth harmonics leave the observed cosine small. Complete fits that impose the lowest predicted coefficient, allowing the other stellar parameters to change, reduce the minimum expected signatures to 5.25 and 4.72 millimag at a cost of about 5.7 and 6.0 robust loss units against identical free-coefficient fits. Taking limb and gravity-darkening coefficients of 0.4/0.32 increases those costs to 13.5 and 14.4. The compensating parameter shifts leave the precession period around 28.4–28.6 yr. At these eclipse-derived radii, there is modest tension with the simplest version of the equilibrium model. Model incompleteness or nuisance orbital harmonics can influence that comparison. [Equilibrium ellipsoidal approximation and limitations](https://doi.org/10.1051/0004-6361/201015550)

## Replay

Python, NumPy and SciPy suffice. Run with OPENBLAS_NUM_THREADS=1. The initial geometry and the historical timings used for the retrospective test came from the earlier independent OGLE analysis. The shared compact TESS orbital sums and extraction receipt are under ../joint/inputs.

1. Run `python geometry.py`, then `python geometry.py --run cv`, `python geometry.py --run sensitivity` and `python geometry.py --run profiles`.
2. Run `python checks.py` for the dilution grid and V test and `python checks.py --bootstrap 120` for the paired night bootstrap.
3. Run `python predictions.py`, `python rate.py` `python ellipsoidal.py`, `python joint_tide.py` and `python width_boot.py`.
4. See ../widths/README.md for the TESS aperture and uncertainty checks. The original 400-replicate individual-cycle summary is retained in tess_widths.json.

All observational credits belong to the OGLE and TESS collaborations and STScI/MAST. Use results with the observing limits and assumptions described above. [MAST TESScut](https://mast.stsci.edu/tesscut/)
