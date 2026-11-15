# Source-level parallax check

Gaia DR3 returned all 1,367 membership entries for the eight clusters. Each row in `members.csv.gz` pairs a UCC membership probability with a Gaia source ID, the native Gaia fields and the Lindegren zero point in mas. The model depends on magnitude, colour estimator, astrometric solution type and ecliptic latitude. We subtract each valid source's model value from its published parallax. Both catalogue reductions and this calculation use Gaia observations. [Gaia DR3 archive](https://gea.esac.esa.int/archive/) · [Lindegren et al.](https://doi.org/10.1051/0004-6361/202039653) · [Correction software](https://pypi.org/project/gaiadr3-zeropoint/)

Astrometric centres use componentwise medians of stars with membership probability at least 0.5, parallax error below 0.18 mas, RUWE below 1.4 and a valid zero-point solution. For velocities, we use members with errors below 12 km/s whose velocities lie within the larger of 6 km/s and 4.5 median absolute deviations of this velocity sample's median. The baseline uses 12,000 resamplings of astrometric and velocity members, then adds an independent 0.010 mas normal offset per cluster per realization. This offset represents an assumed residual systematic. These straight-line calculations test the impact of star selection and distance calibration. The parent report compares trajectories integrated in a Galactic potential.

| Pair | Closest approach in Myr from present | Separation in pc, median [16%, 84%] | Fraction within summed Jacobi radii |
| --- | ---: | ---: | ---: |
| UPK 90 and UPK 93 | +3.37 | 7.80 [6.98, 8.58] | 1.000 |
| ASCC 99 and UPK 20 | −7.29 | 10.25 [9.45, 11.09] | 1.000 |
| ASCC 114 and UBC 168 | −2.86 | 10.38 [3.40, 21.92] | 0.850 |
| NGC 6633 and UPK 41 | −3.61 | 10.65 [9.84, 12.12] | 0.994 |

The fraction counts draws under these choices and treats present-day Jacobi radii as fixed. Raising the systematic floor to 0.020 mas lowers the ASCC 114 and UBC 168 fraction to 0.63. The CSV also includes a fixed −0.017 mas zero point on the same stars, a membership threshold of 0.8, and relaxed RUWE. The highest-ranked future encounter survives all four perturbations. Age relative to the inferred epoch remains important for UPK 41.

From this directory run `python sensitivity.py --output /tmp/flyby-sensitivity.csv` with NumPy and pandas. To refresh the Gaia columns and numerical zero points through ESA TAP, install `gaiadr3-zeropoint==0.1.0` and run `python fetch.py`. Source IDs must be read as strings or 64-bit integers. [Frozen membership basis](../inputs/ucc_members.csv)
