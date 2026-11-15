# Two measured orbits give HD 73712 a hierarchy problem

In the Praesepe triple HD 73712 (CHR 130, HIP 42542), a 48.717-day spectroscopic pair orbits inside a resolved visual pair. A later speckle observation rules out the sky separation predicted by the preferred 2017 orbit. The two newer visual solutions reproduce that observation, but at their nominal eccentricities the outer star comes uncomfortably close to the inner pair. The newer published solution gives only **2.22–2.93 inner semimajor axes** at closest approach over the tested stellar mass ratios. The 2024 solution brings the outer star **inside one inner semimajor axis**. The measured quantity is the stars' separation on the sky. The close passage is inferred from the fitted orbits.

| Visual solution | Period (yr) | Eccentricity | Predicted separation at 2024.1538 | Periastron in au | Periastron / inner semimajor axis |
| --- | ---: | ---: | ---: | ---: | ---: |
| Tokovinin 2017 | 35.45 | 0.169 | 61.1 mas | 21.77 | 37.49–49.54 |
| Tokovinin 2024, Circular 213 | 20.96 | 0.980 ± 0.027 | 128.5 mas | 0.285 | 0.64–0.84 |
| Tokovinin et al. 2026, Table 4 | 20.83 | 0.930, tabulated error zero | 128.1 mas | 0.925 | 2.22–2.93 |

SOAR measured **128.8 mas** at position angle 160.0° in 2024. The 2017 orbit misses it by 67.7 mas. Its residual grew from −3.9 mas in 2018 to +14.5, +28.1 and +67.7 mas in 2021, 2022 and 2024. These four positions postdate the 2017 fit. Both later solutions used the 2024 visit, so its agreement with them is in-sample. Absolute separations above use Gaia DR3 parallax 5.3522 ± 0.1950 mas. The source has RUWE 4.995, and the angular semimajor axis of the 2024 solution has a 32% formal uncertainty. The last column uses the periods and mass ratios alone. [Earlier orbit and its published aliases](https://doi.org/10.3847/1538-3881/aa8459) · [2024 elements](https://www.astro.gsu.edu/wds/bsl/cir213.pdf) · [Latest SOAR measurements and solution](https://doi.org/10.3847/1538-3881/ae6cf5)

## What is measured

Abt & Willmarth measured two sets of moving sharp lines and a broad-lined component. In their 16 exposures with both sharp-line velocities the components are anticorrelated (Pearson r = −0.949), with separations of 40.4–60.5 km/s. Their published inner orbit reproduces the signed velocity differences with a median absolute residual of 5.46 km/s. On HJD 2448674.829 the published component labels give the opposite sign, a 90.4 km/s residual. The authors described the inner period as tentative. A simple circular sinusoid fitted to the paired differences also permits a 56-day alias. We preserve every row, including that exposure. The spectra establish changing line velocities. The assignment of the sharp pair to visual B follows the published hierarchical interpretation. [Original spectra measurements and orbit](https://doi.org/10.1086/307569) · [VizieR measurement table](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/ApJ/521/682)

For outer-to-inner mass ratio q = M₃/(M₁+M₂), Kepler's law gives the outer periastron in units of the inner semimajor axis as (1−e) × (P_outer/P_inner)^(2/3) × (1+q)^(1/3). We tested q from 0.3 to 2 and mutual inclinations from 0° to 180°. The 2026 periastron is 30–46% of the Mardling–Aarseth empirical stability threshold. The 2024 value is 7–10%. For these periods the threshold corresponds to eccentricities of roughly 0.82–0.87. That threshold is an empirical guide. The mass ratios and orientations span scenarios, without implying a probability distribution. The 2026 table quotes no fitted uncertainty around 0.930. [Stability criterion](https://doi.org/10.1046/j.1365-8711.2001.03974.x)

## Direct gravitational check

We initialized point-mass Newtonian triples at outer apastron with q = 0.7, inner mass ratio 0.86 and inner eccentricity 0.11. We used the period of each published visual orbit. A control with outer eccentricity 0.75 uses the 2026 period. Each group contains twelve configurations spanning three mutual inclinations (0°, 90°, 180°) and four inner true anomalies (0°, 90°, 180°, 270°). The lines of apsides and nodes were fixed by the coordinate construction in the script. We integrated to 0.72 outer periods, past the first nominal periastron, and computed the instantaneous two-body period of the original inner pair from its relative position and velocity.

| Tested outer eccentricity | Original pair unbound or its period changed by >10% at both tolerances |
| ---: | ---: |
| 0.980 (2024) | 12 / 12 |
| 0.930 (2026) | 2 / 12 |
| 0.750 (control) | 0 / 12 |

The 2026 outcomes span 36.8–52.6 days at the tighter tolerance. The ten cases below the 10% threshold and the limited time span leave long-lived configurations open. Every control ended within 0.002 days of its starting value. DOP853 ran with relative tolerances 10⁻¹⁰ and 10⁻¹². Maximum absolute fractional total-energy error across all 72 runs was 2.5×10⁻⁷. Some individual high-eccentricity trajectories diverge between tolerances after close scattering, so the reported count requires the threshold to be crossed at both. These controlled examples describe the response to the assumed geometry. They do not constitute a fitted three-body trajectory or an escape-probability estimate.

## Interpretation and test

Tokovinin had identified the triple, the nearly edge-on geometry and ambiguity among possible outer periods by 2017. Five years later Morales et al. reviewed the blended spectra in the context of Praesepe stars. The addition here is a quantitative physical test of the *subsequent* 2024 and 2026 visual elements against the spectroscopic orbit, with later relative positions testing the old solution. A search of SIMBAD references, the orbit circulars and indexed primary papers using “HD 73712”, “CHR 130”, “HIP 42542”, “48.717”, “stability”, “periastron”, “orbit” and “eccentricity” on 2026-11-14 (US Eastern) found the earlier discussion of orbit ambiguity. Priority for the present comparison remains unestablished. [Tokovinin 2017, target discussion](https://www.ctio.noirlab.edu/~atokovin/papers/orbits2017.pdf) · [Morales et al. 2022, Appendix A](https://doi.org/10.3847/1538-3881/ac7329)

The nearly edge-on orbit leaves a nearly straight projected track. Sparse close-passage coverage, a 180° ambiguity in older speckle position angles, and the fact that the instrument measures light from an unresolved subsystem can all affect eccentricity inference. The available fits do not measure the full orbital covariance. The simplest interpretation is that the current elements inadequately describe the stable outer motion. An improved inner ephemeris and velocity measurements of the *broad-lined star relative to the sharp-lined pair's centre of mass*, coupled with speckle measurements that identify which quadrant contains the fainter source through the predicted 2032 passage, could distinguish eccentricity, orientation and component assignment. Future spectroscopy can directly test whether the 48.7-day inner clock changed across a passage.

## Replay

Run from the repository root.

```sh
python results/hd-73712-orbits/verify.py
python -m pip install numpy scipy
python results/hd-73712-orbits/three_body.py
```

The first command uses the standard library. The second program prints one JSON record per simulation and tolerance, comparable to `three_body_results.jsonl`. `radial_velocities.csv` contains only HD 73712, with the authors' broad-line, sharp-primary and sharp-secondary labels b, p and s. Dates are published heliocentric Julian dates and velocities are km/s. `later_astrometry.csv` contains the four speckle measurements in the comparison, with published epoch, separation in mas and position angle in degrees. Dates in the visual-orbit calculation use the decimal years of the source tables. `sources.json` lists observation tables and SHA-256 digests where available. A VizieR server-generated timestamp changes its response digest on a fresh request. The derived velocity CSV is frozen with this research commit.

The target emerged from an exploratory USNO ORB6 conjunction screen, outside the repository's prior campaigns. Of 4,051 orbital solutions, 3,339 had grades 1–5, periods in days or years of at least 0.2 yr, a numeric primary magnitude, angular semimajor axis in arcsec or mas, inclination, eccentricity and argument of periastron. HD 73712 ranked seventh by approximate conjunction separation divided by a magnitude-based angular-size proxy. The proxy was used only to choose targets. The triple consistency test was selected after reading the source's spectroscopy. No survey-wide frequency or significance follows from this selection. To replay the screen, supply the frozen ORB6 SQL-format file from `sources.json` as an optional positional argument to `verify.py`; the verifier checks the full-file SHA-256 before ranking.
