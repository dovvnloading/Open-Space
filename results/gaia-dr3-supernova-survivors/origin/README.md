# Galactic speeds, disk crossings and a remnant alignment

Seven published objects test the 3D machinery. We use the spectroscopic velocities and uncertainties adopted by El-Badry et al. (2026), full Gaia covariance matrices and the parent repository's phase-space and orbit routines. This keeps the frame and escape calculations consistent with the earlier DESI search. The measured Mg shift for J1109+0001 and the Shen fits are separate checks on the spectroscopic inputs. [Literature RVs](../control/published_kinematics.csv) · [Independent fits](../spectra/README.md) · [Shen spectra](../shen/README.md)

| Published object | Galactocentric speed, median [16th, 84th], km/s | Unbound fraction, McMillan 2017 | First disk crossing, median radius and lookback |
|---|---:|---:|---|
| D6-1 | 1,991 [1,819, 2,164] | 1.000 | 6.51 kpc, 0.62 Myr |
| D6-2 | 1,126 [1,082, 1,192] | 1.000 | No crossing in 256 draws over 50 Myr |
| D6-3 | 2,303 [1,849, 2,874] | 1.000 | 3.69 kpc, 2.22 Myr |
| LP 40−365 | 834 [829, 839] | 1.000 | 6.45 kpc, 4.42 Myr |
| J1109+0001 | 974 [722, 1,449] | 0.950 | 7.72 kpc, 3.60 Myr |
| J1240+6710 | 247 [228, 266] | 0.000 | 4.38 kpc, 25.57 Myr |
| US 708 | 921 [908, 933] | 1.000 | 5.85 kpc, 12.86 Myr |

The reference rows use Bailer-Jones geometric distances apart from US 708. For US 708 we use the 8.5 ± 1.0 kpc spectroscopic distance summarized in the 2026 paper and the revised +835 ± 10 km/s velocity. We approximate that distance likelihood by a lognormal and combine it with the corrected parallax. The alternate distance-model runs illustrate the role of priors, especially for J1109+0001. Its median speed varies from 974 to 1,636 km/s across the three models. The Bailer-Jones photogeometric prior represents ordinary field stars. [El-Badry et al.](https://arxiv.org/abs/2606.11293) · [Velocity and escape draws](states.csv)

The Sun is at Galactocentric radius 8.122 kpc and height 20.8 pc, with Cartesian velocity (12.9, 245.6, 7.78) km/s in Astropy's frame. The 4,096-draw speed summaries condition the proper motions on distance, account for astrometric correlations and include an RV floor of 1 km/s. The Bailer-Jones distance quantiles are approximated by a split lognormal. The third distance prior has an exponentially decreasing space density with scale length 1.35 kpc. We evaluate escape against the McMillan 2017 and Gala v2 Milky Way potentials. These fractions are conditional on the measured velocities representing systemic motion and the adopted masses and distances. For each distance model we integrate 256 backward trajectories for 50 Myr in the Gala potential including a moving LMC, at 0.25 and 0.125 Myr steps. The table uses the finer step. A plane crossing admits a disk origin. The progenitor, its motion and the ejection mechanism require independent evidence. [Orbit ensembles](crossings.csv) · [Existing frame and orbit methods](../../dr1-dr19-unbound-stars/METHODS.md)

### G70.0−21.5

Of 7 × 310 star–remnant comparisons over 0–200,000 yr, the only projected intersection with a catalogue enclosing circle is the previously published D6-2 association. Closest angular approach to the catalogue centre is **0.905° at 93,750 yr**. The tabulated remnant measures 330 × 240 arcmin, so this is inside even its 2° minor semiaxis. This short-baseline check projects the heliocentric rectilinear trajectory with perspective from the measured RV, holding the present shell centre fixed. It uses a 250 yr grid and central astrometric values. The quoted time is a grid minimum, not an age with 250 yr precision. Shell morphology, age, distance and unknown bulk motion affect the interpretation. No chance-overlap significance or new birthplace is inferred. [Green catalogue](https://www.mrao.cam.ac.uk/surveys/snrs/snrs.data.html) · [Shen et al.](https://arxiv.org/abs/1804.11163) · [Pair results](remnant_checks.csv)

Run `python origin/run.py` and `python origin/remnants.py` from the study directory with Astropy, Gala, galpy, NumPy and SciPy installed. The imports resolve against the sibling repository result, and the Gaia and Green receipts identify the input bytes.
