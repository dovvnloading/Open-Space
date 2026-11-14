# Spectroscopy and the outer orbit of TIC 96100077

Gaia recorded large radial-velocity changes in this source. Its combined orbit predicts that a spectroscopist following the same absorption lines would measure an outer-orbit semiamplitude of **26.15 km/s**. The available spectra do not identify the moving star or link it to the 7.93-hour eclipses. The catalogue coordinates are 264.526816°, −49.349215° and the Gaia DR3 identifier is 5946882743944130304.

## What Gaia measured

Gaia DR3 used all 16 of its accepted velocities over 618.14 days and 404 of 439 along-scan astrometric measurements. The source-level robust velocity range is 49.20 km/s. Seven RVS transits required window deblending. The median velocity is 9.59 ± 4.01 km/s. The joint eccentric-orbit fit gives a systemic velocity of 13.22 ± 0.89 km/s. Both quantities use Gaia observations. The 2.77 ± 1.92 km/s in DR2 is a median over five visits from an overlapping mission interval. [Gaia source data model](https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_main_source_catalogue/ssec_dm_gaia_source.html)

The 272.720 ± 0.762 day solution has bits 6, 13 and 56 set. Bit 6 records RV availability. Bit 13 means the RV-only search failed its significance threshold. Bit 56 means a spectroscopic trend or stochastic solution entered the combination with astrometry. The combined fit has eccentricity 0.5741 ± 0.0278 and goodness-of-fit statistic 4.32. The tabulated significance 40.73 pertains to the astrometric semimajor axis. A period determination from an independent velocity series would test the combined solution. [Gaia NSS fields and flags](https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_non--single_stars_tables/ssec_dm_nss_two_body_orbit.html) · [Gosset et al.](https://doi.org/10.1051/0004-6361/202450600)

## An internal orientation check and a velocity prediction

Gaia's combined fit shares period, eccentricity and periastron epoch, while allowing separate astrometric and spectroscopic arguments of periastron. The astrometric angle has a 180° ambiguity. The angular difference after allowing that ambiguity is **+0.81°**, with a formal 16th–84th percentile interval of **−1.83° to +3.91°**. This agreement is a consistency check within one fit. It supports an orbit common to the photocentre and RV signal, while leaving the association with the eclipsing pair open. [Gaia astrometric processing, sect. 5 and appendices](https://doi.org/10.1051/0004-6361/202243969) · [Gosset et al.](https://doi.org/10.1051/0004-6361/202450600)

| Quantity | At catalogue elements | Formal 16th–84th percentiles |
|---|---:|---:|
| Outer inclination | 93.55° | 92.42–94.57° |
| Photocentre semimajor axis | 0.4245 AU | 0.4128–0.4410 AU |
| RV tracer semimajor axis | 0.5378 AU | 0.4921–0.5874 AU |
| RV tracer semiamplitude | 26.15 km/s | 23.62–28.98 km/s |
| Velocity at minimum over a complete orbit | 1.72 km/s | 0.79–2.50 km/s |
| Velocity at maximum over a complete orbit | 54.02 km/s | 48.64–59.82 km/s |

The calculation uses the parallax from the orbital solution and the complete 15-parameter correlation matrix in the documented column order. The tabulated efficiency is 0.0. An archive query finds this exact value for 23,296 of 33,467 AstroSpectroSB1 solutions. The reconstructed correlation matrix for this source is positive definite, with smallest eigenvalue 0.01564. We draw 150,000 local Gaussian perturbations of the reported elements using that matrix. The intervals propagate the local fit errors under the published orbit. They omit alternative periods, component misidentification and systematic errors. [Gaia NSS model](https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_non--single_stars_tables/ssec_dm_nss_two_body_orbit.html) · [Query receipts](archive-queries.json)

## Can the eclipsing pair supply the missing light?

Assume Gaia follows a single star with the eclipsing pair on the other side of the outer barycentre. Also assume the negative timing coefficient measures only light travel time. Then its absolute value s is both the pair-to-single barycentric axis ratio and the single-to-pair mass ratio. The photocentre position is the light-weighted mean of those two positions. Calling the pair's effective fraction of the outer system's G light β gives |1 − (1 + s)β| = 0.7893 at the nominal elements. The absolute value leaves two possible fractions. [Astrometric and spectroscopic binary geometry](https://doi.org/10.1051/0004-6361/202347111)

| Assumed absolute clock scale | Pair fraction, photocentre on RV tracer's side | Pair fraction, opposite side | Pair's outer RV semiamplitude, opposite in phase to the RV tracer |
|---|---:|---:|---:|
| 0.800 | 0.117 (0.078–0.149) | 0.994 (0.963–1.033) | 20.9 km/s |
| 0.965 | 0.107 (0.071–0.136) | 0.911 (0.882–0.947) | 25.2 km/s |

The clock scales are fixed examples from existing TESS reductions. Their uncertainties and their dependence on the Gaia timing template are outside these intervals. Fractions above unity are physically excluded. The interval reports the unconstrained transformation to show where the second branch reaches that boundary. On the high-light branch the pair supplies nearly all G light, raising the question of which components dominate the absorption lines in the RVS band. The associated illustrative dynamical masses use the same two-body assumptions. [Machine-readable calculation](orbit-audit.json)

The two Gaia model minima lose 10.77% and 10.10% of the maximum G-band flux. Both include ellipsoidal variation. Removing each Gaussian eclipse from the model at the corresponding minimum restores 6.91% and 6.25% of maximum flux, respectively. Those two distinct losses sum to 13.16%, an illustrative light requirement if they occult opposite stars of unchanged brightness. The source's fitted eclipse model has reduced χ² = 17.58, and the conversion of Gaussian components into physical occultations is approximate. Comparing the two full depths directly to the astrometric fraction would also mix a maximum-light normalization with the mean of a variable photocentre. [Gaia eclipsing-binary model](https://doi.org/10.1051/0004-6361/202245330)

For a constant tertiary, sampling all inner phases uniformly gives ⟨Fmax/F⟩ = 1.0420 for the Gaia model. Under this approximation the pair's fraction at maximum light becomes 1 − (1 − β)/1.0420. The small-light branches above then give 15.27% (11.50–18.29%) and 14.32% (10.87–17.09%) at maximum light. Both can accommodate the nominal Gaussian losses, albeit with deep intrinsic eclipses. Gaia's actual epoch-dependent astrometric weights and ingress shapes could change that comparison. A light-curve fit that measures the tertiary contribution and spectra that separate the line systems could distinguish the branches.

Gaia's `binary_masses` row assumes an isochrone/luminosity primary mass and the two-body NSS model. It gives primary 1.407 (1.333–1.472) and secondary 1.255 (1.086–1.429) solar masses. The secondary-to-primary G flux ratio is 0.110 (0.059–0.163), or 0.099 (0.055–0.140) of combined G light. The assigned components and flux fraction share the Gaia fit and stellar assumptions. Using these figures as independent masses of three stars would double-count information. [Gaia binary mass fields](https://gaia.aip.de/metadata/gaiadr3/binary_masses/) · [Gaia stellar multiplicity](https://doi.org/10.1051/0004-6361/202243782)

## Archive coverage and prior identification

| Archive or survey | Target check on 2026-11-14 |
|---|---|
| Gaia DR3 and Focused Product Release | Both epoch-RV tables return zero rows for the source ID. The DR3 source also has `has_rvs=false` and `has_epoch_rv=false`. |
| GALAH DR4 and RAVE DR6 | No epoch record by Gaia identifier in either survey. Positional queries also found no matching spectra. |
| APOGEE-2 DR17 and SDSS-V | No DR17 allStar entry within 15 arcsec and no returned target in the 15-arcsec DR19/DR20 Valis queries. |
| ESO | No public reduced-spectrum SSA match around the target. Raw science and reduced-product centre searches within ±0.008° in RA and ±0.005° in Dec returned zero records. The Gaia–ESO DR5.1 stellar RV catalogue also returned zero entries in that box. Different pointing centres or unindexed fibre targets can escape these checks. |
| MAST and Las Cumbres Observatory | The MAST spectrum search returned zero entries within 15 arcsec. LCO's public point-coverage spectrum query returned an empty results page and an estimated total of one, so it supplies no usable spectrum. |
| Survey of Surveys I | A 0.005-arcsec match in the radial-velocity catalogue lists a single contributing survey, Gaia. Its RV includes no independent ground spectrum. |

Gaia publishes a mean low-resolution BP/RP spectrum with 343 samples from 336 to 1020 nm. Gaia MSC's two-template fit gives 5967 K (5646–6154 K) and 5982 K (5509–6913 K). It assumes an unresolved coeval binary and restricts the brighter/fainter BP/RP flux ratio to below five. The broad, overlapping ranges do not count stars. The resolving power of BP/RP, roughly 20–100, also precludes using the mean spectrum to resolve 26 km/s orbital motion. [Sampled source spectrum](data/xp-mean.csv) · [Gaia MSC assumptions](https://gea.esac.esa.int/archive/documentation/GDR3/Data_analysis/chap_cu8par/sec_cu8par_apsis/ssec_cu8par_apsis_msc.html) · [Gaia XP spectroscopy](https://doi.org/10.1051/0004-6361/202556362)

ESA schedules Gaia DR4 for **2 December 2026**. Its 66-month baseline extends beyond DR3, and the planned epoch-velocity release includes processed sources down to G_RVS = 14. This source has G_RVS = 10.20 in DR3. A release-day position/identifier crossmatch can establish source-specific availability. Visits after the DR3 observation interval would test the velocity prediction with observations outside its original fit. [ESA release schedule](https://www.esa.int/articles/news-media/esa-s-highlights-for-2026) · [Gaia epoch-data announcement](https://www.cosmos.esa.int/web/gaia/iow_20260703)

The original 2006 ASAS eclipsing-binary catalogue lists `173806-4920.9` with period 0.330413 d, V = 11.01, amplitude 0.13 mag and semidetached classification ESD. SIMBAD currently lists TYC 8351-1727-1 as a spectroscopic binary. The Villanova TESS EB lookup returns no entry for this TIC. These are bounded checks of prior cataloguing. The testable claim concerns the physical association of the known eclipse clock and the known Gaia outer orbit. [ASAS original catalogue](https://www.astrouw.edu.pl/asas/data/Vars.E) · [TESS EB target lookup](https://tessebs.villanova.edu/96100077) · [Archive responses](archive-queries.json)

Run `python orbit_audit.py` with Python 3 and NumPy to reproduce [orbit-audit.json](orbit-audit.json). The [archive receipts](archive-queries.json) include the query URLs, returned rows or response details, and SHA-256 digests for each archive check.
