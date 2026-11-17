# Euclid Q1 dark strong lenses

**No validated dark deflector emerged. One published grade B candidate has an unusually high, redshift-dependent mass-to-light ratio and detectable central light.** Its identification as a strong lens and the source and deflector redshifts remain unconfirmed. The result adds a mass-to-light measurement and an independent search of the Q1 pixels.

## The high-ratio system

**102158888_2659753895653871052** is at RA 265.975390°, Dec +65.387105° (ICRS, catalogue centre). Its published model has effective Einstein radius 0.936″. After subtracting the lensed-source model, the 0.3″-radius core is detected at 49.0σ in VIS and 22.2σ in H. The H-band light inside the Einstein aperture is 21.91 AB. Assigning all observed H flux within that aperture to the foreground lowers the inferred ratio by 11%. The independently released MGE and Sérsic total H magnitudes agree to 0.013 mag.

| Foreground redshift model | Source redshift | Projected H-band M/L within the Einstein radius, median (16th–84th percentile) |
|---|---:|---:|
| Euclid merged PHZ, median 0.17 | 2, assumed | 105.2 (89.0–123.1) |
| Legacy Surveys DR9, median 0.32 | 2, assumed | 54.4 (35.1–76.3) |

Units are solar mass per solar luminosity in the correspondingly blueshifted H passband. This convention compares the galaxy to the Sun through the same rest-frame filter. The two rows are separate conditional analyses. For the Legacy redshift model, all 25 redshift-nearest comparison lenses have lower ratios at the same assumed source redshift. Their median is 2.64 and their maximum is 18.32. Redshift, source identification, foreground extinction, PSF and lens-model uncertainties limit the interpretation. The system was selected from the tail of this sample. The comparison is descriptive.

![VIS and H images of the grade B lead, before and after foreground-light subtraction](imaging/102158888_2659753895653871052.png)

Only one distinct bright residual knot is available for a colour test. The mass and light centroids differ by 0.214″. The Q1 spectroscopy pipeline's nominal solution for the blended object has low probability and fails the published Hα quality cuts. The 3.5″ DESI DR1 cone contains no spectrum. A dust attenuation of about 3.29 mag in observed H would reconcile the Legacy-based ratio with the comparison median at the stipulated redshifts; the photometry does not measure that attenuation. The anomaly warrants follow-up spectroscopy of both the crescent and the central galaxy. [Posterior samples summarized by band and redshift](inference/conditional_ML_posterior.csv) · [Matched comparison population](inference/redshift_matched_peers.csv) · [Spectroscopy and mass-to-light inference](inference/README.md)

## What the Q1 data constrain

| Search | Observed result |
|---|---|
| Published SLDE Q1 | The main sample has 497 grade A/B candidates. Three additional grade A/B serendipitous systems bring the in-footprint total to 500. All 312 systems with a usable positive fitted Einstein radius have source-subtracted VIS core S/N ≥19.74. These are detections of central flux given the published decomposition. |
| Remaining 188 A/B systems | Every parent received a search for a bright ring with a faint centre using release cutouts or native IRSA fallback imaging. One configuration passed both VIS filters. Its trial centre is detected at 6.4σ in H and 5.4σ in J. Failed or missing catalogue mass models cannot give measured Einstein masses. |
| Secondary catalogues | The comparison includes 367 AgileLens entries, 36 AstroVink positions, and the 72-object bright/low-redshift extension in the spectroscopic and geometry audit. Catalogue entries overlap. The faint-H AgileLens grade B ring has a fitted radius of 0.566″. Its fitted centre is displaced about 2.2″ from the catalogue position and contains 18.3σ VIS and 23.4σ H light. |
| Independent sky search | Thirty-two preselected native VIS strips across the three deep fields yield 0.08337 deg² of eligible search centres. Of 1,258 first-stage morphology triggers, six passed the radial-profile test. Each received four-band imaging, local photometry and common-source tests. No configuration met the discovery standard. |

The published parent selection required a bright extended VIS detection. A dark deflector with bright arcs can enter through those arcs, but many dark deflectors would be absent. The sky search did not require any catalogue detection. Its fixed ring scales, coverage masks and preference for azimuthally extended emission restrict sensitivity. The 63.1 deg² Q1 footprint supplies the survey context; the smaller measured area is the denominator for this independent search.

A separate grade A object at (265.530607°, +65.463319°) illustrates how redshift can manufacture an extreme ratio. The Euclid merged photo-z is 0.02 and the Legacy photo-z is 0.721. At their nominal medians the source-redshift-independent H mass-to-light floor changes from about 131 to 1.92. Conversely, the grade A luminous control at (270.085761°, +65.267680°) has a four-image position fit with effective radius 1.339″ versus the published 1.356″. Its 16 image-band measurements all exceed 3σ and are consistent with a common four-band spectrum in the stated aperture test. [Independent fit and cutout](catalog/README.md)

## Measurements and reproduction

- [Catalogue reconciliation, photometry and independent four-image fit](catalog/README.md)
- [Calibrated image measurements, deflector subtraction and mass-centre registration](imaging/README.md)
- [Mass, light, redshift priors, spectroscopy and conditional posteriors](inference/README.md)
- [Independent sky search, all shortlisted configurations, limits and injections](arc_search/README.md)

This work has made use of the Quick Release (Q1) data from the Euclid mission of the European Space Agency (ESA). The [Q1 release](https://doi.org/10.57780/esa-2853f3b), [Strong Lensing Discovery Engine models and images](https://doi.org/10.5281/zenodo.15025832), and [IRSA MER mosaics](https://doi.org/10.26131/IRSA601) are the primary inputs. The result tables and scripts identify exact products and hashes.
