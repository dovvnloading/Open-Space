# Measurements, selection and replay

The observing units are a catalogue row for screening, a Gaia DR3 source for follow-up, and a dated exposure for RV validation. An exposure-specific pipeline comparison can expose a bad fit. RV agreement over multiple nights tests a different failure mode. The DESI, SDSS and literature directories contain their respective queries, source identifications, controls and receipts. The top-level `SHA256SUMS` covers the frozen contribution.

## Selection

1. We evaluated all 4,140,838 DESI DR1 SpecDis v2.1 rows. Eligible rows have a positive Gaia DR3 ID, a five- or six-parameter solution, 0 < RUWE < 1.4, finite position, proper motion and RV, a positive spectrophotometric distance with 0 < fractional error < 25%, and the survey's single-star flag. This last flag is a catalogue classifier. An undetected binary can pass it. The precise-parallax branch also requires `(parallax − zero_point) / parallax_error ≥ 5`. Nominal speeds use full positions and velocities.
2. The first stage retained rows with nominal speeds ≥550 km/s on both the precise-parallax and SpecDis distances, or ≥900 km/s on SpecDis alone. We followed the remaining 550–900 km/s SpecDis tail separately. Repeated source IDs can occur. The combined screen covers 2,209 eligible rows at or above 550 km/s. `desi_screen.csv` contains 1,020 first-pass rows. `moderate/desi_screen.csv` contains 1,187 additional rows for 1,185 additional sources after excluding first-pass source IDs. Intrinsic, survey-wide completeness is unknown. A star with a poor distance, a rejected fit, or a nominal speed below the threshold can be unbound.
3. The RV validation matches the delivered catalogue RV to a DESI coadd by target ID, Gaia DR3 ID and velocity, requires successful stellar classification, zero RV warnings, 0 < RV error < 20 km/s, and agreement with `c × Redrock z` within 60 km/s. We fetched SPDist by both source and target ID, applying its good-input and parameter-consistency flags. The SPDist uncertainty includes the recommended 0.167 mag intrinsic term. We advanced the nearby survivors and stars with Bailer–Jones photogeometric median speed ≥550 km/s, together with quality-eligible candidates in the separate line-of-sight search. In the extended tail, 69 stars received posterior calculations. The primary radial lead's 25.9% fractional SpecDis error put it outside the initial distance cut. The independent radial-velocity selection recovered it.
4. The independent DESI selection covers 4,640,216 MWS coadds. Its 580 km/s cut uses a projected Galactic-rest-frame velocity, independent of distance. The screening frame and final frame differ slightly. The +659 km/s raw screening value for the leading source becomes +657.0 km/s after the DESI backup correction and the adopted frame transformation. We retained lower-speed repeated observations and inspected bright-galaxy associations as checks. The DR19 search separately selected SDSS-V spectra, joined Gaia IDs, and screened both radial and transverse velocity. It includes a no-parallax radial branch and a deterministic all-speed control. Cut definitions, excluded populations and denominators are in the survey reports. [DESI](desi/README.md) · [SDSS](sdss/README.md)
5. The 21-source literature check compared source IDs and sky positions for 2,708 distinct targets. The catalogue index includes historical fast stars, several Gaia-release searches and recent DESI work. Ten have exact or possible prior associations. Published physical-parameter catalogues can also report a high velocity or positive energy without entering an HVS list, as the DR19 MINESweeper example demonstrates. Completeness outside the listed tables, aliases and positional windows is unmeasured. [Sources, coverage and comparison rules](literature/README.md)

## Phase-space posteriors

We used Astropy's right-handed Galactocentric frame with Sun–centre distance 8.122 kpc, solar height 20.8 pc and solar velocity (+12.9, +245.6, +7.78) km/s. X increases from the Sun towards the Galactic centre, Y follows disk rotation locally, and Z points north. Solar parameters are fixed. The central 68% intervals in the tables are marginal quantiles of the Cartesian samples. The median vector need not have the median speed.

For each star and distance treatment, `infer.py` draws 4,096 samples with a seed derived from the integer Gaia ID and model name. It constructs the five-parameter Gaia covariance, adds floors of 0.02 mas in parallax and 0.02 mas/yr in each proper motion, and treats the RA error as an angular displacement in RA × cos(Dec). The parallax-based models condition the four other astrometric parameters on the draw's true parallax using the Gaussian conditional covariance. The parallax-free models use the marginal four-parameter covariance. A Gaussian RV likelihood includes the reported uncertainty and a 1 km/s floor, or a 2 km/s floor for corrected backup velocities. These are analyst-chosen sensitivity assumptions. Formal probabilities condition on the epoch RV representing the systemic velocity. All source IDs are read as strings or 64-bit integers.

DESI parallax zero points come from SpecDis. The DR19 branch evaluates the Lindegren five-/six-parameter correction using magnitude, colour and ecliptic latitude. Two contaminated six-parameter comparison stars fall outside the calibrated pseudocolour range and are explicitly flagged. For distances, we used SpecDis alone, SpecDis with a Gaussian corrected-parallax likelihood, a parallax-only likelihood with an exponentially decreasing space-density prior of scale 1.35 kpc, and published Bailer–Jones geometric and photogeometric posteriors. SpecDis and valid SPDist measurements enter as positive lognormal approximations. Their stellar training data include Gaia, so correlated modelling errors remain possible. Bailer–Jones, StarHorse, GSP-Phot and MINESweeper posteriors enter individually through log-space approximations to the reported quantiles. We do not multiply those posteriors by the parallax likelihood a second time. The GSP-Phot and MINESweeper formal uncertainties can be much smaller than disagreements between atmosphere models. [SpecDis](https://data.desi.lbl.gov/doc/releases/dr1/vac/mws-specdis/) · [SPDist](https://data.desi.lbl.gov/doc/releases/dr1/vac/mws-spdist/) · [Bailer–Jones](https://doi.org/10.3847/1538-3881/abd806)

We added a no-parallax single-dwarf colour–magnitude check for six DESI sources, with 20% fractional distance uncertainty, increased to 25% for dereddened g−i < 0.4. The flow through `candidate_models.csv` shows how assigning a dwarf, a subgiant, or a brighter unresolved system can change the escape classification. These rows quantify conditional hypotheses about luminosity. The RRc control uses its published systemic RV and class-specific absolute-magnitude scale. The DR19 giant has a separate gravity-based parallax-free estimate. [Photometric input measurements](literature/followup/dwarf_distances.csv) · [DR19 distance diagnostics](sdss/README.md)

## Potentials and origins

For each draw we evaluate escape to infinity from the draw's actual 3D position, using `sqrt(-2 Phi)` in the isolated McMillan 2017 potential (galpy) and Gala's `MilkyWayPotential(version='v2')` (2022 model). The two models have different halo and disk parameterizations. For the two radial prospects we also evaluate the published heavy-halo Irrgang Model III. Its fit required that another fast halo star be bound. No probability weights over potential families or binary states have been assumed. [McMillan](https://arxiv.org/abs/1608.00971) · [Gala model](https://gala.adrian.pw/en/stable/supporting/define-milky-way-model.html) · [Irrgang et al.](https://doi.org/10.1051/0004-6361/201220540)

For candidate origins, a kick–drift–kick integration follows 256 posterior states backward for 300 Myr in Gala v2. A moving Hernquist LMC of mass 1.5 × 10¹¹ solar masses and scale length 17.14 kpc accelerates the stars and the MW origin. Each draw varies the LMC centre, distance, proper motions and RV. We test minimum 3D separations of 1 kpc from the Galactic centre and 1 and 5 kpc from the LMC centre. We interpolate the first plane crossing and subtract the local circular velocity to estimate the disk-relative crossing velocity. The model omits dynamical friction, a Galactic bar, time-dependent MW mass and structural parameter uncertainty. These approximations limit identification of a physical launch site. Zero events in 256 draws corresponds to a one-sided 95% binomial upper limit of about 1.16% under the adopted model. Halving the 0.25 Myr step reproduced median closest approaches and crossing times for the leading source to the reported precision. [LMC potential scale](https://academic.oup.com/mnras/article/507/4/4997/6346557) · [S5-HVS1 measurements for the positive control](https://academic.oup.com/mnras/article/491/2/2465/5612212)

## Replay

All commands in this contribution, including survey screening and archive replay, require Python 3.12. From this results directory, a fresh Python 3.12 environment can replay the posteriors and orbit tables using the frozen input rows. The scoped `.gitattributes` preserves the exact deposited file bytes in fresh checkouts, including Windows checkouts with automatic line-ending conversion enabled. `verify.py` checks the hashes, sample integrity, positive control and headline counts. The survey-specific reports provide upstream query replay. Install Poppler's `pdftotext` to rebuild the literature index from its source papers.

```sh
python -m pip install -r requirements.txt
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python verify.py
python infer.py
(cd moderate && python propagate.py)
python orbits.py
python candidate_tables.py
python heavy_halo.py
python systemic_rv.py
python convergence.py
python positive/check_s5.py
python photometry.py
python plot_rv.py
```

The `*_16`, `*_50`, `*_84` columns give marginal posterior percentiles. `p_mc17` and `p_gala22` are fractions of samples above escape speed. In the orbit tables, plane-crossing quantiles use only draws with a crossing. `p_plane_R_lt25` uses all draws as its denominator. `inference.csv` includes audit subjects and two moderate-screen priorities, so concatenate the two inference files only after deduplicating `(source_id, model)` pairs.

To rebuild the archive-level DESI selection, obtain the three large FITS files from `inputs.json`, verify their SHA-256 hashes, and set `SPECDIS`, `SPDIST`, and `BACKUP_RV` to the **absolute paths**. The Gaia TAP, NOIRLab and survey endpoints are public. Exact bytes and SHA-256 receipts accompany the frozen queries.

```sh
python screen_desi.py "$SPECDIS"
python desi/mws_fetch.py desi_screen.csv desi/queries
python fetch_gaia.py desi_screen.csv gaia
python fetch_bj.py
python check_distances.py "$SPDIST"
python make_queue.py "$SPECDIS" "$BACKUP_RV"
python make_sdss_queue.py
(cd moderate && python prepare.py "$SPECDIS" && python ../desi/mws_fetch.py desi_screen.csv desi/queries && python ../fetch_bj.py && python ../check_distances.py "$SPDIST" && python review.py "$BACKUP_RV" && python ../fetch_gaia.py queue.csv gaia --key source_id)
```
