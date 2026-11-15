# Kepler-359’s apparent three-planet commensurability

Individual Kepler transits place the phase combination λb − 6λc + 5λd near 252° and show slow motion over the mission. Across 128 random, joint transit-timing posterior draws, 123 trajectories made multiple full turns or excursions through all angles. The other five encountered a dynamical stopping boundary. This measurement tests whether a four-year near-commensurability identifies a persistent three-body lock. The posterior integrations cannot exclude unsampled configurations or dynamics outside the published coplanar model.

## Measured transit phases

We fit a separate weighted linear ephemeris to the individual Holczer transit times of each planet within each window. Integer transit number supplies the measured orbital phase at transit. The combination of the three fitted phases estimates the transit-phase angle at the window midpoint. It approximates the mean-longitude combination for low eccentricity. The directly measured quantity is the slope of this phase combination. Each half-window slope is expressed per 1,470 days for comparison.

| Window, BJD − 2454900 | Accepted transits b/c/d | Midpoint phase | Drift per 1,470 days |
|---|---|---|---|
| 60–1550 | 40/20/12 | 252.13 ± 0.20° | −13.01 ± 0.93° |
| 60–805 | 17/9/4 | 256.37 ± 0.41° | −20.38 ± 3.16° |
| 805–1550 | 23/11/8 | 249.00 ± 0.18° | −9.47 ± 1.28° |

We include measurements with positive timing uncertainty, zero outlier flag and zero overlap flag. Each independent-planet fit uses its quoted timing uncertainties, scaling its parameter covariance upward when reduced χ² exceeds one. These formal errors describe the fitted linear trends within each window. Dynamical timing variations can produce correlated residuals. The sample was selected for slow motion, so these uncertainties do not estimate a survey-wide detection significance.

## Dynamical test

The integrations use joint osculating samples from Jontof-Hutter and colleagues. The Jontof-Hutter model already measures interactions in the system. The tested three-body angle involves all three planets and uses no longitude of pericenter. For the comparison system Kepler-60 we use the established primitive combination λb − 2λc + λd and the same sampling rules.

| Nominal 300-year experiment | Kepler-359 | Kepler-60 |
|---|---:|---:|
| Random posterior draws | 128 | 128 |
| Angle traversed ≥720° with ≥1 net turn | 118 | 0 |
| Full occupied range below 180° | 1 | 128 |
| Intermediate range | 4 | 0 |
| Encounter or escape stop | 5 | 0 |

All four intermediate Kepler-359 trajectories accumulated multiple net turns by 3,000 years. The one initially confined trajectory, archive row 3749, expanded from 152° at 300 years to 284° at 3,000 years. Over 30,000 years its unwrapped range reached 2,222°, with only −0.30 net turns. It visited all angular phases repeatedly while changing direction. This sensitivity to the integration interval is why net winding or a short interval alone cannot establish libration. Kepler-60's full angular span at 300 years ranged from 12.0° to 53.5°, with a median of 18.7° and a median arc midpoint of 44.9°.

We initialize Jacobi elements at BJD 2455680, assuming coplanar edge-on orbits viewed along the positive y-axis. The input columns give the planet/star mass ratio in Earth/Sun units, the osculating period, e sin ω, e cos ω and first transit time in BJD − 2454900. We convert the transit true anomaly π/2 − ω to mean anomaly. Comparing IAS15 trajectories with the 72 actual Kepler-359 midtimes gives median absolute normalized residuals between 0.66 and 0.96 across 16 draws. This checks the element conventions against the observations used by the source analysis.

WHFast uses 24 steps per inner orbit with angle samples every 45 days. The first 16 random draws per system also run at 48 steps and 11.25-day output intervals. All 32 paired classifications agree. Among the completed nominal draws the largest adjacent sampled angle changes are 40.8° for Kepler-359 and 3.3° for Kepler-60. We mark trajectories unresolved for an adjacent phase step above 90° or sampled relative energy drift above 0.0001. Encounter stops occur at the smaller of the initial adjacent mutual Hill radii or a barycentric distance over ten outer semimajor axes. We do not interpret those stopped trajectories. For the five completed Kepler-359 draws outside the nominal multiple-turn criterion, the extended integration uses 48 steps per inner orbit and samples every 19.5 days. The 30,000-year case uses 45-day outputs with a maximum adjacent step of 2.31° and maximum relative energy error of 4.2 × 10⁻⁸.

## Reproduction

The included input subsets allow offline recomputation. The reference run used NumPy 2.3.5 and REBOUND 5.2.2. Run the following commands in this folder.

```sh
python -m pip install numpy==2.3.5 rebound==5.2.2
python analyze.py --workers 4
python fit_check.py
```

`phase.csv` gives the fitted slopes and `orbits.csv` contains each random row's result and its numerical diagnostics. The two longer integrations are in `extended.csv` and `long.csv`. For circulating solutions the midpoint of an unwrapped range is arbitrary and should not be interpreted as a libration center.

To rebuild the archived inputs, run `python prepare.py`. The code verifies the original posterior archive SHA-256, retains each original zero-based row index, and selects 128 rows without replacement for each system with seed 20261115. The VizieR subset includes all returned rows for both systems. Holczer tabulates only planet c of Kepler-60, so the direct three-planet drift is computed only for Kepler-359. Input SHA-256 hashes are `56f6de3e4c1c3c682c6451e29801b3b48f7e5ee69ae939bca78fc907337b6311` for `posterior.csv` and `035a90468d3e89a8b36996d40ac634cdb9c12d587f0a64b81860f6509ae43a19` for `transits.csv`.

## Sources

- Holczer et al. (2016), individual Kepler transit times, Table 3. [VizieR catalog](https://cdsarc.cds.unistra.fr/viz-bin/cat/J/ApJS/225/9)
- Jontof-Hutter et al. (2021), posterior samples and model conventions. [Paper](https://arxiv.org/abs/2101.01202) and [archived posterior](https://doi.org/10.5281/zenodo.4422053)
- Goździewski et al. (2016), the three-body resonant configuration of Kepler-60. [Paper](https://doi.org/10.1093/mnrasl/slv156)
