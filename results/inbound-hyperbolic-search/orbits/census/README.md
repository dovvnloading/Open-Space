# Gravity refits of the 578 disjoint long-arc catalogue entries

The independently fitted trajectories reached an inbound 250 au crossing for 378 of 578 assigned catalogue entries. Of the 378, 318 had positive reciprocal semimajor axis and 60 had negative values in the gravity-only model. Sixteen of the latter were more than three formal errors below zero. The [denser force-law follow-up](../robust/) tests the strongest signs for activity, data weighting and arc selection.

| Outcome | Objects |
|---|---:|
| Reached inbound 250 au with a fitted state and covariance | 378 |
| Fit produced an endpoint short of an inbound 250 au crossing | 3 |
| Fewer than ten usable ground station-nights | 179 |
| Too few retained positions, too short, or less than 70% of the sampled arc retained | 11 |
| Residual consistency test excluded too many measurements | 5 |
| Measured ground arc shorter than 30 days or fewer than ten positions | 2 |

## Completed

Every frozen SBDB catalogue row with apparent heliocentric excess speed below 1 km/s and reported arc of at least 30 days received a disposition, after excluding the 67 higher-speed rows and the separately selected recent and focused IDs. The primary-data guard required ten valid ground station-nights and an observed span of at least 30 days. Each station/night supplied its median-timestamp record. A deterministic time-ordered sample of up to 96 representatives preserved the endpoints. The 378 completed crossings used 24,712 retained positions in total.

The frozen SBDB heliocentric osculating conic supplied the initial state. The fitter adjusted all six phase-space coordinates against the measured sky positions and checked the fitted directions against an antipodal solution. It integrated backward from each fitted state using DE440 and the variational equations until the incoming 250 au event or the finite integration endpoint. Barycentric 1/a is positive for bound solutions. The quoted one-sigma errors are local linearized uncertainties conditional on the force model. The three endpoints short of an inbound 250 au event are C/1879 M1 (Swift), 205.8 au / C/1879 Q1 (Palisa), 242.8 au / C/1899 S1 (Giacobini), 249.7 au. [Independent fitter](../refit.py) · [Selection and replay](replay.py)

The [summary table](evidence/summary.csv) has one row per assigned SPK ID. Each [compressed shard](evidence/) contains the exact MPC80 sample, sample and raw input hashes, query receipts, fit statistics, state, covariance and event details, including failures. The frozen CmtObs snapshot came from workflow run 44029645452. These sixteen long-arc shards ran at code revision `6532c1ecf6abdec5467914582669fa9dc47c5479`. [Computation](https://github.com/dovvnloading/Open-Space/actions/runs/44035652465)

## Limits

This pass uses ground optical station/night samples from CmtObs, a local Gaussian covariance and a Newtonian gravity model. Active comets can depart significantly from this model. Selecting on a nominally positive heliocentric energy also does not imply a positive original energy. The underconstrained and missing-ground-data entries remain open. The separate [short-arc pass](../census_short/README.md) documents the assigned catalogue rows with a shorter or missing reported arc.
