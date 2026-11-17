# Fourteen recent arcs that timed out in the blind fit

Thirteen of the fourteen objects have nominally bound inbound gravity solutions. The fourteenth is the known interstellar comet 3I/ATLAS, recovered at 57.9684 km/s. All fourteen independently fitted trajectories crossed 250 au with negative radial velocity. The weakest bound result is C/2024 R4 at +4.24 ± 3.97 × 10⁻⁶ au⁻¹. Its formal uncertainty admits a nearly parabolic original orbit.

## Completed

The input was the subset with a timed-out initial fit and frozen station/night positions in the [current MPC retry](../../batch/evidence/recent_stratified.json.gz). The independent fitter retained 2,468 of 2,491 sampled positions. The observation arcs span 335 to 4,906 days. The last column gives barycentric reciprocal semimajor axis at the actual inbound 250 au crossing with propagated formal one-sigma uncertainty.

| Object | SBDB SPK ID | Sampled | Retained | RMS per coordinate (arcsec) | 1/a ± σ (10⁻⁶ au⁻¹) |
|---|---:|---:|---:|---:|---:|
| C/2017 K2 (PANSTARRS) | 1003517 | 180 | 170 | 1.10 | 25.78 ± 0.68 |
| C/2021 G2 (ATLAS) | 1003820 | 180 | 178 | 0.55 | 29.05 ± 1.37 |
| C/2022 E2 (ATLAS) | 1003844 | 180 | 179 | 0.56 | 78.52 ± 1.15 |
| C/2022 N2 (PANSTARRS) | 1003861 | 180 | 176 | 0.99 | 8.01 ± 1.85 |
| C/2022 QE78 (ATLAS) | 1003889 | 180 | 178 | 0.45 | 23.54 ± 2.61 |
| C/2022 R6 (PANSTARRS) | 1003891 | 180 | 179 | 0.49 | 31.72 ± 2.61 |
| C/2023 H5 (Lemmon) | 1003935 | 180 | 180 | 0.57 | 22.86 ± 2.32 |
| C/2023 R1 (PANSTARRS) | 1003943 | 180 | 178 | 0.53 | 16.52 ± 1.88 |
| C/2024 G6 (ATLAS) | 1003998 | 180 | 179 | 0.61 | 46.97 ± 3.12 |
| C/2024 J3 (ATLAS) | 1004001 | 180 | 180 | 0.48 | 57.45 ± 2.69 |
| C/2024 R4 (PANSTARRS) | 1004028 | 180 | 180 | 0.47 | 4.24 ± 3.97 |
| 3I/ATLAS (C/2025 N1) | 1004083 | 180 | 180 | 0.69 | -3782807.12 ± 15.25 |
| C/2025 V2 (Rankin) | 1004102 | 180 | 180 | 0.59 | 49.35 ± 4.53 |
| C/2025 A1 (Lemmon) | 1004115 | 151 | 151 | 0.61 | 142.80 ± 3.06 |

All fourteen selected-observation hashes agree with the source selection. The SHA-256 of its compressed JSON file is `847fe716fee14893b6f8396cb5b590a27e9fc79fffea7115b7f9c07ede37c08e`. The [runner output](evidence/results.json) carries the individual hashes, force model, fitted states, covariances, residual statistics, signed radial velocities, and 50 and 250 au event times. Query receipts document the starting ephemerides, station list and planetary ephemeris. The computation is in [replay.py](replay.py). [Workflow run](https://github.com/dovvnloading/Open-Space/actions/runs/44034891617)

## Scope

These are gravity-only fits to the published station-night sample. The formal errors are conditional on its astrometry, weighting, and force law. The closest bound result has only 1.07σ clearance from zero. Nongravitational forces, activity-induced shifts and a full-archive fit can move small energies enough to change their sign. The predefined trigger for model and holdout follow-up was a non-control negative energy indicator `1/a < −2.5σ`. No non-control object met that trigger in this run. The remaining record in the source file had an acquisition failure and supplied no timed-out selection for this independent pass.
