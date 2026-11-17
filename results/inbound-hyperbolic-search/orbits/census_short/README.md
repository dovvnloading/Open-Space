# Remaining short or missing published arcs

The frozen catalogue supplied 1,615 disjoint IDs in this stratum. Fourteen measured-ground trajectories reached the 250 au inbound boundary with a finite, positive formal energy uncertainty. Four further fits did not reach that boundary. Two of the fourteen show a nominally negative 1/a beyond three formal errors. C/1996 J1-A has an independent [pre-fragmentation parent fit](../fragment/README.md) consistent with zero energy. C/1886 J1 (Brooks) was already reported hyperbolic in Brady's 1965 integration. [Brady, Table II](https://adsabs.harvard.edu/pdf/1965AJ.....70..279B)

| Final disposition | IDs |
|---|---:|
| Actual inbound 250 au measurement | 14 |
| Fit without an inbound 250 au crossing | 4 |
| No supported ground optical station-night in matching MPC data | 1,455 |
| No unique matching packed-designation record in CmtObs | 116 |
| Too few distinct ground station-nights | 22 |
| Too few distinct ground positions | 4 |

## Completed

This stratum comprises frozen SBDB heliocentric e ≥ 1 entries with apparent excess speed under 1 km/s and reported arc under 30 days or missing. The [long-arc pass](../census), higher-speed set, latest MPC inputs and individually targeted fits have disjoint selection from this pass at the frozen-catalogue ID level. The source snapshot matches the one used for the [blind Find_Orb fits](../../batch).

The primary-data guard required six nonduplicate ground-field sky positions, three distinct station/UTC-night groups and a measured span of two days. The fit retained at least 70% of the sampled time span, with the same minimum position, group and span requirements. Up to eight positions per station-night preserved angular rate information, followed by deterministic time-order thinning at 96 positions if necessary. Six initial-state coordinates were varied independently of the sign of the energy. All successful crossings include the modeled time-of-crossing derivative in the formal covariance. [Fitter and selection](../census/replay.py) · [Stratum entry point](replay.py)

The [summary table](evidence/summary.csv) gives one disposition per SPK ID. Four [frozen shard archives](evidence) include the selected MPC80 positions, input hashes, source and station receipts, initial-state epoch, covariance and true boundary states. [Workflow run](https://github.com/dovvnloading/Open-Space/actions/runs/44037019982)

## Limits

The archive includes spacecraft-coded and alternative-format records outside this ground optical fitter. An insufficient-ground-data disposition leaves that orbit unresolved. Sparse arcs can admit distinct orbit branches and non-Gaussian energy posteriors. For example, C/2009 K3 has nominal original 1/a −0.453 ± 1.335 au⁻¹, and C/2006 M2 has −0.052 ± 0.177 au⁻¹. Their apparently large nominal excess speeds have little statistical support. No rate or sky-population upper bound follows from these catalogue dispositions.
