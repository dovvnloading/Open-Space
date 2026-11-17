# IXPE Mrk 501: a conditional repeatability screen

The later 2026 IXPE observation did not cross the predeclared threshold for a change in the selected broadband Stokes vector. The fixed-center comparison gives Δ(q,u) = (0.04436484, −0.00885483), χ² = 0.88363247 and conditional p = 0.64286776. A predeclared XY-centroid sensitivity gives Δ(q,u) = (0.06147492, −0.01675437), χ² = 1.78333708 and p = 0.40997113. Both exceed the predeclared p threshold 5.7×10⁻⁷. The p-values are conditional on the stated statistical and calibration assumptions; this result does not establish physical constancy, a population limit, or a new discovery.

## Question and observing units

Does the measured 2–8 keV source normalized Stokes vector (q,u) change between IXPE Mrk 501 baseline 05004899 and the separate later observation 05004901? The baseline is the aggregate of five segments 05004801–05 from 18 June through 12 August 2026; the later observation spans 28–30 August. NASA documents the [as-run dates and the meaning of suffix 99](https://ixpe.msfc.nasa.gov/for_scientists/asrun.html). The five baseline segments and their 99 product must not be counted twice. Good-time intervals were compared, and the test uses disjoint intervals of detector events. Different detector units are simultaneous measurements of each visit, not independent visits. The roughly 55-day baseline aggregates variable visits and can average away short changes.

We selected DU1 and DU3, calibrated level-2 event Q/U, event energies 2 ≤ E < 8 keV, and excluded the particle flag at STATUS array index 7, following the [HEASARC post-2025 example](https://heasarc.gsfc.nasa.gov/docs/ixpe/analysis/contributed/ixpe_protractor.html). We used a 60″ source aperture at ICRS 253.4675695168°, +39.7601691321° and a concentric 120–270″ annulus, applying the aggregate detector exposure-map ratio to local background. A five-iteration position-only centroid rule defines a sensitivity extraction; it uses the same events and is not an independent test. DU2 was a diagnostic and is excluded from the primary estimate; its calibration history is described in the [IXPE archive notes](https://heasarc.gsfc.nasa.gov/docs/ixpe/archive/).

The unweighted event estimator uses detector-specific, time-selected ixpeobssim effective areas and modulation factors at measured energy. Source and scaled background I/Q/U sums are combined before taking ratios; the covariance includes the event-sum I/Q/U covariance and ratio propagation with the scale fixed. The photon-weighted measured-energy estimator is not a response-folded true-energy estimate. The statistic is the quadratic form of the two-vector difference using the sum of the two statistical covariance matrices; under the conditional two-dimensional Gaussian approximation, p = exp(−χ²/2). No threshold was tuned after the later visit.

| Extraction | Baseline (q,u) | Later visit (q,u) | Conditional p |
| --- | --- | --- | --- |
| Fixed catalog position | (0.00377650, −0.04078458) | (0.04814133, −0.04963942) | 0.64286776 |
| XY-centroid sensitivity | (0.00125975, −0.03224243) | (0.06273467, −0.04899681) | 0.40997113 |

## Limits

The statistical covariance does not include a validated shared calibration systematic, uncertainty in background normalization, or adaptive-centroid selection. The baseline and later visit share the instrument, analysis, and calibration. The baseline’s source centroid was about 23″ from the radio catalog position in event sky coordinates; the cause is unassigned. The finite aperture, point-spread wings, detector and sky background gradients, exposure corrections, and energy redistribution can affect the comparison. A restricted response-matrix study and a separate shared-vector, power-law forward fit were development checks, not a validated universal bound or the primary estimator. The approximate Gaussian p-value is a screening statistic; neither result is a calibrated physical constraint on polarization constancy. Existing IXPE work established X-ray polarization of Mrk 501, e.g. [Liodakis et al. 2022](https://doi.org/10.1038/s41586-022-05338-0) and [2022–2023 follow-up](https://arxiv.org/abs/2407.11128).

## Evidence and reproduction

The [frozen pre-unblinding packet](frozen/conditional_packet.json), [baseline sufficient statistics](frozen/frozen_baseline.json), and matching frozen code preserve the bytes executed from the worktree at commit `0c4fa420585ae7a851c9a282e36795f10ccadadb`. Nine worktree files, including the hashed centroid rule, differed from that commit’s Git blobs only by line endings. Their [Git-blob versions](git-blob-differences/) are deposited separately; the frozen packet hashes identify the executed worktree bytes, so a fresh checkout of that commit alone does not reproduce every packet hash. [Protocol](frozen/PROTOCOL.md) contains historical development sections; its final conditional-packet section specifies the test. The packet itself retains its pre-unblinding “proposed” status for byte integrity. The [local authorization record](evidence/pre-unblind-authorization.json) was recorded before the local acquisition; it is not an externally registered commitment or cryptographic timestamp. The [baseline receipt](evidence/baseline-receipt.json), [later-visit receipt](evidence/receipt.json), [comparison](evidence/comparison.json), and [run log](evidence/run.log) preserve the observed input hashes and output. The recorded receipt SHA-256 proves consistency with acquired bytes, not an independent archive signature.

A standard-library check reconstructs both Stokes ratios and their covariances from the deposited event sufficient statistics, checks frozen file hashes, and recomputes χ² and p. It does not independently re-read raw events:

```sh
python results/ixpe-mrk501-repeat/check_offline.py
```

For an event-level replay, download the exact inputs from the URLs in both receipts and verify their SHA-256 values. Raw baseline support totals about 1.43 GB and is not deposited here; the later event, GTI and exposure-map files are also external. Install ixpeobssim 33.2.1 with NumPy, SciPy, Astropy, matplotlib, regions and skyfield; the packet records SHA-256 for each required time- and detector-specific response in ixpeobssim’s `caldb/ixpe/gpd/cpf` directory. Set `IXPEOBSSIM_DATA` and `IXPEOBSSIM_AUXFILES` to suitable writable campaign directories before executing. PyXspec was not used. The following **post-unblinding adapter** supplies paths to the frozen code and compares the full replay output to deposited evidence to 1e−9; choose a new output path:

```text
python results/ixpe-mrk501-repeat/replay_raw.py --baseline BASELINE_RAW --holdout HOLDOUT_RAW --caldb IXPEOBSSIM_CALDB_GPD_CPF --output NEW_OUTPUT.json
```

The exact native replay and standalone arithmetic check passed on Windows with Python 3.14 and ixpeobssim 33.2.1. The adapter was written after unblinding and did not determine the recorded result. Native files are public via [HEASARC IXPE](https://heasarc.gsfc.nasa.gov/docs/ixpe/archive/), and [ixpeobssim](https://ixpeobssim.readthedocs.io/en/latest/) is distributed by its authors; retain their licenses and the mission acknowledgement when reusing the data or software. No third-party raw event data are redistributed in this result folder.
