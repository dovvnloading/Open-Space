# IMAP magnetic depletion screen

**No candidate passed the frozen five-minute screen.** During the first 12 UTC hours of 2026-07-31, 133 of 144 five-minute bins passed magnetic-data quality and cadence gates. Of 142 overlapping three-bin comparisons, 127 had three valid bins and 103 had both controls at least 3 nT. None then met the required 1 nT drop below the mean of the two adjacent controls. The overlapping triplets are correlated; these counts are not independent trials. The result is limited to this cadence, scale, product, threshold and interval, and does not show that shorter or smaller magnetic depletions were absent.

The predeclared question was whether the first qualifying localized magnetic depletion would coincide with increased proton scalar thermal pressure. A candidate would also have had to lie at least 20% below **each** control. No candidate was selected, so no proton prediction or test occurred. No SWAPI physical values or MAG values after 12:00 UTC were read. The SWAPI fits have some dependence on MAG direction; even a subsequent pressure result would not establish total pressure balance, a new source or a new physical mechanism. The older literature on solar-wind magnetic holes includes [Kasper et al.](https://doi.org/10.1029/2006JA012116).

| Gate | Count |
| --- | ---: |
| Five-minute bins valid | 133 / 144 |
| Three-bin comparisons attempted | 142 |
| All three bins valid | 127 |
| Both controls at least 3 nT | 103 |
| Mean-control drop at least 1 nT | 0 |

The exact input is public NASA SPDF [IMAP MAG L2 normal-rate GSE 2026-07-31 v002.0003](https://spdf.gsfc.nasa.gov/pub/data/imap/mag/l2/norm-gse/2026/imap_mag_l2_norm-gse_20260731_v002.0003.cdf), SHA-256 `f587e5ea2969418a54738580f8289ba2ef32a1d755b6af377c6ad3a8dc774727`. The first 86,400 records cover [00:00,12:00) UTC; 81,321 passed quality, interference and valid-value checks. A valid five-minute bin requires at least 540 retained records, edge coverage and no retained-sample gap above two seconds. The [development output](evidence/development.json) preserves every bin's counts and gates. The thresholds are engineering choices, not calibrated statistical significance; [IMAP MAG instrument documentation](https://doi.org/10.1007/s11214-026-01277-8) gives instrument characteristics.

The [frozen code and protocol](frozen/README.md) preserve exact bytes from commit `0695ac72304f0078306f52a45b67c346ed75aeec`, including the unexecuted, separately guarded SWAPI test. The [manifest](evidence/manifest.json), [authorization](evidence/root-authorization.json), [command](evidence/command.json), [execution receipt](evidence/execution.json), [result](evidence/development.json), and [input receipts](evidence/acquisition.json) preserve provenance. These local records are not external preregistration or cryptographic timestamping. An [independent implementation](evidence/independent.py) reproduced all bin and gate counts from the same MAG CDF; this is an independent calculation, not another observation.

## Offline replay

With Python 3.14 and NumPy, run `python results/imap-mag-swapi-screen/replay.py` to independently recount all bin audits and selection gates from the deposited [MAG projection](evidence/mag-projection.npz). The projection contains only first-12-hour times, field vectors and the combined quality/validity mask. It was made after the screen using the frozen reader by the [adapter](make_projection.py), so this recount does not independently validate native quality filtering. The [projection receipt](evidence/projection-receipt.json) records its SHA and the source hash. To regenerate from the verified external MAG CDF, install cdflib and run `python results/imap-mag-swapi-screen/make_projection.py MAG_CDF NEW_OUTPUT.npz`; the adapter verifies the source hash. Original synthetic CDF and model tests in `frozen/` are implementation checks, not measured false-alarm rates. The frozen Git and authorization guard refers to its original checkout and is not directly runnable from this publication copy; the standalone adapter and replay are explicitly post-screen.

Data are credited to NASA's [IMAP mission](https://imap.princeton.edu/) and the instrument teams; the original raw CDF files remain external to this repository.
