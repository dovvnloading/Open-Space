# MMS spacecraft 1 planar-boundary screen

**No candidate passed the frozen screen.** In a 2026-09-30 MMS burst, 15,194 overlapping one-second windows had the required ten-second quality and cadence context; 11,125 also had both field endpoints at least 3 nT. None showed a rotation of 45° or more. These are overlapping windows, not independent trials. The largest eligible rotation was approximately 13.14°; that value is descriptive and did not set a new threshold. This is a bounded negative screen, not evidence that discontinuities or reconnection were absent.

The test asked whether the first qualifying sharp rotation on MMS1 could be modeled as a steady plane advected with the MMS1 ion flow, predicting timings and field changes on MMS2–4. No rotation qualified, so no plane prediction was formed and no held-out spacecraft field or flow values were read. An eventual plane match alone would not establish reconnection, electron diffusion, or a new physical mechanism. Selection is conditional on burst scheduling, instrument flags, thresholds and a single observation.

| Gate | Windows |
| --- | ---: |
| One-second windows tested | 16,351 |
| Ten-second quality and cadence context | 15,194 |
| Both endpoints at least 3 nT | 11,125 |
| Rotation at least 45° | 0 |

The input is MMS1 FGM burst L2 `mms1_fgm_brst_l2_20260930132803_v5.575.0.cdf` and FPI/DIS burst ion moments `mms1_fpi_brst_l2_dis-moms_20260930132803_v3.4.0.cdf`, from the [NASA SPDF MMS archive](https://spdf.gsfc.nasa.gov/pub/data/mms/mms1/). The frozen interval includes FGM native indices 33–16511, 2026-09-30 13:28:03.421179269–13:30:12.157399823 UTC. The native files contain 16,512 FGM and 866 DIS records; 16,479 and 858 fall in the frozen interval, with 16,473 FGM samples passing aligned quality and valid-range checks. FGM flags must be zero, as must DIS error flags, with nearest DIS time within 80 ms. The native CDFs remain external; [independent evidence](evidence/independent.json) records their full SHA-256 values.

The [frozen code and protocol](frozen/README.md) are exact bytes from commit `6ebb8aa4dbc40ea4afcfa9a7e9ac8618864824d5`, reviewed before spacecraft 1 access. The [manifest](evidence/manifest.json), [authorization](evidence/root-authorization.json), [result](evidence/development.json) and [execution receipt](evidence/execution.json) retain the local ordering and hashes. These local records are not an external preregistration or cryptographic timestamp. The [independent implementation](evidence/independent.py) independently reproduced the counts from the same native inputs, not independent observations.

The [acquisition receipts](evidence/source-receipts.json) retain the exact archive URLs, sizes, full-file hashes and range-query headers for all eight frozen MMS files. A range-prefix digest is distinct from a full-file digest. The other spacecraft supplied time, quality and position metadata; their field and flow values were not analyzed.

## Offline reproduction

With Python 3.14 and NumPy, run `python results/mms-planar-screen/replay.py`. It independently recounts the first four gates from the deposited compressed [MMS1 projection](evidence/mms1-projection.npz). This projection was generated after the screen by the [adapter](make_projection.py), which checks the original CDF hashes and uses the frozen SC1 input filtering; its [receipt](evidence/projection-receipt.json) hashes the projection and gives input lineage. The projection contains cropped times, magnetic vectors and a combined validity/quality flag; it does not independently establish the native filtering. To regenerate it, install `cdflib` and NumPy and run `python results/mms-planar-screen/make_projection.py PATH_TO_MMS1_FGM PATH_TO_MMS1_DIS NEW_OUTPUT.npz`; the adapter refuses inputs whose names or hashes differ. The NPZ hash may depend on serialization details, so the deposited projection is the fixed replay artifact. Original frozen tests in `frozen/` are synthetic software checks, not observations or calibrated detection-efficiency tests.

The field and plasma measurements are mission products from NASA's [MMS mission](https://lasp.colorado.edu/mms/sdc/public/); see the [FGM product guide](https://spdf.gsfc.nasa.gov/pub/data/mms/documents/FGM_DataProductsGuide2016Apr20.pdf) and [FPI dataset guide](https://lasp.colorado.edu/mms/sdc/public/datasets/fpi/). The original frozen code's local Git and campaign-path guards refer to its original checkout and cannot be run directly from this publication copy without reconstructing those inputs and authorization; the standalone replay above is explicitly post-screen.
