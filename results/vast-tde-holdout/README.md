# Sealed VAST holdout packet — review only

This extends the exact f323ab5 development packet; it does not alter it. Actual final sessions remain sealed until separate root review and authorization. `holdout_cli.py` checks the original manifest and code bytes, the exact development summary and all three development file hashes, plus a new reviewed code manifest and authorization record. The authorization must specify `phase: holdout`, `authorized: true`, the reviewed `git_head`, `manifest_sha256`, and `development_summary_sha256`. Only frozen control IDs and the two fixed sessions are eligible; failure produces an inconclusive record. Parquet predicates limit returned rows, not physical decoding of compressed column chunks.

| Source | Baseline mJy/beam | Variance | Interpretation |
|---|---:|---:|---|
| AT2020nov | 1.8474939960179828 | 0.010393621333138812 | prior radio source; association unresolved |
| AT2018hyz | 10.455068728099794 | 0.839456533117947 | known variable control; χ²/dof 385.63/8 |
| ASASSN-15oi | 4.884219194581238 | 0.19231210807883914 | known variable control; χ²/dof 264.93/8 |

For each future session, `z = (future_flux - baseline_mean) / sqrt(future_variance + baseline_variance)`; both z values must reach 3 for a conditional screen crossing. Future flux and variance are not known before access. The known controls fail the constant model strongly. Calibration, selection, blending and non-Gaussian tails remain uncalibrated, so neither result establishes physical significance.

## Development-only pseudo-tests

Before computation, 20 metadata control IDs per target were fixed by SHA-256 of their unpadded decimal IDs (`pseudotargets-frozen.json`, SHA-256 `aef29ad1b93ac5f960e61e46fb4ba67528c48d49e4bacfc48582641a72867cc9`). The last two development sessions were pseudo-holdouts; earlier development sessions trained the same algorithm, excluding each pseudo-target from controls. Independent +0, +0.5, +1 and +2 mJy/beam offsets were applied to each pseudo-target's two test sessions. Failures were retained.

The initial run returned 0/60 measured because its caller supplied controls the production phase guard correctly rejected; its output remains in the D: campaign. The caller was corrected to use only frozen controls, without changing selection, data or thresholds. The second run measured 1/20, 19/20 and 20/20 pseudo-targets for AT2020nov, AT2018hyz and ASASSN-15oi. Other cases lacked source-session coverage. Crossings at offsets 0, 0.5, 1 and 2 mJy/beam were 0/1, 0/1, 0/1, 0/1; 0/19, 2/19, 9/19, 17/19; and 0/20, 4/20, 16/20, 18/20. These share real backgrounds and controls; they are not an independent false-alarm calibration or an event-detection efficiency. AT2020nov's 1/20 coverage severely limits interpretation. `pseudo-results-v2.json` has SHA-256 `72f28886e7fd5907070f284aaa8a01997b3554fdd54bbb5e23615aef512e8874`; both runs are in `D:/Open-Space/campaigns/vast-tde-dev-f323ab5-21ccdd9-20261117`.

Run `python results/vast-tde-holdout/test_holdout_guard.py` for synthetic pre-read guards. From a clean committed checkout, `python results/vast-tde-holdout/holdout_cli.py --make-manifest D:/Open-Space/campaigns/vast-tde-holdout-REVIEW.json` creates an exclusive review manifest. Do not execute holdout processing without separate authorization.
