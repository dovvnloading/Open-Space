# Independent referee audit — GWTC-2.1/3/4.0 strong-lensing pair census

**Verdict: scientifically acceptable as a bounded, calibrated *summary pre-screen* only.** No validated lensed pair or full-posterior catalogue-wide exclusion. No discovery-level label or remarkable-findings log entry. This audit deliberately distinguishes the published LVK null conclusion from what this pipeline itself has demonstrated.

## Independent numerical replay

The standalone `referee_replay.py` does **not** import the analysis script. It independently reconstructs each frozen event's mass–spin summary score and lens-time-delay bin, checks all unique pairs against the compressed deposited table, and reruns 1,500 catalogue-stratified shuffles with the frozen seed. It verifies the indexed Zenodo file MD5 and length against the per-download receipt, recomputes the SHA256 of each extracted NPZ, and checks the original HDF5 SHA256 was recorded. The runner's original HDF5 SHA256 is **not** independently verified by a second download. The per-run receipt is `referee_checks.json`.

- 141 unique eligible events: 34 GWTC-2.1, 23 GWTC-3 and 84 GWTC-4.0. This is a **selected** subset, not every published GW detection.
- 9,870 distinct unordered pairs (141×140/2). All pair delays and scores exactly replay from frozen inputs (numeric tolerance 1e-10).
- Strongest *summary proxy*: GW190720_000836 / GW190728_064510, separation 8.27539 days and score 4.49145 (a log-score/density, **not** log Bayes factor).
- Empirical max-statistic background: 956/1,500 shuffled trials at least as high; `(956+1)/(1500+1)=0.6375749500333111`. This is a family-wise *summary proxy* false-alarm fraction, **not** the false-alarm probability of a full-posterior overlap BF.
- The frozen posterior sample files were individually verified against Zenodo checksums before extraction. See `outputs/posterior_receipts.json` for original HDF5 SHA256 and extracted NPZ SHA256. **Sampling selection is not an independent catalogue population sensitivity study.**

## Published controls: external conclusion is not a reproduction

LVK's O4a lensing study reports no credible lensed pair; its leading posterior-overlap pair was GW231113_200417 / GW231223_202619, whose joint Hanabi log10 Bayes factor was **−3.3** against lensing. Our simpler time-delay/mass/chi-eff proxy ranked it **second** among the 3,486 O4a pairs, **not first**. That is qualitative ranking agreement only, not reproduction of LVK's overlap Bayes factors or joint-analysis odds. One earlier O3 control (GW190421_213856/GW190910_112807) has source PE and can be checked; two published O3 control pairs require GW190424_180648, which is outside the eligible frozen records here. Mark those unavailable, never silently pass them. The official published Bayes factors are **external** evidence, not numbers we have recomputed.

## Sky, frames and neglected evidence

The independently binned equal-area sky samples of the top **proxy** pair show little common high-probability sky region on a ~5-degree grid, even though detector-frame component masses have overlapping credible intervals. The diagnostic changes with pixelization and smoothing: it is **not** a stable sky BF, an image-parity test, or a formal posterior-odds conclusion. Use the pipeline's declared HEALPix estimator for its own limited posterior comparison, and do not overinterpret its approximate prior.

Mass consistency under strong lensing must use **detector-frame** component masses; source-frame masses depend on inferred distance and redshift and can differ between magnified images. Spin vector components, inclination and sky coordinates require a consistent common PE prior and proper waveform/frame conventions. Approximate KDE priors, smoothing, nonuniform selection and target pre-selection are unresolved.

## Required before any discovery-level or full-catalogue null

1. Exact PE prior reweighting and full-posterior overlap BFs for a frozen *complete* eligible catalogue, including quantitative LVK control reproduction with waveform/systematic robustness.
2. A lens-population and detection-efficiency calculation, including joint selection and robust full-posterior trials/global FAP, not transfer of the summary score's background FAP.
3. Direct Morse phase/image type compatibility and a reduced coherent strain likelihood from GWOSC for survivors; detector-frame mass and spin frame-consistency confirmation.
4. Credible foreground galaxy / EM lens candidates in overlapping sky region, delay geometry, alternate lens models, and independent external validation.

**Disposition:** No claimed lensed discovery. Publication-ready only as a transparently incomplete methods/negative-control screen. No revision of `docs/REMARKABLE_FINDINGS.md` is justified.

References: [LVK O4a lensing study](https://arxiv.org/abs/2512.16347) · [Zenodo released O4a lensing controls](https://zenodo.org/records/18163632) · [Haris et al. posterior-overlap method](https://arxiv.org/abs/1807.07062).
