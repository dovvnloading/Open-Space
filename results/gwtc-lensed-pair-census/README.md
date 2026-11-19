# GWTC O1–O4a strong-lensing pair census (new dataset)

**Claim (pre-registration; no candidate yet):** This investigation tests whether any pair of independently catalogued gravitational-wave events can be defended as two gravitationally lensed images of one merger. A large posterior-overlap score is **not** a discovery. The intended null is no pair surviving an empirical false-alarm calibration and the published LVK lensing controls.

**Status:** Dataset manifest, code, frozen inputs and independent review are **pending**. Do not cite this pre-registration as a search result. New analysis only; no reuse of `results/gwtc-successive-mergers`.

## Frozen plan

- Primary releases: GWTC-4.0 O4a corrected PE v3 (Zenodo [17602505](https://zenodo.org/records/17602505)), GWTC-3 v3 ([22685054](https://zenodo.org/records/22685054)), and GWTC-2.1 v2 ([6513631](https://zenodo.org/records/6513631)). Event names, file URLs, Zenodo MD5 and locally verified SHA256 will be deposited with the results. Prefer detector-frame masses; lens magnification corrupts luminosity distance and inferred source-frame masses.
- Screen all frozen catalogue pairs with released samples where resource constraints permit, using `B_overlap = ∫p_1(θ)p_2(θ)/π(θ)dθ` (Haris et al. 2018), with sky and intrinsic factors; use a galaxy-lens delay distribution and independently shuffled real pairs or simulated unlensed pairs to calibrate the global false-alarm probability. Record explicit coverage and failure gates, never silently substitute summary statistics for full posteriors.
- Verify O4a published LVK controls from [GWTC-4 lensing release](https://zenodo.org/records/18163632) and [paper](https://arxiv.org/abs/2512.16347), including its top O4a pair GW231113_200417 / GW231223_202619: the published joint-analysis log10 Bayes factor is −3.3, **not lensed**. Prior LVK conclusions are null.
- For survivors check phase/Morse parity, time-delay and mass-frame consistency, strain (only if acquired from GWOSC), and credible lens galaxies in the sky overlap. A posterior-only screen cannot claim Morse-phase compatibility.
- Data and computation: one GitHub Actions job triggered by changes to source files on `research/gwtc-lensed-pair-census`. Generated `outputs/` are excluded from its push filter to prevent loops; runner uses network, pins hashes and commits bounded receipts. Failures and unfinished controls remain explicitly pending.

## Outputs and checks

| Item | Status |
| --- | --- |
| Frozen release and event manifest / verified hashes | Pending runner |
| Posterior overlap / delay population / empirical false alarms | Pending runner |
| Published LVK controls | Pending runner |
| Detector/source frame, phase and strain, EM lens checks | Pending; no claim without evidence |
| Independent referee review and reproducibility | Pending |

The labeling policy in [`docs/REMARKABLE_FINDINGS.md`](../../docs/REMARKABLE_FINDINGS.md) requires independent validation of a discovery-level result. No `[REMARKABLE]` claim or log entry is authorized by a screen or a null.
