# HARPS search for planets and brown dwarfs

**No previously unreported planet or brown dwarf qualified as a candidate.** The blind search fitted 2,192 peaks across 1,775 eligible HARPS-RVBank target series. It recovered 42.8% of injected circular 30-day orbits at 10 m/s. One low-amplitude signal recurred in two reductions of HD 17970, but the data could not distinguish a separate orbit from a nearly coincident harmonic. The survey places no general limit on the number of companions still present.

## What was searched

| Velocity series | Input spectra | Catalog target series | Series searched | Spectra on searched series | Keplerian fits | Injection trials |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| HARPS-RVBank, SERVAL | 252,615 | 5,239 | 1,775 | 200,827 | 2,192 | 1,068,672 |
| April 2026 CDS correction, retrievable subset | 50,538 | 983 | 375 | 41,856 | 458 | 227,616 |
| ESO Phase 3, DRS | 289,843 | 6,488 | 1,870 | 231,748 | 2,283 | 1,125,600 |

Each eligible target series had at least 12 observing nights spanning 30 days. The search tested periods from 1.5 days to that series’ baseline, fitted the strongest peak on every eligible series and could fit two further components. The three searches share HARPS exposures. Counts in each row stand on their own, and different archive names can refer to the same star. [RVBank target ledger](revised/stars.csv) · [ESO search and limitations](eso/README.md) · [RV methods](rv-method.md)

The author’s Ver.02b release spans BJD 2452936.5403–2459577.8410. Quality cuts kept 236,751 spectra, excluded 15,818 flagged by SERVAL and excluded 46 with invalid RVs or uncertainties. The cadence cuts excluded 3,453 target series for too few nights and 11 for too short a baseline. [Exposure ledger](revised/epoch-ledger.csv.gz) · [Input inventory](revised/ingest.json)

The April 2026 CDS correction was available for 50,538 consecutive Table 4 rows. All matched a unique author epoch. Of 50,537 finite RV pairs, 49.9% agreed within the declared precision. The comparison failed the velocity-scale agreement test, and the full-catalog search retained the author velocities. The other 202,077 author epochs lack a verified corrected counterpart in this run. We therefore ran a separate search using the corrected velocities directly. It retained 47,273 spectra and supplied usable S measurements on 35,900 exposures. The retrievable prefix follows catalog record order, which limits comparisons with the full sample. [Corrected-data search](corrected/README.md) · [Original CDS page receipts](corrected/inputs.csv) · [Period comparison by release](corrected/reduction-comparison.csv) · [CDS correction audit](revised/correction-audit.json) · [Per-field discrepancies](revised/correction-changes.csv.gz) · [CDS release history](https://cdsarc.cds.unistra.fr/viz-bin/ReadMe/J/A%2BA/683/A125?format=html)

## Companion checks

The final S-index review examined all 1,775 eligible target series and matched 167,480 public DACE exposures. It tested S-index for 2,035 of 2,192 fitted peaks across 1,639 target series. Of all peaks, 112 passed the internal RV checks, 1,949 were rejected and 131 remained inconclusive. The all-target S checks reproduced the same 112 passes found in the initial activity review. Among those passes the external screen found 85 prior-period matches, rejected 16 and left 11 for individual adjudication. Each of the 112 has a written disposition with sources. [Final S-index checks](followup/all-activity.csv) · [Coverage and query outcomes](followup/all-activity-summary.json) · [Source-by-source coverage](revised/source-coverage.csv) · [Review of every internal pass](survivor-review.csv)

Of the 20 internal passes in the corrected subset, 15 matched prior periods, two were rejected by the external screen and three required closer inspection. One of the three was GJ 479 at 11.2969 days. It lies inside the 11.279–11.302-day interval of an earlier proposed planet on a star with a measured rotation period near 24 days. The ESO fallback yielded zero internal passes and lacks the SERVAL activity channels. [Review of all 20 corrected-subset passes](corrected-survivor-review.csv) · [Corrected subset’s unresolved fits](corrected/inconclusive-periods.csv) · [Tuomi et al., Appendix A.30](https://arxiv.org/abs/1906.04644)

For each fit the ledger gives the period, semiamplitude, eccentricity, conditional false-alarm probability (FAP), available minimum mass, and activity, rotation, sampling and multiplicity checks. The analysis fits independent zero points across the 2015 fibre change and tests amplitude consistency between eras and between chronological halves. It tests FWHM, BIS, CRX, dLW, H-alpha and S-index, plus R′HK when available. Missing coverage stays visible. [RV fits and activity statistics](revised/signals.csv) · [Initial coverage by indicator](revised/activity-coverage.csv) · [Catalog and mass checks](revised/vetting.csv) · [Individual literature adjudications](review-notes.md)

Of 440 target series with a residual RV trend of at least three formal errors, 352 matched Hipparcos–Gaia accelerations and 168 had a four-component proper-motion anomaly with p < 0.001. No orbit met the criteria for reporting a true companion mass. The corrected subset yielded 89 trends, 78 astrometric matches and 35 anomalies at the same threshold, with zero adopted dynamical masses. A slope and an anomaly alone leave the orbit and companion mass degenerate. [Acceleration ledger](revised/accelerations.csv) · [Joint RV–astrometry method](vet-method.md)

The strongest later-date match was a known stellar binary. Gaia DR3 independently reports an SB1 period of 14.77334 days for BD−02 2880 and a secondary-mass lower limit of 0.214 solar masses. [Published Gaia orbits](followup/gaia-nss_two_body_orbit.csv) · [Gaia mass bounds](followup/gaia-binary_masses.csv)

The test on later ESO dates froze 251 RV waveforms before checking data from BJD 2459578 onward. Of those, 176 had no later exposures and 59 had fewer than 12 later nights. Only BD−02 2880 reached the signed-amplitude comparison. Its 14 nights occupied nine independent sign blocks, giving a best possible family-adjusted p of 0.4902. No waveform achieved the 0.05 family threshold. [All later-date outcomes](revised/holdout.csv) · [Holdout method](holdout-method.md)

### Two unresolved signals

**HD 17970 (GJ 3187).** The blind Keplerian fit gives a period of 72.4799 days, K = 1.783 m/s, eccentricity 0.295, conditional minimum mass 0.02896 Jupiter masses and FAP 1/256. It passed the six activity screens, including S-index on 59 nights, and the fibre-era and sampling checks. Published chromospheric activity gives an indirect rotation estimate of 49.7 ± 4.5 days. We found no independently measured rotation period or later ESO nights for this star. [Rotation and mass sources](vet-method.md) · [Later-date check](revised/holdout.csv)

On 68 nights with paired spectra, the circular 72.48-day amplitudes were 1.72 ± 0.21 m/s with SERVAL and 1.42 ± 0.26 m/s with the ESO DRS. A second fitted period of 362.132 days has a fifth harmonic at 72.426 days. Once that period and its free fifth harmonic were included, adding a distinct 72.48-day term worsened BIC by 7.3 in each reduction. It also worsened predictions for seasons left out of the fits. The nearly coincident frequencies and shared exposures leave the physical origin open. [Harmonic tests](revised/alias-review-summary.csv) · [Paired nights](followup/deep-paired-nights.csv) · [Fits and season comparisons](followup/deep-models.csv)

**HD 126525 (GJ 550).** The blind fit found 1,012.43 days at K = 5.86 m/s near the published 960.41-day orbit. A two-period fit improved in-sample BIC by 23.7 with broad period bounds. One eccentricity reached the fitted ceiling of 0.80, and the companion period moved from 946 to 816 days when the bounds changed. Eight nights qualified for the chronological prediction test and no independent later ESO nights were available. The observations leave a second orbit unresolved. [Published-orbit comparisons for four hosts](followup/prior-orbit-comparison.csv) · [Barbato et al., Table 15](https://boris-portal.unibe.ch/server/api/core/bitstreams/009da830-c2c1-4fea-9cb4-243d39837200/content)

## Sensitivity and limits

At 30 days and a semiamplitude of 10 m/s, the blind search recovered 6,073 of 14,200 injected circular orbits (42.8%). The corrected subset recovered 1,436 of 3,000 injections in the same cell (47.9%). Eight phases per target sampled the effect of cadence. The search recovered six of seven prespecified published orbital periods above the FAP and model-improvement thresholds. It recovered all three controls at HD 40307 and all three at HD 69830. It missed HD 4308 b. Some recovered controls also triggered activity flags, illustrating the cost of the final cuts.

The RV FAP uses the maximum over each frequency scan in 255 sign-randomized realizations of 14-day residual blocks. Its smallest reportable value is 1/256. It is conditional on the star, earlier fitted components and block-noise model. A separate Gaussian grid test includes a survey-wide trials penalty. Long-lived activity can violate the noise assumptions. Injection recovery uses the observed cadence and residuals and measures the initial periodogram stage. It cannot establish the fraction of all stars with planets. [Completeness grid](revised/completeness-grid.csv) · [Per-target injection cells](revised/injections.csv) · [Published controls](revised/published-controls.csv) · [Control fits by archive name](revised/published-control-series.csv)

## Input checksums and reproduction

The author’s RV archive has SHA-256 `9218ebd833f8971dcf304c7a6d9835de1c988dc5faae131f3eb939c7b9682586`. The ESO input used for the separate DRS survey has SHA-256 `19b2df8dd8816cca5facffd5614fd5ea6571f633310e8e7e6441eb39f2f7c138`. Each downloaded CDS page, successful external query and archived output has its own receipt. [RVBank and activity inputs](revised/inputs.csv) · [Catalog inputs](revised/catalog-inputs.csv) · [ESO holdout inputs](revised/holdout-inputs.csv) · [Output SHA-256 manifest](revised/SHA256SUMS) · [Installed versions](revised/packages.txt)

The full author-release search (run 44985728833) used source commit `fe930a83b8a88bfbcf9db310a33e858e9ced30ad`. The workflow records the exact commands for ingest, the 12 search shards, catalog checks, holdout and publication. Its push trigger accepts an edit to the workflow file on `survey/harps-rvbank-companions`. GitHub also accepts manual dispatch once the workflow exists on the default branch. The follow-up and corrected-subset search used commit `fc5a5d90b1e34ff3190fd2396acaf61a8ef45dca` in run 44988368934. The normalized RVBank and ESO inputs are archived with that run’s results, alongside their checksums. [Follow-up workflow run](https://github.com/dovvnloading/Open-Space/actions/runs/44988368934) · [Follow-up checksums](followup/SHA256SUMS) · [Corrected-output checksums](corrected/SHA256SUMS). The first author-release search and ESO fallback remain separately reproducible. [Survey workflow](../../.github/workflows/harps-rvbank-survey.yml) · [Full RVBank run](https://github.com/dovvnloading/Open-Space/actions/runs/44985728833) · [First author-release search](evidence/README.md) · [Source definitions](sources.md)

To check the deposited tables from the repository root, run

```sh
for stage in revised followup corrected; do
  (cd "results/harps-rvbank-survey/$stage" && sha256sum --check --quiet SHA256SUMS) || exit 1
done
python -m open_space.verify
python -m unittest discover -s tests
```

For the scientific replay, use the same workflow on the research branch. The registered manual invocation is

```sh
gh workflow run harps-rvbank-survey.yml --ref survey/harps-rvbank-companions -f scope=full
gh workflow run harps-rvbank-survey.yml --ref survey/harps-rvbank-companions -f scope=followup
```
