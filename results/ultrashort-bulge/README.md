# Ultrashort bulge microlensing

**No new Earth-mass lens or angular Einstein radius has been established.** Four of 2,429 OGLE bulge alerts have catalog timescales below half a day. Independent fits leave all four with a longer solution, a competitive transient or binary lens, or insufficient coverage to measure the finite-source radius. A separate test on the published terrestrial-mass candidate OGLE-2016-BLG-1928 recovers a 39.1-minute Einstein time and a normalized radius of 3.21 under a uniform, unblended-source model. The discovery authors measured 41.5 minutes and 3.39. This event is an external validation control.

This is a research checkpoint from 17 November 2026. [OGLE fits and measured coverage](ogle/README.md) · [KMT photometry, contaminants and cross-matches](kmt/README.md) · [Control inputs and fit](controls/earth-recovery.json)

## Measured cases

| Alert | Result |
|---|---|
| OGLE-2025-BLG-0867 | The binary-lens fit improves chi-square by about 285 relative to the finite-source single lens in the OGLE peak window. |
| OGLE-2025-BLG-1151 | The last eight observations sample the rise, followed by a 0.683-day gap. A timescale of at least half a day costs less than one unit of chi-square. |
| OGLE-2025-BLG-1201 | A finite-source fit gives 0.075 d and rho 1.32. Imposing rho ≤ 1 costs only 1.9 in chi-square. The binary and flare fits both improve on the finite-source single lens. |
| OGLE-2026-BLG-0843 | The continuous OGLE peak-window fit favors 3.80 days. The constrained short fit costs 180 in chi-square. |
| KMT-2024-BLG-0456 | Two observatories confirm the peak. The unconstrained finite-source model uses tens of times the reported baseline flux. With a source-flux ceiling of 1.15 times the photometric baseline, the best tested binary has chi-square 267 versus 449 for the finite-source single lens in the same window. |
| KMT-2025-BLG-1507 | Under the same ceiling, the binary gives chi-square 178 versus 338 for the best tested finite-source single lens. |
| MOA-2025-BLG-0198 | The alert team lists a cataclysmic variable. The public difference photometry records repeated brightenings in 2008, 2011, 2016 and 2025. |

## Frozen sources and ongoing checks

The 2024 and 2025 KMT catalogs list 3,441 and 3,348 public event IDs. The 2026 index lists 2,248 positions and preliminary fits. KMTNet marks the 2026 season proprietary. The 2025 and 2026 OGLE catalogs list 1,500 and 929 bulge alerts. MOA and PRIME contribute 1,669 public bulge alert records across 2024–26. Gaia's all-sky archive returned 89 alerts inside the stated bulge rectangle in 2024 and zero new alerts there in 2025–26. Different surveys can list the same event.

The PSPL, FSPL and restricted binary template fits, ingest receipts and model-selection tests are being assembled. Candidate-specific checks remain incomplete. An event alert list has no monitored-star denominator. A survey occurrence rate or an upper bound on rogue planets per star would require blind-search efficiency and monitored star-time, in addition to cadence.

## Model implementation

The shared fitter imports magnitude cleaning, flux conversion and point-lens magnification from the existing [OGLE work](../ogle-dark-lenses/lens.py). The grid profiles independent baselines and positive source fluxes per stream at the native epochs. VBMicrolensing 5.6.1 supplies the extended-source and binary magnifications. The KMT baseline-capped pass is a calibration sensitivity test. Four binary template shapes at the grid stage and local multistart binary optimization on selected events give alternative fits. Their parameter coverage cannot exclude an arbitrary distant stellar companion. Peak-window and full-baseline chi-square values have different error normalizations and should only be compared within the same analysis.

