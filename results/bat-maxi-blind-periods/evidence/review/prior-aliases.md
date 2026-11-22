# Published orbital periods and one-day aliases

For each daily-survey period, the proposed carrier is the published orbital frequency or its second harmonic. Once-per-day sampling folds this carrier onto the lower frequency $|1 - f_{\rm carrier}|$. The orbit-product fits use their native sampling, quoted flux errors and piecewise means. The daily-frequency and carrier gains reduce the same weighted chi-square, with identical local frequency refinement. A positive difference favors the carrier in that interval. These are descriptive fits to the same underlying monitor observations used in selection.

Offsets use the epoch's Rayleigh width, 1 divided by its span. The frequency predicted from the printed orbital period can evolve across the baseline. An offset larger than one resolution element is reported explicitly. The output JSON includes conditional gains and sampling correlations when the carrier or its orbital fundamental was the stored winner. Near-unit correlations make the corresponding sinusoidal components difficult to distinguish.

## H 1538−522

The published orbit is 3.72831 d. Harmonic 2 has a period of 1.864155 d and predicts a daily reflection at 2.1571998 d. [Hemphill et al. (2019)](https://ui.adsabs.harvard.edu/abs/2019ApJ...873...62H/abstract)

The measured daily period is 2.1572002 d. The reflected carrier has a larger gain in 6 of 6 available comparisons. Its fitted maximum is within one Rayleigh width of the printed orbital prediction in 6 comparisons.

| Monitor | Epoch | Daily gain | Carrier period (d) | Carrier gain | Carrier − daily | Orbital gain | Carrier offset (Rayleigh) |
|---|---|---:|---:|---:|---:|---:|---:|
| BAT | full | 17.029 | 1.8642647 | 1303.2 | +1286.2 | 4316.6 | 0.25 |
| BAT | early | 63.95 | 1.8641547 | 2860.9 | +2797 | 4492.8 | 0.00 |
| BAT | late | 1.2226 | 1.8641547 | 84.602 | +83.379 | 661.55 | 0.00 |
| MAXI | full | 17.447 | 1.8641547 | 3325.8 | +3308.4 | 9940.8 | 0.00 |
| MAXI | early | 9.0389 | 1.8641547 | 1947 | +1937.9 | 5000.1 | 0.00 |
| MAXI | late | 13.827 | 1.8641547 | 1392.7 | +1378.9 | 4824 | 0.00 |

BAT favors the orbital fundamental at 3.7278694 d over its reflected harmonic and the daily period. The daily/fundamental sampling correlations are 0.615 and 0.162. Adding the orbital fundamental to the daily fit gains 4473.8, and adding the daily term to the orbital fit gains 174.21. The specific harmonic's conditional gains were not stored.
The daily-data nuisance-fit status is “requires interpretation” (frozen grid frequency fails after conditional nuisance projection, stronger fitted waveform supplies at least half the observed candidate amplitude).

MAXI favors the orbital fundamental at 3.7283094 d over its reflected harmonic and the daily period. The daily/fundamental sampling correlations are 0.007 and 0.002. Adding the orbital fundamental to the daily fit gains 9942.8, and adding the daily term to the orbital fit gains 19.484. The specific harmonic's conditional gains were not stored.
The daily-data nuisance-fit status is “requires interpretation” (frozen grid frequency fails after conditional nuisance projection, stronger fitted waveform supplies at least half the observed candidate amplitude).

## SMC X-1

The published orbit is 3.8919232 d. Harmonic 2 has a period of 1.9459616 d and predicts a daily reflection at 2.0571254 d. [Falanga et al. (2015)](https://ui.adsabs.harvard.edu/abs/2015A%26A...577A.130F/abstract)

The measured daily period is 2.0572321 d. The reflected carrier has a larger gain in 6 of 6 available comparisons. Its fitted maximum is within one Rayleigh width of the printed orbital prediction in 6 comparisons.

| Monitor | Epoch | Daily gain | Carrier period (d) | Carrier gain | Carrier − daily | Orbital gain | Carrier offset (Rayleigh) |
|---|---|---:|---:|---:|---:|---:|---:|
| BAT | full | 18.279 | 1.945986 | 3092.4 | +3074.1 | 14970 | 0.05 |
| BAT | early | 45.815 | 1.9458661 | 5381.8 | +5336 | 16292 | 0.10 |
| BAT | late | 8.3232 | 1.9458661 | 422.24 | +413.92 | 1919.8 | 0.10 |
| MAXI | full | 41.292 | 1.9458661 | 6669.3 | +6628 | 16564 | 0.16 |
| MAXI | early | 36.652 | 1.9458661 | 3767.9 | +3731.2 | 9022.3 | 0.08 |
| MAXI | late | 24.194 | 1.9458661 | 2925.2 | +2901 | 7598.4 | 0.08 |

BAT favors the orbital fundamental at 3.8922118 d over its reflected harmonic and the daily period. The daily/fundamental sampling correlations are 0.544 and 0.251. Adding the orbital fundamental to the daily fit gains 15872, and adding the daily term to the orbital fit gains 920.59. The specific harmonic's conditional gains were not stored.
The daily-data nuisance-fit status is “requires interpretation” (additional significant stronger families remain beyond the audit limit, frozen grid frequency fails after conditional nuisance projection, stronger fitted waveform supplies at least half the observed candidate amplitude).

MAXI favors the orbital fundamental at 3.8917322 d over its reflected harmonic and the daily period. The daily/fundamental sampling correlations are 0.004 and 0.001. Adding the orbital fundamental to the daily fit gains 16565, and adding the daily term to the orbital fit gains 42.457. The specific harmonic's conditional gains were not stored.
The daily-data nuisance-fit status is “requires interpretation” (frozen grid frequency fails after conditional nuisance projection, stronger fitted waveform supplies at least half the observed candidate amplitude).

## Her X-1

The published orbit is 1.7001676 d. Harmonic 1 has a period of 1.7001676 d and predicts a daily reflection at 2.4282295 d. [Staubert, Klochkov & Wilms (2009)](https://arxiv.org/abs/0904.2307)

The measured daily period is 2.4282312 d. The reflected carrier has a larger gain in 6 of 6 available comparisons. Its fitted maximum is within one Rayleigh width of the printed orbital prediction in 6 comparisons.

| Monitor | Epoch | Daily gain | Carrier period (d) | Carrier gain | Carrier − daily | Orbital gain | Carrier offset (Rayleigh) |
|---|---|---:|---:|---:|---:|---:|---:|
| BAT | full | 271.5 | 1.7001668 | 20177 | +19905 | 20177 | 0.00 |
| BAT | early | 252.53 | 1.7001668 | 35844 | +35592 | 35844 | 0.00 |
| BAT | late | 89.351 | 1.7001668 | 2272.5 | +2183.2 | 2272.5 | 0.00 |
| MAXI | full | 0.71576 | 1.7001668 | 7807.5 | +7806.8 | 7807.5 | 0.00 |
| MAXI | early | 2.6644 | 1.7001668 | 4003.5 | +4000.9 | 4003.5 | 0.00 |
| MAXI | late | 3.9773 | 1.7001668 | 3820.5 | +3816.5 | 3820.5 | 0.00 |

BAT full-span target/carrier sampling correlations are 0.709 and 0.496. Adding the carrier to the daily fit gains 25322, and adding the daily term to the carrier fit gains 5416.7.
The daily-data nuisance-fit status is “retained at fixed grid frequency”.

MAXI full-span target/carrier sampling correlations are 0.012 and 0.000. Adding the carrier to the daily fit gains 7808.7, and adding the daily term to the carrier fit gains 1.9216.
The daily-data nuisance-fit status is “retained at fixed grid frequency”.

## LMC X-4

The published orbit is 1.408379 d. Harmonic 1 has a period of 1.408379 d and predicts a daily reflection at 3.4487057 d. [Falanga et al. (2015)](https://ui.adsabs.harvard.edu/abs/2015A%26A...577A.130F/abstract)

The measured daily period is 3.4488233 d. The reflected carrier has a larger gain in 6 of 6 available comparisons. Its fitted maximum is within one Rayleigh width of the printed orbital prediction in 6 comparisons.

| Monitor | Epoch | Daily gain | Carrier period (d) | Carrier gain | Carrier − daily | Orbital gain | Carrier offset (Rayleigh) |
|---|---|---:|---:|---:|---:|---:|---:|
| BAT | full | 30.783 | 1.4083594 | 3034.1 | +3003.3 | 3034.1 | 0.08 |
| BAT | early | 30.321 | 1.4083594 | 4948.8 | +4918.5 | 4948.8 | 0.04 |
| BAT | late | 10.327 | 1.4083594 | 288.11 | +277.78 | 288.11 | 0.04 |
| MAXI | full | 5.3853 | 1.4083594 | 1188.4 | +1183.1 | 1188.4 | 0.06 |
| MAXI | early | 13.776 | 1.4083594 | 697.93 | +684.16 | 697.93 | 0.03 |
| MAXI | late | 7.1719 | 1.4083594 | 507.19 | +500.02 | 507.19 | 0.03 |

BAT full-span target/carrier sampling correlations are 0.322 and 0.220. Adding the carrier to the daily fit gains 3194.5, and adding the daily term to the carrier fit gains 191.17.
The daily-data nuisance-fit status is “retained at fixed grid frequency”.

MAXI full-span target/carrier sampling correlations are 0.006 and 0.006. Adding the carrier to the daily fit gains 1188.6, and adding the daily term to the carrier fit gains 5.5806.
The daily-data nuisance-fit status is “retained at fixed grid frequency”.
