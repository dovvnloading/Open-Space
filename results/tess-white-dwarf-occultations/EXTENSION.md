# Independent bright-catalogue extension

Freeze the full 157-row bright parent. Remove the nine Sector 109 detector sources, the positive controls, nonunique or position-only TIC matches, and sources lacking a finite hydrogen-model mass of 0.1–1.44 solar masses and radius of 0.002–0.2 solar radii. Derive the radius from catalogue mass and log g. For each remaining source choose its newest publicly obtainable Sector 55–108 TESScut visit at observer-epoch coordinates. Query coverage at 2025.0, check each listed visit from newest first at its sector midpoint, and verify target position, detector, and 200-second cadence. If the service reports an error, log the attempted visit and try an earlier listed one. After a 200-second cube is delivered, evaluate the same Gaia crowding cut, unique TIC check, and quality-zero aperture extraction used in Sector 109. Do not replace a usable observation because of crowding, a variable light curve or its search outcome.

Use one visit per independent target. Apply the Sector 109 detector and the same four fixed period/radius cells, random phase and random conditional transit inclination, nominal Gaia aperture fraction and catalogue stellar model. Run 50 injections per target per cell with seed 337109 independent of the initial cohort. Count every original compatible periodic trigger, of known or unknown origin, in its applicable cell. All candidates get an independent assessment. Report finite-Monte-Carlo one-sided confidence as in the original cohort. Preserve this cohort separately. No discovery controls or repeated stars enter the denominator. Scope conditional on the models and the selected bright archive-covered sources.

Cycle 5 began in Sector 56. The frozen availability range also admitted Sector 55 to the metadata checks. Its 10-minute FFIs cannot pass the 200-second cadence gate. Replay requests start at Sector 56 and the independent audit verifies that every included visit lies in S56–108.

## Measured sensitivity

Of the 140 requested stars, 39 passed, 89 failed the Gaia crowding cut, four had nonpositive median aperture flux, and eight had no available qualifying visit. The accepted visits yielded 417,350 finite quality-zero cadences. The detector used 417,173 after excluding short observing segments. Each cell has zero original period-compatible periodic triggers among these 39 stars.

| Period (d) | Radius (Earth radii) | Recovered / injected | Primary 95% transiting-host ceiling | Secondary 95% transiting-host ceiling |
|---:|---:|---:|---:|---:|
| 0.25 | 1 | 1,338 / 1,950 | 19.37% | 11.77% |
| 0.25 | 11.2 | 1,898 / 1,950 | 11.11% | 8.16% |
| 1.40794 | 1 | 1,449 / 1,950 | 17.42% | 10.82% |
| 1.40794 | 11.2 | 1,896 / 1,950 | 11.13% | 8.17% |

The primary calculation uses simultaneous per-star Clopper–Pearson lower bounds on detection efficiency with a total error probability of 0.025. It inverts the resulting Poisson-binomial survey count at 0.025. Each upper bound has at least 95% pointwise coverage under the assumed common independent host rate, stellar and dilution models, and phase and inclination distributions. At all inclinations the primary calculation allows 100% occurrence in every cell. [Fixed-cell trial ledger](extension/extension_recovery_trials.csv.gz) · [Primary bounds](extension/extension_limits.json)

The separately reported zero-count calculation bounds the *sum* of efficiencies with a one-sided Hoeffding inequality (error probability 0.005). For a host fraction $f$ and a lower bound $L$ on that sum, the arithmetic–geometric mean inequality gives $P(0)\le(1-fL/39)^{39}$. Inverting at 0.045 gives the secondary ceilings. Weighting each star by its geometric transit probability gives a 98.37% isotropic-host ceiling at 0.25 days and 11.2 Earth radii. In the other three cells this procedure allows 100% occurrence. The 6-hour model assumes a detached circular orbit and does not model tidal survival. Both sets of bounds apply to the selected bright stars and the stipulated individual cells. [Secondary calculation](aggregate_sensitivity.py) · [Secondary bounds](extension/weighted_bounds.json)

## Reproduction checks

The frozen-record audits can run without archive access. They check cohort membership, all 9,000 fixed-cell trials across both cohorts, the compatible-host counts in the extension and both confidence calculations.

```bash
python audit_recovery.py sector109
python audit_extension.py extension
```

A second pixel reduction reacquired the first, middle and last included extension stars sorted by numeric Gaia ID. For TIC 5636651 (S85), TIC 405122007 (S82) and TIC 76774929 (S92), the FITS TIME, FLUX and QUALITY arrays and the complete extracted light-curve files reproduced their original SHA-256 digests. The blind unmodified searches and all 750 individual trial outcomes agreed, including the recovered ephemerides for matched continuous-period trials. Replaying the same trials tests reproducibility and adds zero independent injections to the occurrence calculation. The receipts retain the newly generated and original FITS file hashes. [Replay script](replay_sensitivity.py) · [MAST receipts and replay counts](extension/pixel_replay.json)

```bash
OPENBLAS_NUM_THREADS=1 python replay_sensitivity.py replay-cache replay-report.json
```
