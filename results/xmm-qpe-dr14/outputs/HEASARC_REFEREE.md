# NASA HEASARC QPE frozen-archive referee audit

**Attempted:** 670/670 observations. **Verified:** 637 observations, **910** source EPIC-pn 0.2–2 keV light curves of 967 eligible detection rows.
**Preliminary flags:** 6; **validated novel QPEs: 0**; no [REMARKABLE].
**Independent source-masked field particle-curve comparisons:** 910. **Scaled off-source aperture repeat flags:** 1.
**Poisson count-domain proxy recoveries:** 836/1820. Not a spatial/event-list injection, independent population FPR, or validated sensitivity.
**Missing verified pn curves:** 33 observation IDs / 57 detection rows. No missing product is a negative.

## Batch receipts

| Batch receipt | Attempted | Status |
|---|---:|---|
| [heasarc_batch_cohort_000_080.json](heasarc_batch_cohort_000_080.json) | 80/80 | BATCH_FINISHED_UNVETTED |
| [heasarc_batch_cohort_080_200.json](heasarc_batch_cohort_080_200.json) | 120/120 | BATCH_FINISHED_UNVETTED |
| [heasarc_batch_cohort_200_320.json](heasarc_batch_cohort_200_320.json) | 120/120 | BATCH_FINISHED_UNVETTED |
| [heasarc_batch_cohort_320_495.json](heasarc_batch_cohort_320_495.json) | 175/175 | BATCH_FINISHED_UNVETTED |
| [heasarc_batch_cohort_495_670.json](heasarc_batch_cohort_495_670.json) | 175/175 | BATCH_FINISHED_UNVETTED |

## Trigger-by-trigger assessment

### 0144570101 / 101445701010002
FITS-position-verified at (206.90164, 17.567958); recurrence [1.833] hours, field particle-flare coincidences 1.
**Disposition:** KNOWN_SEYFERT2_AGN_AND_PARTICLE_FLARE_NOT_NOVEL_QPE. Gaia DR3 counterpart rows: 1; significant astrometric star flag: False.
SAS pile-up/soft off-source event-list re-extraction, full SIMBAD/4LAC and eRASS1/TNS/ATel are outstanding unless separately documented. QPE vs AGN/TDE remains unconfirmed for any unknown flags.

- Published study: https://arxiv.org/abs/0801.4246
- Published study: https://arxiv.org/abs/1803.10976
- Published study: https://doi.org/10.1051/0004-6361/202142444
### 0841482101 / 108414821010004
FITS-position-verified at (205.435623, -0.975584); recurrence [5.583, 9.25] hours, field particle-flare coincidences 1.
**Disposition:** UNVETTED_FLAG_NOT_QPE_CANDIDATE. Gaia DR3 counterpart rows: 1; significant astrometric star flag: False.
SAS pile-up/soft off-source event-list re-extraction, full SIMBAD/4LAC and eRASS1/TNS/ATel are outstanding unless separately documented. QPE vs AGN/TDE remains unconfirmed for any unknown flags.

### 0851180501 / 108511805010001
FITS-position-verified at (195.500647, 27.782786); recurrence [5.5, 3.75] hours, field particle-flare coincidences 1.
**Disposition:** UNVETTED_FLAG_NOT_QPE_CANDIDATE. Gaia DR3 counterpart rows: 1; significant astrometric star flag: False.
SAS pile-up/soft off-source event-list re-extraction, full SIMBAD/4LAC and eRASS1/TNS/ATel are outstanding unless separately documented. QPE vs AGN/TDE remains unconfirmed for any unknown flags.

### 0864560101 / 108645601010001
FITS-position-verified at (195.500813, 27.782451); recurrence [1.333, 5.0, 2.917, 5.417, 3.25, 5.583, 3.917] hours, field particle-flare coincidences 1.
**Disposition:** UNVETTED_FLAG_NOT_QPE_CANDIDATE. Gaia DR3 counterpart rows: 0; significant astrometric star flag: False.
SAS pile-up/soft off-source event-list re-extraction, full SIMBAD/4LAC and eRASS1/TNS/ATel are outstanding unless separately documented. QPE vs AGN/TDE remains unconfirmed for any unknown flags.

### 0904640201 / 109046402010001
FITS-position-verified at (195.500454, 27.782836); recurrence [1.667, 4.917, 1.417, 4.583, 1.583, 4.583, 1.417, 4.833, 2.667] hours, field particle-flare coincidences 1.
**Disposition:** UNVETTED_FLAG_NOT_QPE_CANDIDATE. Gaia DR3 counterpart rows: 1; significant astrometric star flag: False.
SAS pile-up/soft off-source event-list re-extraction, full SIMBAD/4LAC and eRASS1/TNS/ATel are outstanding unless separately documented. QPE vs AGN/TDE remains unconfirmed for any unknown flags.

### 0904640301 / 109046403010002
FITS-position-verified at (195.50059, 27.782729); recurrence [4.917, 1.5, 4.667, 4.417, 6.083, 3.333] hours, field particle-flare coincidences 1.
**Disposition:** UNVETTED_FLAG_NOT_QPE_CANDIDATE. Gaia DR3 counterpart rows: 1; significant astrometric star flag: False.
SAS pile-up/soft off-source event-list re-extraction, full SIMBAD/4LAC and eRASS1/TNS/ATel are outstanding unless separately documented. QPE vs AGN/TDE remains unconfirmed for any unknown flags.


## Referee decision

**KEEP DRAFT / DO NOT MERGE.** Verified HEASARC pn source curves are observational evidence; unverified source rows are not QPE nulls. A source-masked particle-flare curve and scaled source BACK1–3 are useful controls, but not a full independent soft-band event-region background check. Count-domain injections do not provide trial-adjusted false-alarm probability or event-list completeness. No validated novel QPE and no remarkable entry.
