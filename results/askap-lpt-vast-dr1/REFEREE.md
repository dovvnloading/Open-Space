# Independent referee: VAST Extragalactic DR1 long-period-transient screen

**Decision: accept as a reproducible *catalogue-level screening result*, NOT as an LPT discovery, candidate confirmation or a calibrated null. No [REMARKABLE] marker.**

## Independently checked

- The frozen [primary-input receipts](output/network_receipts.json) give image/source/measurement SHA-256s with three successful 200 downloads. [Referee assertions](output/referee.json) independently parse CSV/JSONL and reconstruct all cuts for **all 149** source rows (not just 60 leads). The byte/manifest/source-image joins, source uniqueness, rank order and recheck against 693 epoch rows passed. The published release advertises 6,449,202 measurement rows, whereas the downloaded Parquet metadata contains **6,449,203**; the one-row discrepancy was deliberately **not** corrected or upgraded to a discovery.
- The release's **2,945** dated and uniquely identified image/field records were frozen in [frozen_images.jsonl](output/frozen_images.jsonl) before source screening. Six predeclared known LPT coordinates have **zero** intersections with the approximate field-square envelopes. Therefore no eligible positive LPT recovery or completeness/sensitivity test is possible from this dataset alone.
- Independently regrouping the 693 ranked-epoch CSV rows yields **14** light-curve patterns with at least two unflagged >=8-sigma Selavy detections and at least one forced <3-sigma measurement. Four other sources have a single pair of observations less than one hour apart (minimum **13.105 minutes**); there is no measured period, folded time series or repeated minute-to-hour pulse train.
- The latest [PSRCAT v18](https://doi.org/10.25919/qkn0-vm25) download, SHA-256 `bac32471fef3e60cfd734f1db85503a9bd0a565018cb98edff54475f87eb642f`, supplied **4,262** parseable pulsar positions. A crossmatch of all 60 top-ranked sources finds **15** within **2 arcseconds**; three of the fourteen apparent on/off patterns are **known pulsars**: `13564869` / PSR J2222-0137, `14593764` / PSR J2144-3933, `13593032` / PSR J0034-0534. These are false positives **for novelty**, not necessarily spurious radio signals. There are **11** other ATNF-negative on/off patterns, unvetted for star/AGN/other known radio-source classes.
- Gaia DR3 5-arcsec cone queries succeeded for **12/12** sampled top-ranked sources, with positional matches near radio IDs `14199978` and `14192287`; photometric and proper-motion-aware classification was **not** attempted, and the other 48 leaders were not Gaia-queried.

## Conservative surviving pattern ledger (not discoveries)

| Source ID | Good >=8σ Selavy epochs | Forced <3σ epochs | Brightest good peak (mJy/beam) | ATNF within 2 arcsec |
| --- | ---: | ---: | ---: | --- |
| `13564869` | 7 | 1 | 20.467 | J2222-0137 |
| `14133427` | 3 | 2 | 13.995 | none identified |
| `13805206` | 4 | 4 | 12.765 | none identified |
| `13299183` | 9 | 1 | 9.548 | none identified |
| `13399649` | 5 | 1 | 7.541 | none identified |
| `15304967` | 3 | 4 | 7.019 | none identified |
| `14593764` | 4 | 3 | 7.004 | J2144-3933 |
| `13767816` | 3 | 2 | 6.263 | none identified |
| `13593032` | 4 | 1 | 6.236 | J0034-0534 |
| `14015818` | 8 | 1 | 6.231 | none identified |
| `14164377` | 2 | 3 | 6.162 | none identified |
| `13629007` | 2 | 5 | 6.111 | none identified |
| `13497959` | 2 | 5 | 6.046 | none identified |
| `14476792` | 9 | 1 | 5.849 | none identified |

## Adversarial limitations

The source catalogue was already pre-filtered for isolated, detected, relatively compact sources and is not complete for intermittent LPTs. The top-60 ranking is a severity/effort cap, **not** a random or held-out sample. The variance and max/min cuts can select known pulsars, scintillating AGN, solar-/RFI-contaminated snapshots and flare stars. The DR1 fluxes are separated by weeks to months and cannot deliver a ten-minute-to-hours periodicity proof. Native FITS cutouts, per-image rms checks, full RACS I/Q/U/V polarisation and spectral indices, Legacy/unWISE deep counterpart exclusion, multi-survey pulsar/literature vetos, Swift/eROSITA X-ray searches, cadence-dependent recovery simulations and a measured scintillation/flare-star false-positive distribution are **not completed**. The contextual **84/117** AGN+star classification in the published sample is not an estimated selection FPR. With 3 known pulsars among 14 apparent on/off patterns, the **minimum known-object contamination fraction is 21.4%**, but that too is *not* a population FPR.

**No discovery, no all-sky null, no false-positive-calibrated search sensitivity.** The research record is ready to merge only with this bounded disposition, full receipts, and no [REMARKABLE] tag. Scientific follow-up requires a separate investigation and a true minute/hour cadence.
