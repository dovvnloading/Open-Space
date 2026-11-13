# Which stellar signals were already reported?

Eleven of 120 exploration targets passed the two-sector sinusoidal coherence rule. All eleven have previous variability evidence. The strongest apparent new rotation signal illustrates a trap. TIC 33881153 has no VSX or Gaia DR3 variable within 42 arcsec and no TESS-SVC entry. A later, much larger rotation survey had already measured its period over 27 sectors. Open-Space recovered the first harmonic. [Catalog measurements and provenance](../results/stellar-catalogs.json) · [Boyle et al. TESS All-Sky Rotation Survey, accepted release](https://zenodo.org/records/19917941)

## Catalog comparison

Positions are taken from the discovery FITS headers. The search queried the October 2026 VSX VizieR snapshot and both the Gaia DR3 variability summary and classifier tables. It also checked the exact TIC ID in the 84,046-object TESS Stellar Variability Catalog (SVC) and 1,046,317-object TESS All-Sky Rotation Survey (TARS). All 36 cone queries and 24 exact-ID queries completed across the eleven exploration objects and the separate control. A dash marks a successful empty exact-TIC lookup. Periods are in days. [AAVSO VSX](https://cdsarc.cds.unistra.fr/viz-bin/cat/B/vsx) · [Gaia DR3 variability](https://cdsarc.cds.unistra.fr/viz-bin/cat/I/358) · [MAST TESS-SVC](https://archive.stsci.edu/hlsp/tess-svc) · [MAST TARS](https://archive.stsci.edu/hlsp/tars)

| TIC | Measured | SVC | TARS | Interpretation |
| ---: | ---: | ---: | ---: | --- |
| 31780319 | 1.32242 | 1.32408 | 1.3231 | Period reproduced. TARS flags a possible binary. |
| 33881153 | 4.78536 | — | 9.6079 | Half the published period within 0.4%. |
| 50380257 | 0.67439 | — | 1.3485 | Half the independent MAST eclipsing-binary catalog period of 1.348285 d. |
| 55527074 | 2.28048 | 2.29568 | — | Same timescale as SVC. |
| 234338669 | 0.459124 | 0.459153 | 0.4595 | Both surveys report this modulation. |
| 234495456 | 0.380722 | 0.380574 | 0.3806 | Both surveys report this modulation. |
| 238194430 | 0.699785 | 0.700267 | 0.7001 | Both surveys report this modulation. |
| 279511712 | 1.65648 | 1.64967 | — | V0450 Car, an SPB star in VSX at 0.33 arcsec, listed at 1.6518 d. |
| 382256692 | 0.148970 | — | — | BB Dor, a nova-like cataclysmic variable, seen at its published superhump period. |
| 425975034 | 0.371215 | — | 0.3719 | TARS flags two contaminating sources. The VSX eclipsing binary BB Tuc lies 36.26 arcsec away. The signal's host requires localization. |
| 426029868 | 4.68480 | 5.60519 | — | The later sector reproduces the SVC period. The early value has negligible power in that sector. |
| 100100827, control | 0.470403 | 0.470567 | — | WASP-18. VSX lists the known 0.941453 d orbital period and Gaia's classifier reports EP. |

For BB Dor, VSX lists an orbital period of 0.154095 d. The shorter optical period is a separate, previously reported signal. Bruch measured 0.1489422 d in TESS data, just 2.4 seconds from this run's 0.1489703 d. Finding the true literature match required distinguishing a photometric superhump from an orbital period. TESS-SVC, TARS, Bruch's analysis and this pipeline can use overlapping TESS observations. Their agreement demonstrates reproducibility across reductions or methods. [Bruch, MNRAS 519, 352, section on BB Dor](https://academic.oup.com/mnras/article/519/1/352/6855266) · [Spectroscopic orbital measurement](https://academic.oup.com/mnras/article/422/3/2332/1045205) · [TESS eclipsing-binary matches](../results/tess-eb-check.json)

## Later-sector and harmonic checks

The discovery grid fits periods up to 5 days. In the audit we widened it to 0.125–12.5 d, ten samples per independent frequency element, a floating mean and the same 30-minute PDCSAP medians. On the exact discovery period and its double we recorded the fraction of variance explained. Short sector baselines limit period resolution, and spot patterns can evolve.

| TIC | Sector | Role in discovery | Strongest period (d) | Variance at original period | Variance at double period |
| ---: | ---: | --- | ---: | ---: | ---: |
| 33881153 | 1 | Discovery | 10.21 | 19.1% | 57.8% |
| 33881153 | 2 | Two-sector validation | 9.81 | 25.0% | 41.1% |
| 33881153 | 3 | Later observation | 8.73 | 12.7% | 59.1% |
| 33881153 | 28 | Later observation | 5.52 | 11.9% | 13.6% |
| 33881153 | 29 | Later observation | 4.65 | 32.5% | 1.3% |

TARS reports 9.6079 ± 0.1031 d for TIC 33881153 and flags multiple periods, a possible binary, and a doubled adopted period. Sectors 1–3 support the longer timescale. In sector 29 the shorter harmonic dominates. Four of 32 validation-sector block shuffles passed the initial discovery rule. The TESS pixel follow-up was weak or unstable in sectors 1 and 2. Taken together, the data support previously cataloged variability with changing waveform. [Five-sector measurements with FITS hashes](../results/stellar-harmonic-33881153.json) · [Stellar pixel measurements](../results/pixels-stellar.json)

TIC 426029868 produced a weaker 4.68480 d signal, with semiamplitudes of 0.356 and 0.232 parts per thousand. TESS-SVC reports 5.60519 ± 0.4681 d in sector 13. Our sector 13 reanalysis peaks at 5.59397 d, explaining 14.34% of variance. The frozen 4.68480 d explains 0.068%. This later sector reproduces the published SVC timescale. Since SVC used the same sector, it supplies independent observing dates relative to sectors 1–2 and a cross-check of the published reduction. [Three-sector measurements with FITS hashes](../results/stellar-harmonic-426029868.json)

## Running the checks

The callable supports both displaced pixel centroids and targets with TIC IDs. Supply the TIC ID only for an identified target. SVC and TARS use exact IDs. VSX and Gaia report angular separations for all catalog matches within the supplied radius. VSX uses catalog J2000 coordinates, Gaia uses epoch 2016.0, and separations use the supplied position without proper-motion propagation. A positional match inside two nominal TESS pixels is an association to inspect. Period comparisons report raw ratios and the nearest factor among 0.5, 1, and 2. They apply no automatic physical classification.

```python
from pathlib import Path
from open_space.variable_catalog import search

result = search(66.29706643, -76.19357287, Path("cache"), radius_arcsec=42,
                tic_id=33881153, period_days=4.7853551153771825)
```

The first exact-ID check downloads the two public releases, verifies their SHA-256 and row counts, and builds shared SQLite indexes. Their compressed or uncompressed archive files together occupy about 240 MB, and the indexes add about 130 MB. Successful cone responses retain their exact URL and response digest. A new cache queries the current VizieR release, whose timestamp can change even for identical scientific rows. A cached run reuses the original bytes. Inspect the individual statuses and the aggregate status before treating a query as empty.

```bash
python -m open_space.variable_catalog results/stellar-source-cohort.json --cache cache --output stellar-catalogs.json
python -m open_space.stellar_harmonics --tic 33881153 --period 4.7853551153771825 \
  --sectors 1 2 3 28 29 --discovery-sectors 1 2 --catalog-dir cache/inventories \
  --cache cache/fits --output harmonic-33881153.json
python -m open_space.stellar_harmonics --tic 426029868 --period 4.684800735743462 \
  --sectors 1 2 13 --discovery-sectors 1 2 --catalog-dir cache/inventories \
  --cache cache/fits --output harmonic-426029868.json
```

The [complete 121-target input](../results/stellar-source-cohort.json) is a JSON object with a `targets` list from the two-sector stellar search. The worker also uses this lookup after a stellar signal crosses its screening threshold. The survey catalogs publish specific monitored samples with different selection criteria. A clean catalog cone is insufficient evidence for an unreported source.
