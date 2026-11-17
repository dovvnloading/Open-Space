# Repeat-epoch compact-source variability in Abell 2744 and MACS0647

No mutually identified compact source doubled or halved in flux at a signed difference significance of at least 7σ in the three matched-filter comparisons below. The population consists of sources detected separately in both epochs. A transient that vanished entirely, appeared for the first time, was deblended differently, or lay on a resolved host can fall outside this test.

| Field | UT observing dates | Filter | Mutual matches | Compact at ≥5σ in both epochs | ≥20σ compact calibration controls | Median formal 5σ AB magnitude at eligible positions, early / late |
|---|---|---|---:|---:|---:|---|
| Abell 2744 | 2024 July 3 / October 23 | F200W | 223 | 101 | 82 | 29.16 / 28.09 |
| Abell 2744 | 2024 July 3 / October 23 | F444W | 172 | 152 | 136 | 28.52 / 27.74 |
| MACS0647 | 2022 September 23 / 2023 January 8 | F200W | 681 | 447 | 418 | 28.44 / 28.43 |

The October Abell 2744 pair uses the JWST pipeline's observation association c1021. The two visits compared here sample only their common catalog coverage. Counts refer to filter entries. Abell 2744 sources can appear in both bands.

The calculation uses aperture-corrected fluxes and the pipeline's extended-source flag. Mutual nearest neighbors within 0.5″ establish a median relative translation. A second mutual match uses 0.12″. Median late/early flux ratios of bright compact controls set relative gains. Robust scatter of their log-ratios was 3.55%, 3.17%, and 2.71%, respectively. The difference uncertainty adds the two formal flux errors and this fractional scatter times the lower epoch flux in quadrature. The screen requires compact morphology and positive ≥5σ photometry in each epoch. The AB magnitudes are formal flux-uncertainty diagnostics for eligible positions. No injection or blank-field area calculation was performed on these catalogs.

Run `python catfetch.py`, then `python catmeasure.py` with NumPy, SciPy and pandas. The six exact STScI archive products, byte counts and SHA-256 hashes are in [catalog_inputs.csv](catalog_inputs.csv). [Per-comparison results](catalog_screen.json) · [Threshold crossings](catalog_flags.csv)

Prior Abell 2744 work includes [TNS AstroNote 2025-217](https://www.wis-tns.org/astronotes/astronote/2025-217) and [PLATE-23a](https://arxiv.org/abs/2608.20781). MACS0647 contains [previously reported lensed stars at z = 4.758](https://academic.oup.com/mnrasl/article/527/1/L7/7281004). There are no flags to cross-identify from this test.
