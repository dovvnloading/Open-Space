# Mars gravity field: 2026 independent cross-release audit (work in progress)

**Status:** Research checkpoint; no scientific discovery claimed at this stage. Created 2026-11-15.

## Question

Does the September 24, 2026 public release of JPL's MRO120G Martian gravity field reveal spatially coherent, operationally consequential changes relative to MRO120F (2020) and an independent GSFC solution? Which changes survive normalization, reference-frame, regularization, uncertainty and correlated-fit controls?

**Observing unit:** degree/order spherical-harmonic coefficients of Mars' static gravity potential (not independent raw radio measurements). This is a *model-revision comparison*, not new independent spacecraft observations.

**Predeclared falsification:** reject any "new gravity anomaly" interpretation if changes reflect different model constraints, reference orientation, GM, tide conventions, or correlated inversion assumptions; require independent field or tracking-data evidence before asserting physical change.

## Provenance

- NASA PDS Geosciences: MRO derived gravity [release notes](https://pds-geosciences.wustl.edu/mro/mro-m-rss-5-sdp-v1/mrors_1xxx/errata.txt), last updated 2026-09-24.
- JPL model MRO120G (2026), [table](https://pds-geosciences.wustl.edu/mro/mro-m-rss-5-sdp-v1/mrors_1xxx/data/shadr/jgmro_120g_sha.tab) and [label](https://pds-geosciences.wustl.edu/mro/mro-m-rss-5-sdp-v1/mrors_1xxx/data/shadr/jgmro_120g_sha.lbl).
- JPL model MRO120F (2020), [table](https://pds-geosciences.wustl.edu/mro/mro-m-rss-5-sdp-v1/mrors_1xxx/data/shadr/jgmro_120f_sha.tab) and [label](https://pds-geosciences.wustl.edu/mro/mro-m-rss-5-sdp-v1/mrors_1xxx/data/shadr/jgmro_120f_sha.lbl).
- Independent GSFC GMM3 [table](https://pds-geosciences.wustl.edu/mro/mro-m-rss-5-sdp-v1/mrors_1xxx/data/shadr/gmm3_120_sha.tab).
- Primary paper: Konopliv et al. (2026), *Updated Mars Orientation and Gravity Field and Indication of Decadal Changes in Rotation Rate*, DOI [10.3847/PSJ/ae9ba7](https://doi.org/10.3847/PSJ/ae9ba7). Do not rediscover effects already published there.

## Nonoverlap

At branch creation the repository had zero open PRs. The existing `results/` directory covers Mars impacts, solar-system occultations, Martian crater relationships and other topics, but contains no MRO120G gravity-field analysis. This work creates only new files; it will not modify existing investigations.

## Planned deliverables

1. Reproducible downloader/SHADR parser with SHA-256 receipts and label checks.
2. Degree-band and spatial sensitivity audit with independent baseline, explicit uncertainty caveats and controlled tests.
3. Only retain claims materially supported by actual archive data, with a clear counterexample/next test.
