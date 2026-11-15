# Independent archival searches

Starting revision: `f85370fcd39f596c439bf58675ef0eb3d3880af0`. The main result directory, open pull request #31, and active branch names were checked before selecting these questions. The concurrent Cepheid, pulsar recurrence, ultracompact companion and brown-dwarf alignment studies are outside the selection.

| Question | Observing unit | What would defeat the claim? |
| --- | --- | --- |
| Did an independently recorded hard gamma-ray transient accompany an archival fast X-ray transient? | One X-ray localization and event window compared with an independently timed detector stream | A time convention mismatch, an incompatible direction, high accidental-match rate, or an existing identification |
| Are there close approaches between numbered asteroids with a resolvable gravitational signal? | A pair of small bodies and a stated time interval | A miss on numerical ephemerides, uncertainty larger than the effect, or prior publication of the same event |
| Can archival X-ray data establish periodic eclipses in optically periodic binaries? | One optically localized binary observed over one or more X-ray visits | The fold is due to source blending, observing gaps, background flares, chance alignment or a published measurement of the same effect |

Only an analysis with source identifiers, replayable inputs, measured evidence and a literature check should be promoted to a discovery candidate. Distinguish an observed event, a proposed physical association and a prediction for future observations.

## Completed investigations

- [FXT 981217 and the untriggered BATSE burst](../results/fxt-981217-batse/README.md) is a historical recovery and catalog discrepancy. The 2007 table cites an earlier joint analysis. The packet includes the original prior-art row, joint detector light curves, time-shift controls and a third-spacecraft coverage test.
- [M33N 1996-12a and a 544-day Mira](../results/m33n-1996-12a/README.md) compares the six historical nova positions with the later Mira catalogue, modern photometry and shifted positions. The 1997 co-identification was already indexed in SIMBAD.
- [Comet 407P and asteroid (75547)](../results/comet-407p-asteroid-75547/README.md) compares the three-dimensional encounter in two orbit services with the ZTF images around closest approach.

## X-ray binary screen

A separate screen linked 9,006 bright 5XMM-DR15 stacked-source rows (`ep_cts > 750`, `sum_flag <= 2`, `gaia_match_prob > 0.9`) to Gaia DR3 eclipsing-binary candidates by the 5XMM Gaia source ID. It returned 122 X-ray rows for 121 Gaia IDs. One striking flare in ObsID 0692810301, source 002, coincides with Kiso J201946.11+385133.0. It rises from about 0.03 to 0.28 EPIC-pn counts/s in 1,500-second bins. Pasquato et al. had already visually labeled this exact observation and source a flare (catalogue row 602, classifier probability 0.9913). EXOD DR1 had separately detected the event. We closed this lead. Another source, J021721.6−034129, showed a broad variation during a visit spanning 1.5 optical cycles. Its rate stayed high during the primary optical eclipse, and the low phase changed substantially between cycles. We found insufficient evidence for a new periodic X-ray eclipse. [5XMM-DR15](https://xmmssc.aip.de/cms/catalogues/5xmm-dr15/) · [Gaia DR3 eclipsing binaries](https://doi.org/10.1051/0004-6361/202245330) · [Published flare catalogue](https://vizier.cfa.harvard.edu/viz-bin/VizieR?-source=J/A%2BA/708/A224) · [EXOD catalogue](https://vizier.cfa.harvard.edu/viz-bin/VizieR?-source=J/A%2BA/697/A50)
