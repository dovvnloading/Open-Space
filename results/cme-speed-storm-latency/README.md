# A source-pinned 1996–2025 analysis of CME speed and storm-peak latency

**Status: new Open-Space exploratory archival association, independently replayed; not an operational forecasting model or claim of discovery priority.**

This work adds a previously untouched catalog-level heliophysics population to Open-Space. The upstream storm table and article were published by [D'Awarave and Vemareddy (2026), *Solar Physics* 301, 85](https://doi.org/10.1007/s11207-026-02666-6). They assembled intense storms and solar sources; our new analysis asks a *different* question: among severe ICME-related storms, is the **delay from shock detection at L1 to the Dst minimum** associated with an earlier solar remote-sensing observable, the CME plane-of-sky speed? The authors' existing findings and earlier shock/sheath literature are acknowledged below; do not claim that CME speed or storm timing were unknown.

The source was pinned to upstream commit `97747daa2ab3256986bdeafa707dd8214da71911`, rather than an unversioned live feed. The CSV SHA-256 is `8794251f244af391ab275ae6e10d5ca36c37c591461ca7dd8d3d4d31d2de0eb2`. All results and sampling decisions are in [`summary.json`](summary.json); exact code is [`examples/cme_speed_dst_delay.py`](../../examples/cme_speed_dst_delay.py). [GitHub Actions replay #43770574081](https://github.com/dovvnloading/Open-Space/actions/runs/43770574081) completed successfully, using a fresh download of this pinned input.

## Quantitative finding

From **139** archived severe Dst events, exclude **19 CIR-driven** records without an ICME shock timestamp; of the remaining **120** shock-associated events, **113** have valid positive CME speeds.

| Sample | Severe-storm cases | CME speed/latency Pearson r | Median shock→Dst minimum with CME <700 km/s | Median with CME ≥1000 km/s | Slow minus fast |
|---|---:|---:|---:|---:|---:|
| **All, 1996–2025** | **113** | **−0.379** | **22.75 h** (51 cases) | **10.93 h** (40 cases) | **11.82 h** |
| Earlier era, 1996–2010 | 68 | −0.310 | 21.70 h | 10.88 h | 10.82 h |
| Later era, 2011–2025 | 45 | −0.497 | 25.60 h | 12.77 h | 12.83 h |
| 72-hour-separated shocks only | 104 | −0.385 | 22.75 h | 10.67 h | 12.08 h |
| Exclude unusual speeds or >48h delays | 109 | −0.429 | 22.75 h | 10.72 h | 12.03 h |

Bootstrap (4,000 event replicates, seed 20261115): 95% interval of **−0.511 to −0.241** for Pearson r, and **+4.67 to +15.78 hours** for the slow-minus-fast median gap. A **calendar-year block bootstrap**, protecting against repeated storms in the same year, gives **−0.484 to −0.277** for r.

A separate **10,000-permutation null** shuffles Dst delays only *within five-year periods* to preserve coarse solar-cycle composition. Its 95% r range is **−0.202 to +0.137**; the observed −0.379 is outside it (two-sided empirical **p = 0.0001**). Partial correlation after removing linear associations with the final storm Dst, peak Bz, quiet-region flag, and five-year period is **−0.281**. **Those adjustment variables are retrospective; they are never presented as predictors available before a storm.**

This association is a meaningful cross-era *timing effect*, not a demonstration that CME speed can predict whether a storm will become severe.

## Interpretation, falsification and limitations

1. **Source selection:** The sample was assembled *because* storms reached Dst around/below −100 nT. No nonstorm CME denominator is present, so calculating probability of any geomagnetic storm from this table would be a selection-bias error.
2. **Retrospective source association:** Even if CME sky-plane speed is measurable before Earth arrival, selecting the actual CME responsible for an observed storm in this curated catalog uses hindsight. The timing result is not operational forecast skill.
3. **Known physics:** Shock/sheath and ejecta-driven Dst delays have been studied for decades. A shorter delay for fast CMEs could reflect earlier strong southward Bz or sheath dynamics, geometry and the distribution of solar sources; this analysis does **not** isolate a mechanism or prove causality.
4. **Exploratory cutpoints:** The <700 and ≥1000 km/s groups were chosen after inspecting archive values. Their impressive difference is descriptive. The permutation hypothesis is the unbinned speed/latency correlation. The earlier/later eras are robustness checks **not blind held-out prospective confirmation**.
5. **Uncertainties:** L1 shock and hour-averaged Dst timestamps have subhour-to-hour uncertainty; shocks/ejecta can interact; extreme events can cause overlapping geomagnetic storm phases. A 72-hour shock-gap sensitivity retains direction and magnitude but does not remove all event dependence. Median CIs are event-bootstrap intervals, not full observational systematic errors.
6. **Missingness and scope:** Seven ICME records lack CME speed. Source-specific selection and missing speeds can bias estimates; no extrapolation to CIR storms or all CMEs is justified.

## Reproduce

```sh
python -m unittest discover -s tests -p test_cme_speed_dst_delay.py -v
python examples/cme_speed_dst_delay.py --fetch \
    --bootstrap 4000 --permutations 10000 --output /tmp/cme-latency.json
```

The script needs only Python standard library, downloads 15 KB from an **immutable, pinned GitHub commit**, records SHA-256, and writes machine-readable results. Source raw CSV is **not** vendored.

## Prior literature

- [D'Awarave & Vemareddy (2026), 139 severe storms and solar source attribution](https://doi.org/10.1007/s11207-026-02666-6), also [institutional abstract](https://prints.iiap.res.in/handle/2248/8995?show=full).
- [Zhang et al. (2007), solar and interplanetary sources during 1996–2005](https://doi.org/10.1029/2007JA012321).
- [Lugaz et al. (2017), timing and geoeffectiveness of shocks/sheaths](https://pmc.ncbi.nlm.nih.gov/articles/PMC5882492/).
- [Annotated independent severe-storm list, NASA CDAW](https://cdaw.gsfc.nasa.gov/CME_list/dst100/), useful for a future source-definition comparison.

**Next scientifically decisive test:** preregister fixed CME-speed cutpoints on the *independent* NASA CDAW severe-storm event set or a genuinely later era, validate shock and source identifiers, and obtain an unselected denominator of Earth-directed CMEs. Until then, the discovery is a reproducible, sizable archival association rather than new physics or forecasting capability.
