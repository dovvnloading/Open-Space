# No significant neutrino–transient association in the frozen samples

Thirteen archival WISE flare–neutrino pairs involve **10 distinct IceCube alerts**. Independently randomizing alert right ascension gave a median of **10** matched alerts (95th percentile 15). Permuting alert times within year, quarter and tier gave a median of **12** (95th percentile 13). The count's unadjusted upper-tail probabilities are **0.6061** and **1.0000** respectively. Both joint-search, trial-corrected upper tails are **1.0000**. The smallest 90% area among the matching alerts is also consistent with chance. These coincidences provide no evidence of a new neutrino source.

For 2024 through **2026-11-17 12:09:02 UTC**, the native maps contained **zero registered optical, X-ray or gamma-ray pairs**. The sole pair admitted by the older GCN coordinate bounds was **IC240424A / GRB 240419B**, 4.560 days apart. IceCube explicitly reported the burst inside its 90% region in the original circular. That recovers a published positional coincidence. The public coordinate extrema alone cannot reproduce the region in scrambled skies. No source-level significance is assigned. [IceCube GCN 36283](https://gcn.nasa.gov/circulars/36283) · [Burst and afterglow audit](optx/grb240419b/README.md)

## Scope and counting

The IceCat-1 V4 table and native likelihood files contain 348 alert-grade events. Eight carry cosmic-ray vetoes. Of the 340 eligible tracks, seven have public IceCat-2 replacements. The public IceCat-2 release is a seven-event preliminary release. GCN supplied 71 distinct 2024–26 Gold/Bronze notices after collapsing revisions. Two were retracted. The latest alert in the frozen inventory is November 13, 2026. Of 69 active GCN tracks, 23 have released multi-order probability maps and 46 have only revised coordinate extrema in the archive. There are **409 distinct eligible tracks**, 363 with native regions. The historical table ends on 2023-10-14 and the separate GCN selection starts on 2024-01-01. This design leaves the intervening weeks unsampled. [IceCube receipts and maps](ice/README.md)

| Search | Frozen eligible EM epochs | Alert–epoch pairs passing time cut and having native geometry | Pairs inside native 90% regions | Pairs inside legacy bounds requiring separate verification |
|---|---:|---:|---:|---:|
| 2024–26 optical TDEs and nuclear flares | 73 | 613 | 0 | 0 |
| 2024–26 EP / Swift XRT and BAT | 362 | 66 | 0 | 1, already reported by IceCube |
| 2024–26 Fermi-LAT weekly flare peaks | 232 | 74 | 0 | 0 |
| 2024–26 primary FLAIRES sample | 0 | 0 | — | — |
| Earlier LAT weekly flare peaks | 1,229 | 1,047 | 0 | 0 |
| Earlier WISE FLAIRES peaks | 822 | 22,033 | 13 | 0 |

The modern calculation retains nine late-2023 archive events for possible negative lags, with complete relevant quarter strata. All 23 modern maps with optical/X-ray/gamma time support are from 2026. The independent bounds audit covers the 46 other active GCN alerts, including 2024–25. A position outside the publisher's coordinate extrema is excluded; a position inside them needs the actual contour or an explicit publisher membership statement. The historical comparison uses all 340 active archive maps, with 23 later native maps also loaded to preserve the cohort and any boundary overlap. Counts exclude verified control aliases and host positions. [Modern counts](analysis/modern/coverage.csv) · [Archival counts](analysis/archival/coverage.csv) · [All archival pairs](analysis/archival/native_pairs.csv) · [Modern GCN bounds audit](analysis/modern/legacy_unresolved.csv)

## Frozen samples and tests

- Optical inputs combine the TNS daily mirror, ZTF BTS, TDECat, Rubin/ZTF LIONS and a published Fink nuclear-flare selection. The 73 sources include 69 classified TDEs and four additional nuclear flares. Survey discovery and peak dates trace different physical epochs. The public catalogues cannot establish completeness for unclassified nuclear transients. [Inputs and selection](optx/README.md)
- The X-ray sample includes 192 EP-WXT detections, 140 independent Swift BAT triggers and 34 serendipitous Swift XRT detections. Four identified common events count once. The EP sample comes from public GCN reports. Brightness, pointing, and reporting determine which events enter these samples. [X-ray event list](optx/xray_events.csv)
- The Fermi LCR freeze covers 1,415 variable DR2 blazars or blazar candidates and 1,325,903 reported weekly bins. Of these, 612 sources had at least ten bins meeting the quality cuts. The flare threshold used the median of **qualifying detections**, which can suppress sensitivity to faint states. At best the input light curves reach the 2026-11-06 weekly midpoint. [LCR selection](gammair/README.md) · [Offline audit of all fits](gammair/FIT_AUDIT.md)
- The corrected FLAIRES catalogue has 823 peaks dated 2010–22. One named control host was masked. Its release date in 2026 adds no contemporary MIR coverage. Of 17 modern 2MRS supplementary peaks, 16 remain after masking the controls; none produced a pair. A separate post-hoc TDECat check found seven 2024 W2 brightest visits, six at the end of their observing series. The one bracketed maximum, AT2022fpx, passes the coordinate-extrema screen for IC240307A (14.25° from the best fit), with native membership undetermined. It is outside the statistical sample. [MIR selection](gammair/README.md) · [Supplemental check](analysis/modern/fpx.md)

The [first branch commit](https://github.com/dovvnloading/Open-Space/commit/09feaebda0171b7589c74e1fd4e66ecde177818a) fixes the modern search. The [bounds audit](BOUNDS_AUDIT.md) and [historical extension](ARCHIVAL_PLAN.md) were registered before inspecting the relevant discovery pairs. The latter followed the first modern screen and the historical control replay. Both primary nulls use 9,999 trials with inclusive ties and a plus-one tail estimate. The four modern families and two statistics cost eight tests. A claim across both epochs costs twelve. For archival MIR localization specificity, the unadjusted probabilities are 0.2015 (RA) and 1.0000 (time). Its largest inverse area is 0.3998 deg⁻², corresponding to 2.50 deg². The initial MC cannot resolve a discovery-level tail. [Archival statistics and exceedances](analysis/archival/statistics.json) · [Modern statistics](analysis/modern/statistics.json)

The modern time permutations produced invariant count and specificity statistics for all three exposed families. The archival MIR count varies, but its specificity statistic is invariant under time shuffling. A ±180-day MIR window is wider than a shuffle quarter. These trials show limited temporal discrimination. The RA test holds the catalogue positions and epochs, alert declinations and timestamps, and localization shapes fixed. A three-harmonic RA uniformity check gives diagnostic p=0.982 (340 archive events) and p=0.470 (23 recent maps). It cannot establish exact exposure or survey completeness. No flux, rate, or diffuse-emission upper limit is derived here. [Null histograms](analysis/archival/null_histograms.csv) · [Exposure diagnostics](analysis/exposure.json)

## Positive controls

| Counterpart | Original reported selection | Native IceCat-1 V4 | Native preliminary IceCat-2 |
|---|---|---|---|
| TXS 0506+056 / IC170922A | Published high state recovered. Closest qualifying LCR peak is 7.371 d before the alert, outside the registered ±7 d window. | Inside | Inside |
| AT2019dsg / IC191001A | Recovered by the authors' original selection | Outside | Outside |
| AT2019fdr / IC200530A | Recovered by the authors' original selection | Outside | Outside |
| AT2019aalc / IC191119A | Recovered by the authors' original selection | Inside | Outside |

The original optical procedure reproduces exactly three matches over 43 alerts × 63 flares. The newer localizations change their spatial interpretation. The TXS weekly bin coincident with the alert is 2.91 times its qualifying-bin median. Keeping the registered threefold threshold and ±7-day window documents a real loss in sensitivity to a known flare. [Historical replay](controls/README.md) · [Numerical containment](ice/native_control_membership.csv) · [IceCube reconstruction comparison](https://pos.sissa.it/501/1224/pdf)

## Check status and replay

Completed work includes freezing the inputs, deduplicating alerts and named episodes, native-region membership, the independent bounds audit, both nulls, historical controls, and source-by-source follow-up. The source populations and finite interval are defined in this directory. Fractional source attribution, detector livetime modelling, a new reconstruction of the 46 older GCN contours and completeness for unpublished or unclassified transients remain open. NEOWISE's observing mission ended in 2024. [Final NEOWISE data release](https://irsa.ipac.caltech.edu/data/WISE/docs/release/NEOWISE/)

Run from this directory with Python 3.12, NumPy, Astropy, astropy-healpix and Matplotlib. The primary crossmatch is offline. The replay copies the frozen tree to a new disposable directory, runs the core analyses there and checks the frozen hashes in both trees. It verifies deterministic reproduction with the same code and inputs, not independent survey selection, native-map reconstruction or physical validation. The snapshot-specific READMEs document how to regenerate each survey projection. `SHA256SUMS` fingerprints the committed evidence and scripts.

```sh
sha256sum -c SHA256SUMS
python analysis/replay_copy.py /path/to/new/disposable-replay
python controls/historical.py controls/inputs
```
