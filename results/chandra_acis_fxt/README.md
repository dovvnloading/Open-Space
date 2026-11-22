# Chandra ACIS fast X-ray transient search — bounded null result

**Null result: no new extragalactic fast X-ray transient was confirmed in the audited, deposited subset.** Real fast flares and unresolved signals were recorded; neither the date-eligible full archive nor the intrinsic transient population has been ruled out.

## What was actually measured

The archive census contains **6,608 date-eligible ACIS imaging ObsIDs** first made public from **2021-01-01 to 2026-11-21**, inclusive. The stopped analysis includes **6,259 ObsIDs with complete event-file scanning**, every prescribed time scramble and injection draw, and **100,245,154.25 seconds (100.245154 Ms)** of eligible 60-second flare-start time. Of these, **5,892** completed the required catalog and photometry checks; **367** have incomplete follow-up. The included event searches produced **33,802 raw real trigger rows**, **341 distinct instrument-stage passes**, but **zero confirmed *new extragalactic* events**.

| Trial length | Measured eligible start time |
|---:|---:|
| 60 s | 100.245154 Ms |
| 180 s | 98.022579 Ms |
| 600 s | 90.266042 Ms |
| 1,800 s | 69.520332 Ms |
| 3,600 s | 42.617302 Ms |

**Not searched for this published result:** 349 of the 6,608 dated ObsIDs; by group: **72 (34), 181 (35), 184 (35), 185 (35), 186 (35), 187 (35), 188 (35), 189 (35), 190 (35), 191 (35)**. Group **72** (34 ObsIDs) completed its original science step but its original GitHub evidence upload failed; an independent repeat was audited separately and *not* silently substituted into the original inference. Nine other groups, **181 and 184–191** (315 ObsIDs), completed original runners and have archived artifacts, but were not included in the deposited measurement table. They are therefore explicitly **not searched/included for this bounded result**, not claimed as absent astronomical data. The individual **349 ObsIDs** appear under `not_searched_for_bounded_result.groups` in [the stop ledger](STOP_LEDGER.json). Of 206 primary Actions jobs, 204 succeeded; the two failures were group 72's artifact transport and the final report job. No additional search was launched to close these gaps.

The 14 searched observations with unverified release dates form a **separate**, non-date-eligible supplement: 369.4 ks of eligible start time, 175 real trigger rows and zero retained events. They are not added to the 6,259 or 100.245154 Ms totals. [Original source run](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443) · [Per-ObsID measured exposure and byte receipts](evidence/primary_ingest/44985879443-1/PROGRESS.md) · [Archived source-group audit](evidence/primary_ingest/44985879443-1/coverage.json)

## Ten strongest follow-up cases at the stop decision

Only these ten previously flagged, high-priority or held-out cases were adjudicated for final publication using already deposited evidence; the remaining instrument-stage passes were **not** promoted to discoveries. *Unresolved* means the signal or emitter is not conclusively classified, **not** that it is a new extragalactic transient.

| Original ObsID and trigger | Stop decision | Evidence and limitation |
|---|---|---|
| 30892 `real/0` | Unresolved | Highest prespecified held-out seed score; 0/24 scrambled exceedances, but two-test-adjusted FAP is 0.08 and no secure host, redshift, or independent extragalactic association exists. |
| 23037 `real/0` | Excluded | CCD 8 row/frame-correlated events dominate the apparent burst; original held-out adjusted FAP is 0.16. Not a credible astrophysical transient. |
| 23473 `real/0` | Excluded | 9/24 held-out controls reach its seed score; adjusted FAP is 0.80. No defensible new-source claim. |
| 27089 `real/0` | Unresolved | Ten dispersed trial photons, a substantial annular count and quiet flanks; no host in checked catalog cones, no independent significant repeat, and no support for association with GW170817. |
| 22558 `real/0` | Unresolved | Compact multiframe flare and gradual tail; G=18.06 Gaia positional overlap without secure distance makes a stellar origin possible. |
| 22929 `real/0` | Unresolved | Eight photons in eight frames and quiet flanks; 13 covered independent visits show no trial-corrected repeat, and no host or redshift was established. |
| 27939 `real/0` | Unresolved | Thirty-five-photon sky-localized flare, but red optical objects and Gaia candidates in a crowded Galactic WR 112 field allow an ordinary stellar origin. |
| 32068 `real/0` | Unresolved | Forty-one-photon flare in the Galactic Lighthouse PWN field; nearby Gaia star with parallax and incomplete optical screening prevent extragalactic attribution. |
| 24264 `real/0` | Unresolved | Flare in a crowded LMC field has five ambiguous Gaia overlaps and no reliable source membership or distance. |
| 28479 `real/0` | Unresolved | Sharp 180-second dither-consistent rise, but persistent diffuse Perseus emission prevents confirmation of an independent point source or galaxy host. |

The fixed prospective held-out analysis assigned 3,351 ObsIDs, found 3,345 evaluable, and retained three real triggers versus 22 total survivors in its 24 matched control catalogs. The best held-out event, **ObsID 30892**, has **0 of 24 control exceedances**, but the protocol's two-test-adjusted FAP is **0.08**; its host and distance are unestablished. The other two, ObsIDs 23037 and 23473, have adjusted values **0.16** and **0.80**. None passes the combined significance, detector, and host-association requirements for a newly identified extragalactic event. [Preregistered rules](HELD_OUT.md) · [Original held-out evaluation](evidence/held_out/README.md) · [Source identity and counts](evidence/held_out/survivors.json)

**One-line extragalactic verdict: none of the ten establishes a plausible *new extragalactic* fast X-ray transient after all available checks.** The prior detailed work on the reproducible 78-photon ObsID 23324 flare and the targeted 37-photon ObsID 32380 flare remains in [priority-source notes](PRIORITY.md): the former had already been identified as variable in that same Chandra visit, and the latter favors stellar activity.

## Negative controls, completeness and limits

The fixed moderate-count 128-ObsID cohort provided 1.5912 Ms of eligible 60-second start time; two of 18 real trigger rows passed instrumental tests, and each was equalled or exceeded in one of 24 time scrambles (permutation estimate **0.08**). Injections recovered **939 / 1,297** attempted synthetic events (72.4%) **in that cohort and its injected signal distribution**, not over the whole archive. The separate undated-release sample recovered **139 / 165** injections. Each deposited observation retains its source URL/checksum, full raw detections, rejected rows, controls and injection attempts. [Cohort and controls](evidence/provisional/moderate-count-v1/44985879443-1/summary.json) · [Methods and reproducibility](METHODS.md)

The [machine stop ledger](STOP_LEDGER.json) preserves **all 341** original instrument-stage passes by rule: **10 individually reviewed**, of which **2 were excluded** and **8 unresolved**; the other **331** retain final-stop status **unresolved**, including events with earlier provisional dispositions. Thus **339 are unresolved at stop**. Their underlying observation ledgers, earlier vetting decisions, candidate diagnostics and exact event-file receipts remain available for independent review. The stop decision did **not** reclassify them as astrophysical non-events.

**Scope limits:** 349 dated ObsIDs are outside the measured sample, 367 included ObsIDs have incomplete catalog/photometry checks, the release date is unresolved for 14 additional fields, and no complete population-wide sensitivity function has been calibrated. The trial-count and injection controls do **not** justify an intrinsic volumetric rate bound. This is a bounded **null result for confirmed new extragalactic sources**, not evidence that fast X-ray flares do not exist.

[Strong-source field notes](PRIMARY_EVEN_NOTES.md) · [Odd-shard review and source identities](PRIMARY_ODD_NOTES.md) · [Original files and checksums](evidence/primary_ingest/44985879443-1/progress_inputs.csv.gz)
