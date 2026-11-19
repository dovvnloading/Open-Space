# Coverage, selection denominators and unresolved experiments

**Unit**: one 2023-02-13 KM3NeT event, MJD(UTC)=59988.0533299, ICRS RA=94.3°, Dec=-7.8°, 99% circle r=3.0°; event angular posterior approximated by published circular summaries, not native likelihood. All bounds in [README.md](README.md) frozen before screening.

## completed as of initial independent readback (2026-11-19 UTC)

| Parent catalogue and time gate | Parent counts and status | Coincidences | Physical chance probability |
|---|---|---|---|
| Literature PoS/2026 AJ 17 blazar candidates | All 17 coordinates independently reprojected. 2 inside 50%, 3 inside 68%, 8 inside 90%, 17 inside 99%; maximum table rounding discrepancy 0.041394° | 17 **KNOWN**, not new | **NOT ESTIMABLE** from position-selected list. #8 PMN J0606-0724 p=0.0026 is authors' **pre-trial** radio flare probability; global p missing |
| TNS public objects daily 2026-11-19 whole-sky ZIP (source SHA256 b0adf81c18415b3f352949cbe315837f759bd226d13c5a33ea323e96ed339f5a), discovery ±30d | 213,510 valid positions/dates; 2,938 epoch-eligible rows | **0** in 3° circle in successful second attempt | p upper-tail for ≥0 is 1.000 (trivial); **p lower-tail for ≤0 is 0.399 under 999 control shifts** despite a random-sky mean of 3.128 and all-sky uniform prior 2.013. This demonstrates highly nonuniform sky coverage; *no source probability* |
| TNS entries explicitly classified TDE, 30–180d (additional to 30d sample) | 6 epoch-eligible rows | **0** within 3° | p upper-tail for ≥0 is 1.000; no physical association p |
| CHIME/FRB Catalog 2 v20251231 (source SHA256 5108ada779d279a2547d9f9e73ae25bfdd40d8496d6ba7255ec29c6629057a48), ±1d | 5,045 sub-burst rows, 5,008 usable positions and times; 0 within ±1d (nearest catalogued burst 1.8409d from event) | **0** in 3° | p for count ≥0 is 1.000; no 2023 Feb 13 coverage/exposure test; no physical chance p |

**Final independent runner:** [successful checksum and independent spherical/999-shift replay](https://github.com/dovvnloading/Open-Space/actions/runs/44443328806); 2,944 frozen temporal parent positions reproduce all three family counts and sky-shift distributions. **Source receipts and caveat:** [GitHub Actions evidence for the successful TNS+CHIME pass](https://github.com/dovvnloading/Open-Space/actions/runs/44442696625/artifacts/14492385686) (second run). First TNS attempt returned HTTP 403 and is separately receipted; the second succeeded with HTTP 200. This is a **catalogue scan**, not a completeness-calibrated nondetection. TNS discovery dates can differ by weeks from onset/peak; no forced ZTF/ATLAS fluxes were retrieved; CHIME transient sky exposure is time-dependent. A 3° circular counterpart membership does not replace the exact KM3NeT or CHIME localization posterior.

## Required searches left pending: do not classify as no-match

| Preregistered family | Exact gate | State / reason |
|---|---|---|
| NASA Swift GRB table (HTTP 200, source SHA256 5f78bc4c0bef0a5c68c64b9d57e359b0ec7d681167df59e6ae90075333dd5a07), 3,900,078 HTML bytes | ±1d and ±1h | **COMPLETED Swift subset**: 1,658 tabulated events with date and BAT direction (2,046 listed total, some incomplete); **0 bursts in ±1 day**; spatial p for count ≥0 is 1.000 but all-mission rate/exposure is untested; [audit artifact](https://github.com/dovvnloading/Open-Space/actions/runs/44442949469/artifacts/14492201221) |
| Non-Swift GCN/GBM GRBs | ±1d and ±1h | PENDING complete GCN circular + Fermi GBM all-sky sample, revised regions and coverage |
| Swift MAXI/BAT/XRT non-GRB transients | ±7d primary, ±30d sensitivity | PENDING time-dependent exposure and source-episode/cross-instrument dedup |
| ZTF/ATLAS forced photometry | ±30d primary, TDE/nuclear ±180d | PENDING full field footprints, forced-source light curves and pre-event non-detections |
| New Fermi-LAT blazar light-curve repository | ±30d ±180d | PENDING differential weekly bins from frozen latest public release and eligible flux thresholds; gamma flares in original papers are KNOWN |
| OVRO 15 GHz and VLBA | ±30d ±180d | PENDING independently selected full monitored sample and available dates, not the 17 handpicked blazars; PMN J0606-0724 flare is KNOWN |
| VLASS / ASKAP RACS epoch changes | same source cross-epoch; ±180d contemporaneity check | PENDING matched beams, astrometry, parent source counts, field coverage; Filipović et al.'s 1,052 ASKAP sources + three named sources are KNOWN |
| NEOWISE W1/W2 epoch-flare | ±180d; cadence audit | PENDING native epoch fluxes and seasonal detection limits; the collaboration already checked WISE/NEOWISE candidate blazars |
| Galaxy clusters / nearby star-formers and persistent galaxies | positional 3° | PENDING full parent catalogue, distance/redshift, X-ray/radio and dust biases; no short-time chance model applies |
| 8-family global null and physical neutrino model | family maximum over shifts and trials | PENDING all complete parent families, source selection/exposure masks and independent population prior; **no discovery-level association p can be given** |

Uniform source densities are diagnostic priors: TNS transient sky distribution and angular exposure are very non-uniform near the Galactic plane. The 999 sky shifts hold event declination fixed, omit offsets within 6.5°, and preserve the actual listed source coordinates/epochs. This *does not* preserve optical survey completeness or CHIME observing cadence. Source-specific 0 p-values are **never inferred** from absence, failed HTTP, or a finite Monte Carlo.

## Falsifiable next checks

1. Retrieve the native sky exposure mask for TNS discovery feeds (ZTF/ATLAS limiting magnitude, temporal coverage), and obtain light curves for all *named objects* inside 3° regardless of discovery date; selection on discovery epoch alone misses earlier eruptions.
2. Crossmatch GBM/GCN high-energy transients using 90% native localizations; an unpublished GBM event could be missed by the Swift GRB subset.
3. Use fixed, complete LAT and OVRO parent source lists and a family-wise maximum statistic for simultaneous flares. Demand post-trial result, not the reported 0.26% pre-trial coincidence.
4. Reconstruct annual VLASS/RACS radio baselines and NEOWISE scan-coadd time sampling, including random positions with comparable Galactic latitude/pointing, before accepting any orphan radio/IR flare.
5. Test a neutrino-emission model for any survivor against IceCube 15-year steady/flaring null (https://arxiv.org/abs/2609.10931) and multiphase EM constraints.

**Disposition for this initial dataset:** LIMITED / INCONCLUSIVE. No counterpart identification. No [REMARKABLE] label or log entry. Positive results in the listed missing families are not ruled out.
