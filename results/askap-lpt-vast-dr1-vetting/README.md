# VAST Extragalactic DR1: native-image vetting of 60 strong radio-variable leads

**Result — no new long-period radio transient is validated.** Follow-up to merged [PR #200](https://github.com/dovvnloading/Open-Space/pull/200): of **60** ranked leads (**693** epoch measurements), **15** coincide with known ATNF pulsars, **2** have plausible stellar counterparts without proven radio flares, **6** have ambiguous Gaia positional matches, and **37** remain unclassified because essential polarisation and host constraints are missing. The on/off subset of **14** splits into **3 pulsars, 1 plausible star and 10 unresolved**. This is **measured source-cutout screening, not a discovery or a calibrated null limit**.

## Native images and quality checks

All **60/60 original source FITS cubes** were retrieved from the public [CSIRO VAST Extragalactic DR1 release](https://doi.org/10.25919/nh9d-t846); the **693/693** image extensions have exact DATE-OBS matches to the published measurements and valid WCS. The source-ID lookup comes from the [survey authors' public download dictionary](https://github.com/IrisdeRuiter/VAST_Extragalactic_DR1), verified against pinned Git blob `b0dc867d6baf4a967d84fba781d2b27ed96bbe67`. Its remotely provided pickle is **read as strings only**, never executed or deserialized. Every downloaded native cube has a SHA-256 receipt.

For each epoch the analysis estimates median/MAD local background and noise around the source, central peak significance, neighboring bright peaks, negative troughs and pixel distance from the cutout edge. It records **122** possible bright-neighbor/sidelobe flags across **17** sources, including nine deep-negative flags, and produces **60** inspection montages. These are automated warning signs, **not proof of instrumental artefacts** or independent human inspection. The cutouts are centered, so they **cannot check proximity to the original ASKAP field edge or off-stamp bright-source sidelobes**; those would require full native-field images.

## Independent survey checks

| Public data | Completed evidence | Interpretation limit |
|---|---|---|
| RACS low/mid/high Stokes I | 46 / 47 / 42 positional catalogue matches; **45** two-or-more-band spectral slopes | Noncontemporaneous frequencies/observations can mimic spectral indices; blending is possible |
| RACS low/mid/high Stokes V; VAST Stokes V | **Not measured**; RACS-low2 public DAP source-list inventory did not expose an accessible small V catalogue | No circular-polarisation detection *or upper limit* is claimed; separately staged CASDA images and leakage validation are required |
| Gaia DR3 | 12 matches within 3 arcsec; 2 plausible non-pulsar stellar associations, 6 ambiguous close positional matches | A Gaia position alone is **not** a flare-star classification; strong parallax or very bright optical counterpart strengthens the association |
| Legacy Surveys DR10 (NOIRLab TAP) | Queries attempted for all 60; **no usable returned table** | All 60 are unavailable, **not** 60 optical non-detections |
| unWISE (IRSA) | 33 matches, 26 catalogue nonmatches, one unavailable | Positional matches only, no secure optical/IR host ID |
| Pan-STARRS DR1 (MAST) | 30 position rows; **13** with at least two detections and a valid PSF-band magnitude | Other position rows are not photometric detections; coverage and depth matter |
| eRASS1 | 8 source positions within 12 arcsec | Associations are provisional; nonmatches are not X-ray flux upper limits |
| Swift 2SXPS | 60 queried circles, no catalogue entries | No exposure map or point-source sensitivity constraints; not a confirmed absence of X-rays |

The [RACS-low2 polarisation release](https://doi.org/10.1017/pasa.2026.10225) describes Stokes V observations and leakage cuts, but **the VAST DR1 FITS source cubes used here are Stokes I only**, and the required calibrated V images were not retrieved. Consequently no compact, strongly polarised and reliably counterpart-free survivor has been established. There is no available minute-to-hour cadence time series here: a periodicity test and candidate-specific literature/ATel search are deferred pending a confirmed target with public sub-epoch data. An unexplained source is **not** automatically a new long-period radio transient.

## Reproducibility and evidence

The [single-job GitHub Actions workflow](../../.github/workflows/askap-lpt-vast-dr1-vetting.yml) runs on pushes to the research branch, path-filtered to this result's Python sources; generated `output/**` never retriggers it. All remote access and computation run **in that one job** on GitHub's networked runner. The initial collection downloads source cubes and queries catalogues; subsequent runs can use the committed, receipt-backed evidence for deterministic independent reclassification. Evidence was committed in bounded table and montage batches.

- [Classifications of all 60 sources](output/classifications.csv), with class basis, image flags, spectral slopes, Gaia parallax strength and direct counterpart catalogue query status.
- [Every epoch's image diagnostics](output/epoch_image_audit.csv) (693 records); [60 per-source montages](output/montages/).
- [Source-cube and index hashes](output/image_receipts.json), [individual remote catalogue queries](output/counterpart_receipts.json), [RACS-low2 Stokes-V data-availability probe](output/racs_low2_v_receipts.json).
- [Summary and status counts](output/summary.json), [independent machine referee](output/referee.json), [conservative Gaia/photometry reclassification receipt](output/reclassification_receipts.json).
- [Independent scientific referee report](REFEREE.md) separating verified evidence from untested claims.

**Labelling:** this is a bounded, reproducible follow-up with no validated discovery-level event. No exceptional-finding log entry is warranted.
