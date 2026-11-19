# Dark 1LHAASO / >100-TeV sources: public counterpart census

**Claim (calibrated catalogue screen; no new accelerator identification).** Of the **90** unique 1LHAASO sources in a frozen, hashed **2023-preprint transcription of Cao et al.** (43 marked >100 TeV), **25** lack an association in the paper's TeV-associated column (**10** are UHE). Newer public survey catalogues account for multiple geometrical overlaps, but the data **do not support a new PeVatron class, microquasar or old-pulsar halo identification**. A unique accelerator and a complete multiwavelength negative census remain unestablished.

## Checked findings — with the matched-latitude null

| Counterpart catalogue/class | Among 25 originally TeV-unassociated | Mean in 512 shifted skies | One-sided count probability | What it says |
|---|---:|---:|---:|---|
| Fermi **4FGL-DR4 v35**, GeV | **19** | 9.46 | **0.00195** (null resolution floor) | An interesting *positional* excess, **not** 19 identifications |
| ATNF pulsars, latest public tar | **16** | 10.15 | 0.158 | No established pulsar-population excess |
| **Green 2024 October** Galactic SNRs (310 entries) | **9** | 5.38 | 0.207 | No established remnant-class excess |
| **4HWC (2026)**, 85 TeV sources | **13** | Withheld | Withheld | Prior TeV counterpart candidates; **declination mask approximate** |
| **eRASS1 DR1 hard**, 5,466 X-ray detections | **0 of just 1** approximately covered | Withheld | Withheld | The other 24 sky positions are **not X-ray non-detections** |
| H.E.S.S. HGPS **2018** legacy, 78 entries | **0 of 4** approximately covered | Withheld | Withheld | Legacy veto only; not newest H.E.S.S. |
| MeerKAT SMGPS, ASKAP/EMU, Gaia young clusters, SNRcat full, eRASS1 **soft** | **Pending** | — | — | Cannot assert anything from absent tables |

**Five** of the 25 lack a geometrical match in *all three completely ingested GeV/ATNF/Green classes*, including **four UHE sources**: 1LHAASO J0007+5659u, J0206+4302u, J1959+1129u and J2200+5643u. None of those five matches the 4HWC positional aperture, but **they are not globally confirmed dark**: multiple radio, soft X-ray and cluster catalogues remain untested, and proposed counterparts already appear in the literature (notably the 2026 Peanut pulsar proposal and the original LMXB suggestion).

**Pulsar energy budget:** 55 pair-level ATNF coincidences; only **14** have a usable KM2A 25–100 TeV power-law flux, reported distance, and spin-down power. Of those, **seven** satisfy the necessary Lγ/Ė ≤ 1 energy screen, **seven** exceed it; **41** pairs have insufficient inputs. Energetic possibility is **not** an identity, age/transport model or confirmation.

## Frozen inputs, receipts and scripts

The **2024 ApJS published CDS/VizieR table could not be retrieved** during the November 2026 Strasbourg/CDS outage. The frozen primary input is instead the independent Gammapy transcription of the **2023 Cao preprint** (180 detector-component rows; 90 deduplicated sources), **not represented as an exact 2024-table copy**.

- **Primary bytes:** [frozen/1lhaaso_gammapy_2023_source.csv](output/frozen/1lhaaso_gammapy_2023_source.csv); SHA-256 **0e10ce0f7625d454d0574a6ed3e2d0bb27cddf153739c7d665455c878b6f5841**; source [Gammapy data table](https://github.com/gammapy/gammapy-data/blob/main/catalogs/lhaaso_cat.csv). Published paper [Cao et al. 2024, ApJS 271, 25](https://doi.org/10.3847/1538-4365/acfd29).
- **Exact input URLs, retrieval UTC, file hashes, row counts and failure statuses:** [output/receipts.json](output/receipts.json), covering 4FGL v35, ATNF download/database, official Cambridge Green list, eROSITA-DE hard 2024, HGPS 2018 and [Alfaro et al. 2026 4HWC](https://arxiv.org/abs/2602.00263).
- **Full per-source/per-class chance rates and uncovered flags:** [output/source_class_controls.csv](output/source_class_controls.csv); **all candidate counterparts, distances, sizes and pulsar power ratios:** [output/candidate_pairs.csv](output/candidate_pairs.csv).
- **Class census:** [output/class_census.csv](output/class_census.csv); **partial-core unassociated census:** [output/strict_dark.csv](output/strict_dark.csv); **full prose report:** [output/REPORT.md](output/REPORT.md); **13/13 independent assertions:** [output/referee.json](output/referee.json).
- Scripts: [pipeline.py](pipeline.py) gets and freezes catalogues, [hawc4.py](hawc4.py) validates the 2026 HAWC paper table, [controls.py](controls.py) calculates source matches/controls and 25–100 TeV energetic screens, and [referee.py](referee.py) reopens the original bytes and checks source units, pair distances, finite-N probability bounds and claim limits.

## Null and inference restrictions

- **Unit:** one physical 1LHAASO source, never two KM2A/WCDA rows counted independently. Original dark means *no listed TeV counterpart in the preprint's source association field*, not "unobserved at all wavelengths".
- **Matching radius:** clip(source 95%-statistical localization error + 39%-containment scale + 0.2° allowance, **0.35–1.5°**), plus counterpart extension. This is a reproducible screening choice, *not* a true likelihood-ratio model of morphology or association.
- **Chance rate:** (1 + number of matching controls)/(1 + number of eligible controls), from **512** common deterministic 5–355° longitude shifts at **exactly matched Galactic latitude** (seed **20261119**). For all-sky classes, the class-level p counts how often shifted skies equal or exceed the actual number of matches. A 1/513 p is a **finite-simulation floor**, not 3-sigma evidence for an accelerator.
- **Major confound:** fixed latitude does not match **Galactic longitude, Fermi diffuse background, longitude-dependent coverage, source crowding or gamma-ray exposure**. There are multiple classes, possible many-to-one matches and pre-existing model proposals. Fermi's nominal GeV excess **cannot** be interpreted as 19 physical associations or discovery-level significance.
- eROSITA-DE sky coverage is a rough l=180–360° mask; SMGPS and H.E.S.S. use bounding-box masks; 4HWC uses an approximate declination window. Global p-values for these are withheld. HGPS 2018 is not a newly available 2025/26 survey.
- Pulsar budget uses published/derived Ė and (when reported) distance, the fitted KM2A 25–100 TeV power-law integral and isotropic-equivalent luminosity. The efficiency screen checks a *necessary* power bound only; no age/diffusion/spectral/morphology proof exists.
- **Novelty register:** [literature.csv](literature.csv) records prior 2024–26 models/associations (Xie et al. 2025, the 2026 LHAASO Peanut and MSP J0218+4232 hypothesis, Cao's original 4U 1957+11 suggestion, and 2026 4HWC). It is **finite, not a full per-source ADS/SIMBAD literature search**.

## Completed versus pending independent-referee checks

- [x] Claim-first draft PR; labelling gate read in [docs/REMARKABLE_FINDINGS.md](../../docs/REMARKABLE_FINDINGS.md); **no discovery-level label or findings-log entry**.
- [x] Frozen primary transcription SHA-256, 180 rows, 90 unique positions, 43 UHE flags; 25 original no-TeV and 10 UHE no-TeV.
- [x] Source-resolved 512 matched-latitude sky controls, source/pair output tables and finite-sample/null caveats.
- [x] Fermi DR4, Green October 2024, ATNF public dump, eRASS1 hard, legacy HGPS and 2026 4HWC loaded with original URL/hash receipts; pulsar power/distance/size screen.
- [x] Independent-file review: **13/13** assertions passed; failed/pending catalogues preserved as unknown rather than nondetections.
- [ ] Obtain **Cao 2024 original VizieR catalogue** and compare every source/association/TS100 to the 2023 transcription; resolve possible evolution since preprint.
- [ ] Ingest **eRASS1 soft/main v1.2**, release-mask and point/extended-source classification; currently no full X-ray nondetections.
- [ ] Ingest **MeerKAT SMGPS source and SNR candidate tables**, **ASKAP/EMU**, **full SNRcat** and **Gaia DR3 vetted young clusters** with survey-specific completeness/exposure footprints; newer H.E.S.S. products also require review.
- [ ] Assess per-source TeV morphology/spectra, joint GeV/TeV associations, background-corrected longitude-stratified controls, trial factors, independent follow-up and full ADS/SIMBAD literature before any physical ID.

## Reproduce and merge criteria

[GitHub Actions: one push-triggered, path-filtered job](../../.github/workflows/lhaaso-dark-counterparts.yml) downloads public data **on the Actions runner**, computes controls/energetics, performs its independent referee, uploads an evidence artifact and commits inputs, controls and findings in **three bounded commits**. The workflow filters its generated output/ folder out of push triggers, avoiding commit loops. No network request is required or attempted in the local Python sandbox.

**Disposition:** bounded catalogue evidence only; **no newly identified PeVatron**. Pending datasets, exposure maps and the original 2024 input forbid promotion to a discovery result.
