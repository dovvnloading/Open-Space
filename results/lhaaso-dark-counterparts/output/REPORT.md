# 1LHAASO counterpart census: checked run output

**Result:** positional/energetic candidate screen only; no new accelerator identified.

Frozen Cao 1LHAASO transcription: 90 physical sources, 43 UHE, 25 lacking an original TeV association; 10 of the latter UHE.
A negative association is not a survey non-detection unless that class has a validated footprint.
**Primary-data limitation:** Gammapy transcription of the 2023 Cao preprint is frozen and hashed; the 2024 CDS VizieR bytes could not be verified during the Nov 2026 CDS outage.

## Sky-shifted class controls

| Class | Status | Covered sources | Observed matches | Mean control matches | Tail p |
|---|---|---:|---:|---:|---:|
| GeV_4FGL_DR4 | screened | 25.0 | 19.0 | 9.46289062 | 0.00194931774 |
| pulsar_ATNF | screened | 25.0 | 16.0 | 10.1503906 | 0.157894737 |
| SNR_Green2024 | screened | 25.0 | 9.0 | 5.38085938 | 0.20662768 |
| SMGPS_new_SNR_candidate | catalogue_pending |  |  |  |  |
| SMGPS_known_SNRcat_subset | catalogue_pending |  |  |  |  |
| eRASS1_hard_Xray | approximate_footprint_only | 1.0 | 0.0 |  |  |
| Gaia_DR3_young_OC | catalogue_pending |  |  |  |  |
| TeV_HESS_HGPS_legacy | approximate_footprint_only | 4.0 | 0.0 |  |  |
| TeV_HAWC4_2026 | approximate_footprint_only | 25.0 | 13.0 |  |  |

Individual source/class probabilities: source_class_controls.csv. Pair offsets, conditional chance rates and pulsar efficiency/size tests: candidate_pairs.csv.
Chance probability uses (1 + matched sky-shift controls) / (1 + eligible controls). 512 shared longitude offsets; Galactic latitude held fixed, shifts >= 5 degrees. Source aperture = clip(ePos95stat + r39 + 0.2 deg, 0.35 deg, 1.5 deg); catalogued counterpart radius added.

## Critical limitations

- This chance test matches Galactic latitude but not local longitude, source confusion or exposure in each survey. It is a screening null, not a source-identity posterior.
- eROSITA/SMGPS/HESS survey masks are approximate; their global p-values are withheld. eRASS1 hard is not the 930203-source soft main catalogue. Only the German half-sky is public.
- SMGPS SNR candidates are not the 510599-source compact radio catalogue. The project's prior SMGPS image-index astrometry correction is relevant when using maps.
- 4HWC (Alfaro et al. 2026 arXiv:2602.00263) is the newest HAWC catalogue; its TeV detections and association proposals pre-empt novelty. Declination acceptance is approximate; no global HAWC4 p-value is quoted.
- Fermi DR4 v35 is a corrected 2024 release of a catalogue first published 2023; HGPS is 2018 legacy. Never attribute new detections to a survey they predate.
- ATNF spin-down and distance allow a necessary power/size test, not causality; KM2A 25-100 TeV power-law flux is not a measured bolometric budget.
- Multiple correlated source/class hypotheses and crowded Galactic longitude invalidate a naive discovery significance from a small nominal p.

## Pre-existing 2024-2026 literature claims (finite, not exhaustive)

| 1LHAASO source | Reference | Earlier suggestion |
|---|---|---|
| 1LHAASO J0007+5659u | [paper](https://doi.org/10.1051/0004-6361/202452796) | Discussed among six isolated UHE sources; no secure counterpart |
| 1LHAASO J0206+4302u | [paper](https://doi.org/10.1051/0004-6361/202452796) | Traveling PWN/diffusion modeling for the three-source high-latitude complex |
| 1LHAASO J0212+4254u | [paper](https://doi.org/10.1051/0004-6361/202452796) | Traveling PWN/diffusion modeling for the dumbbell-like region |
| 1LHAASO J0216+4237u | [paper](https://doi.org/10.1051/0004-6361/202452796) | PSR J0218+4232 candidate already discussed; not established |
| 1LHAASO J0206+4302u | [paper](https://doi.org/10.1016/j.scib.2026.09.046) | Peanut complex includes high-latitude 1LHAASO sources; millisecond pulsar J0218+4232 suggested |
| 1LHAASO J0212+4254u | [paper](https://doi.org/10.1016/j.scib.2026.09.046) | Peanut complex; authors report positional chance 3.6e-4 for J0218+4232; source morphology supersedes point matches |
| 1LHAASO J0216+4237u | [paper](https://doi.org/10.1016/j.scib.2026.09.046) | MSP J0218+4232 potential PeVatron already proposed by LHAASO; not newly discovered here |
| 1LHAASO J1740+0948u | [paper](https://doi.org/10.1051/0004-6361/202452796) | PSR J1740+1000 already considered; physical identification uncertain |
| 1LHAASO J1959+1129u | [paper](https://doi.org/10.3847/1538-4365/acfd29) | LMXB black-hole candidate 4U 1957+11 suggested at about 0.23 deg in original paper; NOT a new microquasar identification |
| All 1LHAASO candidates | [paper](https://arxiv.org/abs/2602.00263) | 4HWC 85 TeV sources and 22 unpaired HAWC sources with 1LHAASO; prior HAWC TeV detections not physical accelerator IDs |

## Independent referee

13/13 assertions passed. Discovery-level label: **not permitted**.
Pending catalogue/coverage checks: askap_emu, erass1_dr1_main_soft, gaia_dr3_hunt_ref, meerkat_smgps_compact_510599, smgps_new_snr_candidates, snrcat_full, snrcat_smgps_overlap.
Still requires independent TeV morphology/spectrum, all-region exposure masks, full per-source literature checks, and targeted multiwavelength observations.

Frozen source bytes, normalized catalogue extracts, byte hashes and original URLs are preserved in this folder; see receipts.json.
