# eRASS1 DR1: search for X-ray-dim isolated neutron stars

**Claim (provisional; NOT a discovery): no previously unknown XDINS has yet been validated by this analysis.** This is a search in progress, **not** a calibrated null result. Candidate counts and exclusions must come from the reproducible outputs, not from the existence of an empty list.

## Frozen scope and evidence bar
- Dataset: public eROSITA-DE **eRASS1 DR1, western Galactic hemisphere**, *main* 0.2–2.3 keV and *hard* 2.3–5.0 keV source catalogues; use the exact catalogue release, URL, bytes, SHA-256 and retrieval UTC logged by the runner (the provider lists main **v1.2** as of 2026-01-12). No later DR2 sources.
- Reject low-quality, extended, variable and hard-band emitters; screen for soft, point-like, sufficiently bright and significant detections. Assess sensitivity with injection/positive controls before any completeness claim.
- Require meaningful multi-survey non-detections in Gaia DR3, Legacy Surveys DR10, Pan-STARRS1, DES, unWISE/CatWISE and 2MASS, taking into account **survey sky coverage and depth**. No match is *not* a limiting magnitude; an XDINS-like f_X/f_opt ≳ 10^4 requires an actual optical limit and band definition.
- Screen coronal stars, AGN, radio/IR sources, photometric-redshift galaxies, cluster extent, confusion and spurious detections.
- Compare candidates with known XDINS/INS in the DR1 footprint, Kurpas et al. 2024 (33 published eRASS candidates) and subsequent follow-up work. Inspect XMM/Chandra/Swift archival matches, thermal spectra (kT ~40–120 eV) and timing evidence where available.

## Reproduce and audit
Workflow: `.github/workflows/erass1-xdins.yml`. Analysis: `results/erass1-xdins-search/run.py`. Dependencies: `requirements.txt`. The single push-triggered Actions job must generate `receipts.json`, `screening.json`, `candidates.csv`, `controls.json`, and `checks.json`. Hashes refer to the **raw inputs**; the main catalogue is too large to commit. Generated files are not allowed to retrigger this workflow.

**Completed:** research scope and discovery labelling policy read; catalog and literature endpoints identified. **Pending:** catalogue download and checksum; measured screen, cross-survey coverage/matches, positive controls, archival/literature follow-up, independent referee review. No discovery-level claim or entry in `docs/REMARKABLE_FINDINGS.md` without validation.

Sources: [eROSITA DR1 official catalogues](https://erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/Catalogues_dr1/), [eRASS1 main data model](https://erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/Catalogues_dr1/MerloniA_DR1/eRASS1_Main.html), [Kurpas et al. 2024](https://arxiv.org/abs/2405.12846).
