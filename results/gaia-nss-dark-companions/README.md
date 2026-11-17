# Gaia DR3 dark-companion screen: three orbital leads, no confirmed compact object

**Standing:** orbital measurements and model-dependent companion candidates only. None passes the inclusion bar for a discovery-level compact object; ordinary luminous companions, stellar variability or suspect orbit/primary-mass models remain viable. Read the [independent referee audit](REVIEW.md), [high-mass coverage audit](high_mass_audit.md) and [individual orbital cases](orbit/CANDIDATES.md).

**HD 264291 (Gaia DR3 3378588057203660160)** has a joint Gaia DR3 NSS–LAMOST DR11 period of **1,007.3 (1,001.8–1,013.1) days**. The 100 quality-selected exposures on **18 nights** give a conditional dark-companion mass of **1.77 [1.71, 1.84] solar masses** (16th–84th percentiles) and a nightly RV residual RMS of **0.69 km/s**. The inference assumes that the spectral lines trace the luminous primary, uses an IsocLum primary-mass prior, and equates its orbit with the photocentre. A hierarchy containing two ordinary, roughly 0.885-solar-mass dwarfs remains consistent with the broad-band SED and conditional spectral tests. The outer orbit is supported; a neutron star, a new identification or exclusion of companion light is **not** established. [Joint fit](orbit/lamost_candidate_fits.json) · [Spectral and photometric alternatives](vetting/README.md)

**UCAC4 720-033289 (Gaia DR3 275582835635329024)** has a validated **470.1-day** single-lined orbit and RV mass function **0.640 [0.606, 0.675] solar masses** from seven LAMOST nights. Gaia FLAME's flagged giant primary mass implies a conditional edge-on companion minimum of **3.60 [3.51, 3.69] solar masses**; the independent LAMOST LASP-based giant model gives **2.81 [2.65, 2.93]** at the fixed central mass function. The 3-solar-mass crossing is therefore **primary-model dependent**, and there is no DR3 inclination measurement. GALEX NUV flux and optical Hα wings are substantially more compatible with a luminous companion than the tested zero-light models. It is **not** a black-hole identification. [Orbit](orbit/CANDIDATES.md) · [Stellar models](extension/followup/MASS.md) · [UV and spectra](extension/followup/NUV.md)

**Gaia DR3 2951241529971223808** shows four quality-passing DESI DR1 velocities over **45.87 days** that follow a predicted **276-day** SB1 trend. A nominal FLAME primary yields a conditional edge-on minimum of **3.47 [3.24, 3.72] solar masses**; a 30% lower primary gives **2.87 [2.67, 3.08]**. Gaia photometry also varies on **138.185 days**, almost half the nominal orbital period. Pulsational variability, ellipsoidal effects and ordinary stellar light are not excluded. [Fit and diagnostics](orbit/CANDIDATES.md) · [Photometric receipts](extension/desi_lpv/photometry.json)

## Screen and completed survey checks

| Screen and follow-up | Deposited result |
|---|---:|
| Gaia NSS sources with valid mass function | 345,265 |
| Gaia sources with applicable IsocLum mass | 188,638 |
| Initial screen: conditional companion mass >1.4 / >3 M☉ | 983 / 112 distinct sources |
| FLAME-pass sources above 1.4 / 3 M☉ | 6,797 / 291 distinct sources |
| FLAME additions not already in the initial >1.4 / >3 sets | 6,750 / 271 distinct sources |
| Deduplicated two-screen union above 1.4 / 3 M☉ | **7,733 / 383** distinct sources |
| Union with any clean independent survey RV / same-instrument repeat | 568 / 257 |
| >3 M☉ union with any clean RV / same-instrument repeat | 24 / 12 |
| Source–solution–instrument repeat-RV shape checks | 138, with no skipped comparisons |
| Union members found in at least one ingested prior publication table | 1,627 |

The completed **FLAME-extension survey validation** comprises LAMOST DR11 clean-exposure extraction and covariance-propagated orbit checks, SDSS/APOGEE/BOSS and DESI DR1 epoch checks, deduplicated independent-night coverage, four positive-control cadence sensitivities, published-ID crossmatches, and follow-up for the giant, variable star and two *known luminous* Kepler eclipsing-binary negative controls. Asteroseismology reclassifies **21 of 24** usable screened giant matches below 1.4 M☉ when their primary masses are substituted. The 138 checks are not 138 independent targets; the sample is sparse and phase-biased. [Union coverage](coverage.json) · [Survey shape checks](orbit/survey_union_shapes.csv) · [Cadence depths](orbit/cadence_union.json) · [Seismic test](vetting/seismic_depth.json)

## Pending checks and limitations

Full bulk-source re-acquisition from ESA Gaia (167 MB NSS, IsocLum binary masses, and the FLAME TAP response) is separate from the deposited offline row-by-row replay; input URLs, publisher MD5 where applicable, exact SHA-256 and query receipts are recorded, and an optional CI step attempts the hash-gated upstream rescreen. Of the **383** above-3 M☉ screens, **359 have no quality-passing external survey epoch** and another **12** have no same-instrument repeat. Independent orbit coverage cannot establish darkness or screen completeness. The complete electronic catalogue of Li et al.'s 2026 **1,835-target** Gaia+RV study was not located; exact-ID negative checks and SIMBAD names do not prove novelty. Homogeneous high-resolution spectroscopy, a reliable primary-mass model, direct companion-light constraints and full phase coverage remain pending for the proposed compact-companion interpretations. [Referee checks, hashes, prior-art statuses and specific follow-up](REVIEW.md)

## Reproduce and inspect

Run the **offline, non-mutating** independent referee audit with Python 3:

```sh
python results/gaia-nss-dark-companions/referee_replay.py
OPENBLAS_NUM_THREADS=1 python results/gaia-nss-dark-companions/orbit/fit_lamost.py --out /tmp/hd264291-refit
python results/gaia-nss-dark-companions/referee_replay.py --check-fit /tmp/hd264291-refit/lamost_candidate_fits.json
```

The two original screens can also be rerun after re-fetching their large hashed Gaia inputs; the [Gaia screen](gaia/README.md) and [FLAME program](extension/flame.py) preserve parameter choices. Check [raw LAMOST receipts](lamost/raw/receipts.json), [extension response receipts](extension/lamost/raw/receipts.json), [published reference-table receipts](novelty/inputs.json), [SIMBAD query receipts](orbit/identification_inputs.json), and [positive-control source receipts](orbit/sources.json). No entry has been added to the remarkable-findings log.
