# Referee: unconfirmed 4XMM J134144.5−005832 (ObsID 0841482101)

**Do not call a new QPE.** Frozen DR14 DETID **108414821010004**, SRCID **203057509010014**, position **RA 205.435623°, Dec −0.975584°**, FITS-coordinate separation **0.886″**. This is an *off-axis source in a quasar pointing*, not the observation target LBQS 1338−0038. No association to the target QSO should be inferred.

## NASA PPS time-series evidence

- EPIC-pn merged NET source file **P0841482101PNX000SRCTSR8004.FTZ** from `https://heasarc.gsfc.nasa.gov/FTP/xmm/data/rev0/0841482101/PPS/`, FITS SHA-256 **0dbbf228758e4fd9c367f7ad2174118a962aad94f1fcfc0de4bc9350b1b2a47d**; see original per-GET retries, compressed hashes and source position in [primary batch receipt](heasarc_batch_cohort_495_670.json). The catalog pn pile-up diagnostic is **0.06640594**, but **SAS epatplot not performed**.
- Fixed finder 300s 0.2–2keV NET peak results:
  1. **MJD TT 58880.880208**, 0.14201 cps, local background 0.01917 cps, excess 4.06σ and multiscale integrated 10.61σ.
  2. **MJD TT 58881.112847**, 0.16667 cps, local background 0.01311 cps, excess 4.05σ and integrated 10.19σ.
  3. **MJD TT 58881.498264**, 0.17452 cps, local background 0.00523 cps, excess 4.35σ and integrated 12.20σ.
- Nominal successive separations **5.583h and 9.25h**. The source-masked particle-field background is **100.31 counts/s** at the first peak, automatically flagged as a field particle flare. Its next two values at source peaks are **3.97 and 3.53 cps** and are *not* individually flagged. Thus the three-peak interpretation fails the independent field-veto, while the two apparently clear peaks require true soft-band event-background confirmation, GTI alignment, exposure checks and spectral evolution before interpreting 9.25h as a QPE interval.
- Source-specific scaled off-source aperture BACK1–3 curve had **zero recurrent two-peak flags**, but is not an independently extracted matched-energy event-region light curve or empirical false-positive measurement. Two count-domain surrogate synthetic pairs recovered; this is **not event-list completeness**.

## Host, confounders and novelty

- Frozen SDSS DR18 crossmatch: optical **PhotoPrimary galaxy**, objID **1237671990268723298**, SDSS `SpecObjAll` sciencePrimary spectroscopy **GALAXY**, reliable **z = 0.1476348**, `zWarning=0`. This supports a galaxy optical counterpart but not a complete nuclear-association or AGN-subtype diagnosis.
- Gaia DR3 source **3661858914360014848**: parallax **−0.4154 mas**, significance **−0.635σ**, proper motions **−0.040 ± 0.763** mas/yr (RA), **−0.230 ± 0.504** mas/yr (Dec), RUWE 1.18. No significant astrometric star veto; not a full stellar/variable-star census.
- The *pointing target* LBQS 1338−0038 is a studied quasar, but its position does not coincide with this off-axis X-ray detection; **do not falsely classify this source as the QSO**.
- Current partial audit does **not** establish negative SIMBAD/4LAC matches, eROSITA eRASS1/Swift/Chandra multi-epoch variability or TNS/ATel novelty. See [automated follow-up receipts](heasarc_trigger_catalog_vetting.json) when present.
- **QPE vs TDE vs AGN:** all remain open after partial particle veto. At least two candidate soft peaks survive an initial particle diagnostic, but their thermal spectral behavior, event pattern/pile-up, false-alarm probability and genuine off-source soft-band net extraction are **unverified**.

**Disposition: unconfirmed variability / follow-up priority; zero validated new QPEs. No remarkable label or merge until instrumental and novelty gates close.**
