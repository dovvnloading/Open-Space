# HEASARC pn recurrent-peak veto: XMMU J134736.6+173403 (0144570101)

**Referee classification: KNOWN SEYFERT-2 AGN / PARTICLE BACKGROUND FLARE; NOT A NEW QPE.** This is an **unvetted finder trigger**, not a validated QPE candidate. Never label [REMARKABLE].

## Observational evidence

- Frozen host detection **101445701010002**, SRCID **201445701010002**, 4XMM position RA **206.901640°**, Dec **+17.567958°**; SDSS DR18 galaxy centroid offset **0.0504″**. The measured pn FITS source coordinates are within **1.029″** of the catalogued X-ray position.
- NASA HEASARC mirror PPS product **P0144570101PNX000SRCTSR8002.FTZ**; 659,919 compressed bytes, compressed SHA-256 **0678c259c5dc040d624e43450f695b14cf71767d6e376f2a608ceaef43a7c6af**, uncompressed FITS SHA-256 **43e8f0cddd6aab896c09633f18250103d9bc2c36427ed20fe1fbcce5d188585c**; HDUCLAS2 **NET**, explicitly coordinate verified, RATE1+RATE2+RATE3 soft 0.2–2 keV, **197 good 300-s bins** (14.745h of exposed data).
- Frozen finder identifies two nominal 0.2–2 keV excursions, at **MJD TT 52814.217014** and **52814.293403** (separation **1.833 hours**), with peak NET rates **0.34088 / 0.38031 counts/s**, versus estimated quiet 0.039408 counts/s. This is not a published QPE period or recurrent thermal-eruption determination.
- NASA PPS source-masked **P0144570101PNS003FBKTSR0000.FTZ** (compressed SHA-256 **1e7000188bd0087019d4d41f7c26962b76570f738586fe50af18d2f6057de285**; FITS SHA-256 **fb843955d66b2a0878d45cf17f3778a3bfde28195f17c409c4a7c3102ba035d9**) shows particle-background rates **3.54** and **29.76 counts/s** at the source peaks, respectively. The registered 5-MAD/3× particle-flare threshold was **11.53 counts/s**, so the **second apparent peak overlaps an off-source field particle-background flare**. The first does not. **This fails the independent source-masked particle-field veto for a two-peak train.**
- The scaled source-specific BACK1–3 off-source aperture series produces **0** repeated-peak flags. This is supporting, not definitive: the proxy uses source-curve uncertainty and does not replace independently extracted soft-band background regions, particle GTIs or event grades. SAS epatplot pile-up and soft-temperature evolution were not completed. Two of two count-domain synthetic pair injections recovered, **not event-list completeness**.

## Independent host and literature checks

- The frozen SDSS DR18 primary-galaxy match is nuclear; Gaia DR3 positional counterpart **1244596486324206336** has parallax significance **1.07σ** (not a 5σ astrometric star flag), though this alone cannot provide comprehensive stellar clearance.
- **Carpano et al. (2008)**, [arXiv:0801.4246](https://arxiv.org/abs/0801.4246), reported this 2003 XMM-Newton source and a **6.5× drop within ~1h**, followed by a long low state. It was already a studied anomalous object, not a new flare discovery.
- **Carpano & Jin (2018)**, [MNRAS 477, 3178](https://academic.oup.com/mnras/article/477/3/3178/4956771) ([arXiv:1803.10976](https://arxiv.org/abs/1803.10976)), established with **2008 Chandra** astrometry that the X-rays originate in **SDSS J134736.39+173404.6**, a **Seyfert 2 AGN**. In **29 Swift observations** they reported **23.82±0.07h and 71.44±0.57h quasi-periodic variability**: neither is the 1.833h triggered separation. The historical radio-loud AGN is **not a newly found QPE**.
- Broader SIMBAD/4LAC and eRASS1, TNS/ATel database crossmatches were **not comprehensively completed**. The existing literature and direct particle flare are sufficient to reject novelty, not to claim external catalog-wide negative matches. TDE vs QPE or new thermal-eruption classifications are unsupported.

**Decision:** exclude from validated new QPE candidates; do not infer an archive null or empirical false-positive rate. Full raw download attempts, SHA-256 and FITS metadata are in [batch receipt](heasarc_batch_cohort_000_080.json). This is **one preliminary trigger, zero validated novel QPEs from this source**.
