# Independent referee: VAST DR1 60-lead follow-up

**Disposition: accept as a reproducible, bounded native-Stokes-I image and catalogue follow-up; reject any claim of a discovered long-period transient or of a calibrated negative LPT search.** Reviewed on 19 November 2026 UTC after the last GitHub Actions evidence run.

## Audit of evidence and controls

The independent machine referee in [referee.py](referee.py), supported by an additional [source-class audit](reclassify.py), verifies exactly **60** distinct source IDs, **693** original timestamp-matched image epochs, **14** original on/off patterns and **15** ATNF v18 known-pulsar matches (three of the 14). The 60 DAP native FITS source cubes match their pinned public URL index and carry per-file SHA-256 receipts. All 693 WCS-centered cutouts were assigned unambiguously to measurement epochs. The astronomical judgement below is independent of the pipeline's self-reported success, but **does not substitute for a human visual referee of each FITS extension or full-field calibration images**.

The diagnostic ledger contains **122** epochs with possible bright-neighbour/sidelobe flags from **17** sources (nine epochs also flag deep negative noise). These are *possible* artefacts, **not proven instrumental contamination**, and do not demonstrate real variability. Every cutout was centered away from its *cutout edge*, which says nothing about edge proximity in the parent ASKAP pointing or sidelobes from sources outside the postage stamp. Native brightness alone cannot validate all 14 on/off patterns.

The Stokes-I RACS low, mid, high position catalogues yielded 46, 47 and 42 available source matches, with **45** two-or-more-band slopes. These are noncontemporaneous slopes and can be biased by different epochs, selection and blending, so they do not establish a spectral energy distribution or an AGN/scintillator identity. In particular a steep positive slope is not by itself coherent radio emission.

Gaia DR3 has 12 position matches within 3 arcseconds, but only eight **non-pulsar** associations are within 2 arcseconds. A second-stage referee correctly rejects treating all eight as flare stars: just **two** have strongly stellar supporting evidence (Gaia G<8 or >5-sigma parallax, source IDs 13629007 and 13572974); **six** remain positional associations of undetermined astrophysical nature. Even the two plausible stars have **no proven radio flare association**.

The archive-level direct counterpart queries returned **33 unWISE** position matches; the Pan-STARRS DR1 mean table returned **30 positional rows**, but only **13** satisfy the minimum requirement of at least two detections and a valid PSF-band magnitude. No credible full Legacy DR10 Tractor crossmatch was obtained because the TAP response was not reliably parseable (60 unavailable, not 60 nonmatches). eRASS1 has **eight** source matches within 12 arcseconds; Swift 2SXPS yielded **zero catalogue entries** in the queried 12-arcsecond circles, which says nothing about source exposure or X-ray flux limits. Optical/IR source matches are not guaranteed physical identifications. Proven AGN/scintillator classes: **zero**; confirmed flare-star classes: **zero**; secure artefact classes: **zero**.

### Decisive missing tests

The published VAST DR1 source cubes contain Stokes I, **not Stokes V**. The bounded public RACS-low2 DAP file inventory returned no usable small Stokes-V source file; this is a failure of accessible measurement, **not a polarisation non-detection**. Calibrated RACS-low/mid/high and VAST Stokes-V cutouts and position-dependent I-to-V leakage assessment are outstanding. The sparse catalogue epochs cannot test minute-to-hour periodicity, and no source met the prerequisite of *independently verified compactness, strong circular polarisation and a secure absence of counterparts*. Therefore an accelerated snapshot/visibility period scan and object-specific ADS/ATel claim searches were **not justified by this evidence** and were not undertaken.

## Verdict

**Final classes: 15 known pulsars; 2 plausible stellar associations (flares unproven); 6 ambiguous Gaia associations; 37 otherwise unresolved.** Among the 14 on/off patterns: **3 known pulsars, 1 plausible stellar association, 10 unresolved**. No validated new LPT and no defensible calibrated null rate. Treat the remaining sources as **follow-up leads**, not discoveries; preserve all 60 classification rows and 693 epoch rows. The standalone [referee.json](output/referee.json) validates accounting/anti-overclaim consistency, **not independent human visual or polarimetric verification**.

**Recommendation:** merge the scoped observational audit, report all missing survey products, and make no exceptional-discovery designation or log entry.
