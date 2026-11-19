# Independent referee report — Euclid DR1-F cold T/Y dwarf screen

**Decision: ACCEPT as a reproducible, explicitly incomplete catalogue-screening result; REJECT any discovery claim, temperature record, clean-survey null, sample-completeness or effective-volume interpretation.** The source-hash replay, comparison population and attempted external receipts support the limited claim in [README.md](README.md). This assessment is deliberately stricter than ranking a VIS-dropout as a Y dwarf.

## Checks independently reproduced

The network-free, independent [referee.py](referee.py) read and SHA-256-verified all **32** sequential HATS batch packets and validated the frozen `selection.py` fingerprint. The 32 manifests cover **768 unique partitions**, the advertised **888,152,165 sources**, and **647,197 independently re-tested strict-selected objects** without a duplicate MER ID. It checked the retained selection limits and shortlist IDs. This reproduces the [machine referee report](output/referee.json) status **`INDEPENDENT_RECEIPT_REPLAY_PASS_NO_DISCOVERY_CLAIM`**, not a PSF, spectroscopic or native-image review.

The control replay reports the following conditional numerators and denominators. They must **not** be divided by all frozen literature positions as though the entire footprint were imaged at uniform quality:

| Class | Frozen controls | MER matches ≤3″ | Pass independent strict mask |
| --- | ---: | ---: | ---: |
| Known T8/T9/Y | 256 | **20** | **2 (10.0%)** |
| Known T6/T7 | 266 | 10 | 0 |
| High-z quasars (contamination controls) | 120 | **24** | **3 (12.5%)** |

Both recovered late-T/Y controls have credible CatWISE mirror counterparts. The known T8.5 has a **6.74σ** CatWISE motion component and W1−W2 **3.046 Vega**, and known Y0 has a **1.12σ** component and W1−W2 **3.093 Vega**. Their propagated separations are **0.526″** and **0.996″**. These tests show mirror data access and at least one properly motion-supported positive; they **do not** validate a motion cut on all Y types. The 3/24 quasar passage and 2/20 cold-control passage are powerful falsifiers of purity and completeness.

The independent screen's **647,197** colour-mask rows, including **368,529** with pathologically extreme J−H <−1.4, are *not* astrophysical classification results. The heuristic 250-entry file is a **screening register**, not a 250-object new-dwarf catalog.

## External-evidence disposition

Of **40** highest-priority photometric rows queried, CatWISE2020 at NOIRLab returned **35** cones with rows and **5** with none, with **zero qualifying motion-propagated matches** for the ranked leads. The original IRSA proxy returned 40 HTTP 502 errors. An unWISE alternative was attempted using the original IRSA service and a mirrored full/reduced-column ADQL query, but **all 40 remained service errors**. Consequently a missing CatWISE/unWISE association is neither a WISE forced-photometry upper limit nor an exclusion of an otherwise faint candidate. Legacy DR10 and Pan-STARRS catalogue returns do not establish verified non-detections in ground-image pixels, and 40 empty SIMBAD mirror responses are not a complete literature search. See [query receipts](output/external_receipts.json), [earlier error receipts](output/first_pass_irsa_proxy_failures.json), [second-stage failures](output/second_pass_unwise_schema_failures.json), and [CatWISE positive-control queries](output/wise_control_probes.json).

**No claimed new source has all of**: an independently confirmed point-source NISP PSF, clean VIS and NISP science pixels from multiple exposures, a significant held-out proper-motion track, convincing methane/W1−W2 and J-band photometry with contaminant rejection, a confirmed infrared spectral subtype and a complete publication check. High-z quasars, stationary red galaxies, blended sources, detector artefacts and poor flux estimates remain plausible contaminants.

**Referee ruling for physical quantities:** do not assign Teff or a distance to any unclassified rank from Euclid two-colour measurements. Published-parallax `M_W2` medians by T8/T9/Y0/Y1 are conditional anchors only; they neither uniquely classify a candidate nor establish any sub-500 K object or a cooler-than-WISE-0855 competitor. A very cold Y dwarf may also fall outside this simple Y−J/J−H box. The actual effective survey area, magnitude-dependent limits, selection function and Y spectral-type recovery are not calibrated.

## Merge eligibility

**PASS the reproducibility/documentation bar** as a **bounded, inconclusive full-release *catalogue* measurement**. **FAIL the discovery-level remarkable-findings bar.** The reproducible null is *zero motion-and-methane-supported new objects among the top 40 externally queried photometric priorities*, **not** zero new Y dwarfs over 1,900 deg². Hold all 250 entries as unconfirmed photometric leads and record **zero validated novel T/Y objects**.

Required follow-on observations for a discovery claim: per-source NISP PSF and instrument/pixel-masked image examination, a distinct-epoch astrometric detection, actually measured WISE/Spitzer/JWST colors or spectra with bandpass-appropriate empirical comparisons, full 2026 literature completeness, and injection/recovery controls. Those checks are absent, so **no `[REMARKABLE]` tagging or `docs/REMARKABLE_FINDINGS.md` entry** is warranted.

[Source and runbook](RUNBOOK.md) · [machine checks](output/referee.json) · [control recovery](output/controls_recovery.json) · [full sample summary](output/summary.json) · [40 follow-up receipts](output/external_receipts.json) · [catalogue index and selection](selection.py).
