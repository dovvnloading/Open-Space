# Significant-discovery leads: follow-up results

**23 November 2026.** Base revision: `972cd5169eea2b83d7407df84eb60c167b3c084e`.

This pass advances the selected leads with an additional public JWST epoch,
revised microlensing observing forecasts, fresh checks of affected catalogs,
and concrete follow-up material. The central distinction remains scientific:
the observations and numerical discrepancies can be reproduced, while the
black-hole and transient classifications still require further evidence.

## Most consequential advances

### 1. An additional epoch constrains the GOODS-S brightening

A fresh MAST inventory revealed **JWST program 5324, observation 27, on
4 February 2026**, omitted from the earlier program-limited history of
GOODS-S site 185. A protocol fixed the target, filters and photometry method
before examining the new flux. All four requested filters have usable
mosaic pixels.

In F277W, the 0.15-arcsecond aperture contains **44.46 ± 3.37 nJy** in
February, compared with **70.07 ± 1.81 nJy** in the deposited October epoch.
These errors describe local displaced-aperture scatter, not the full
photometric uncertainty. Seven nearby comparison sources show a median
February/October ratio of 0.933; the nucleus is 0.635. A common scale offset
therefore cannot explain its entire difference.

The useful outcome is a **provisional last-low to first-high interval of
4 February–9 October 2026**, about 247.5 observer-frame days. It does not
identify the onset, peak or physical cause. Native exposures, pixel-quality
checks, relative astrometry and source-specific PSF effects remain the next
validation step. [Measurements, controls and receipts](transients.md#additional-goods-s-185-epoch)

### 2. The black-hole candidate needs spring astrometric planning

The 1,900-point OGLE fit and all 24 rows of its frozen photometric forecast
reproduce. Propagating **all 4,800 stored posterior samples per geometry**
through the annual-parallax trajectory replaces the old rectilinear July
estimate for maximum source-centroid displacement:

| Geometry | Median maximum | Central 95% interval |
| --- | --- | --- |
| Positive impact parameter | 6 May 2027 | 1 March–21 July |
| Negative impact parameter | 11 April 2027 | 21 February–18 June |

These are conditional predictions from saved photometric samples, with
correlation and model limitations recorded in the packet. They describe an
unblended source; the measured unresolved photocenter also depends on
neighbors, their colors and motion. **No lens mass has been measured.**

A geocentric solar-elongation check gives approximate JWST windows of
15 March–5 May and 6 August–27 September 2027. July lies outside the field of
regard. Spring observations should be assessed promptly, and late-summer
observations remain valuable. The astrometric tail can persist after the
light curve returns near baseline, so later images need a lensing model too.
[Full calculation and proposed observing package](ogle.md)

A subsequent run of the official GTVT 1.2 package with a live JWST ephemeris
supported these seasons but exposed a one-day mismatch between its date
labels and position vectors. The packet preserves the unmodified output and
a date-aligned audit. Exact visit boundaries still require APT checks; this
does not change the spring and late-summer observing priorities.

### 3. The catalog corrections remain reproducible

| Lead | Work completed now | Narrow conclusion |
| --- | --- | --- |
| NICER magnetar catalog | Six compact detector/timing/forensics replays | Evidence supports a +69.184-second correction to the published absolute labels; the hardware clock and physical origin of the flashes are separate questions. |
| MeerKAT SMGPS image index | Fresh 693-row Harvard VizieR retrieval, coordinate reconstruction and WCSLIB comparison | The affected index values persist at that mirror. The native image WCS is a separate product. |
| NGC 1261b / Reticulum II | Fresh primary-source crossmatch and coverage checks, plus current Galstreams retrieval | All 12 original branch-b velocity stars match Reticulum II lists; the derived track still has 690 unphysical extrapolated samples. This does not disprove the sky/proper-motion stream. |

The [correction packet](corrections.md) contains precise report text,
attachments and falsifiers. A new [58-target spectroscopy manifest](corrections-ngc1261b-targets.csv)
retains the remaining atlas candidates, exact Gaia IDs and the one known
later Awad spectrum. It is a follow-up selection, not a membership catalog.

### 4. The infrared leads survive, with clearer interpretation limits

All **113 NEOWISE rows** reproduce, including the same **seven strict
candidate flares**. UGC 2653 and CGCG 539-025 remain the most useful nearby
follow-up targets in this packet. The UGC radio contrast is 3.22 times its
statistical error; an illustrative independent 10% epoch-calibration term
reduces that ratio to 2.40. Its SPHEREx color is an outlier relative to 65
controls, rather than a calibrated nine-sigma discovery probability.

The result page, manuscript and message now acknowledge that coronal lines
also occur in interacting supernovae, with SN 2005ip as a primary-literature
counterexample. Spectral profiles, ratios, evolution, localization and radio
behavior need to be interpreted together. [Ranked targets, radio scenarios
and classification tests](transients.md)

## Reviewable outputs

| Material | Location |
| --- | --- |
| OGLE calculations, observing requirements and unchanged photometric test | [ogle.md](ogle.md), [code](ogle_checks.py), [numerical receipt](ogle-evidence/checks.json) |
| JWST intermediate-epoch measurements and infrared follow-up | [transients.md](transients.md), [frozen target protocol](transients-goods5324-protocol.json), [photometry](transients-goods5324-photometry.csv) |
| NICER, SMGPS and NGC correction reports | [corrections.md](corrections.md), [replay receipt](corrections-replay.json), [live-source receipt](corrections-live-checks.json) |
| Revised research-note sources and publication status | [Publication review](publication.md), [existing paper directory](../../papers/rnaas/README.md) |
| Updated observing and notification drafts | [Astrometry request](../../papers/rnaas/messages/02-astrometry-followup-email.md), [infrared request](../../papers/rnaas/messages/03-flare-followup-email.md), [ATel draft](../../papers/rnaas/messages/04-atel-ogle-2026-blg-0058.md) |

## Verification and outstanding work

The packets distinguish replays of the same observations from newly examined
archive data. OGLE's live photometry download returned HTTP 503; no new OGLE
data were substituted into its frozen predictions. The JWST February result
uses newly examined public pixels. Other fresh checks compare current catalog
products with their archived versions. Literature searches are bounded and
do not establish discovery priority.

An additional research agent reviewed the new calculations. That review
caught and resolved extra thinning of the already thinned OGLE sample file,
and a JWST extraction failure path that could leave an obsolete output CSV.
The final calculations, hashes and failure behavior were checked again.
Agreement between agents is an implementation review, not peer review or an
independent astronomical observation.

The next scientific actions are native-exposure validation of the February
JWST point, an instrument-specific astrometric feasibility assessment for
OGLE, calibrated radio/spectroscopic observations of the nearby flares, and
author/maintainer review of the catalog evidence. The packet contains drafts
for those discussions. Telescope allocations, external correspondence,
transient registrations, journal submissions and external catalog changes
have not occurred in this pass.
