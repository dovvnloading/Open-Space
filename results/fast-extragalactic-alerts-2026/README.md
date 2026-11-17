# Fast extragalactic optical transients in the November 2026 public alerts

**No new extragalactic transient is validated in the primary search.** The strongest unresolved Rubin source is 171198361324486700. Its difference flux faded from 5,900 ± 372 to 4,358 ± 307 nJy in the i band in 67.5 minutes on November 12 UTC. The 0.329 ± 0.103 mag change also appears in fixed-aperture photometry. Legacy Surveys resolves no host at the position and places a pre-existing point source 0.11 arcsec away. No distance, rise time or total duration has been measured. Lasair/Sherlock predicts a variable star on the basis of that same point source. [Lasair record](https://lasair.lsst.ac.uk/object/171198361324486700/)

| Primary search, November 1–16 UTC | Rubin | ZTF |
|---|---:|---:|
| Broker objects first detected in the interval with at least two catalog detections | 22,365 | 41,630 |
| Detection rows inspected | 71,725 | 90,769 |
| Objects with two eligible same-band positive measurements 0.02–1 d apart | 1,074 | 389 |
| Objects passing all photometric cuts | 23 | 9 |

Seventeen Rubin sources have a passing pair without the listed image or Solar System flags. All nine ZTF objects have a stellar counterpart, an image residual near a bright star, or a Solar System identification. Two IDs on successive nights trace asteroid (18150) to within 0.14 arcsec at six observed positions. The ZTF public-program verification covers all 923 eligible pairs. [ZTF results](ztf/ZTF.md) · [Horizons residuals](context/horizons_18150.csv)

The opening commit froze the observing interval, independent same-band positive pairs, S/N ≥5 at both epochs, 0.02–1 day between epochs, absolute change ≥0.30 mag and ≥3 combined errors, and an absolute slope >1 mag/day. Both signs of optical change were retained. In practice the primary cohorts use first broker detection of either sign. An old negative detection can therefore keep a newly positive object out of this table. [Frozen protocol](https://github.com/dovvnloading/Open-Space/blob/f351c39959d602b60fcfe85c477f792ca6e59313/results/fast-extragalactic-alerts-2026/README.md)

## Completed

Archived the photometry, screened all eligible pairs, checked the public ZTF exposures, examined survey stamps, and compared optical positions against localized triggers. An illustrative synthetic transient peaking at AB 18 with 2 mag/day rise and decay passed the two-point screen in 4.27% of Rubin and 0.436% of ZTF injections on the recorded alert epochs. Effective area-time and volumetric rates remain unmeasured. [Recovery](recovery/results.csv) · [Trigger comparisons](coincidences/README.md) · [Published-event controls](controls/README.md)

## Pending at this checkpoint

The final object-by-object audit, supplementary older-object tests, and linked cross-survey evidence are being assembled. No photometric class or discovery priority has been established for the compact Rubin source.
