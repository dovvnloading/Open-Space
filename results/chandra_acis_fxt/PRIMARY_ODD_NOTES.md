# Human review of odd-numbered primary shards

Reviewed worker ledgers from 16 odd-numbered primary shards. The 500 ObsIDs represented contain 5978 raw real triggers and 36 real instrument passes with an unresolved catalog classification. The scope excludes the separately reported moderate-count cohort. [Actions run](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443)

Counts in this index describe the available worker records, including partially processed observations. The combined survey report defines which observations and distinct triggers enter inference.

| Shard | Assigned | Ledgers | Complete follow-up | Raw real | Unresolved instrument passes | Incomplete catalog checks, all streams | Evidence |
|---|---:|---:|---:|---:|---:|---:|---|
| 3 | 18 | 18 | 17 | 1234 | 4 | 70 | [Artifact](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14743171356) |
| 5 | 22 | 22 | 22 | 41 | 0 | 0 | [Artifact](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14741045243) |
| 7 | 30 | 30 | 30 | 907 | 2 | 0 | [Artifact](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14740169910) |
| 9 | 31 | 31 | 30 | 653 | 2 | 0 | [Artifact](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14741374743) |
| 11 | 32 | 32 | 31 | 658 | 0 | 1 | [Artifact](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14741114551) |
| 13 | 33 | 33 | 32 | 529 | 2 | 3 | [Artifact](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14742945368) |
| 15 | 33 | 33 | 31 | 85 | 2 | 48 | [Artifact](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14742179009) |
| 17 | 33 | 33 | 32 | 80 | 2 | 35 | [Artifact](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14743130536) |
| 19 | 33 | 33 | 30 | 28 | 3 | 32 | [Artifact](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14743042773) |
| 21 | 33 | 33 | 32 | 406 | 4 | 2 | [Artifact](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14742989228) |
| 23 | 34 | 34 | 33 | 418 | 0 | 1 | [Artifact](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14742858255) |
| 25 | 33 | 33 | 32 | 396 | 1 | 1 | [Artifact](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14744645818) |
| 27 | 34 | 34 | 31 | 91 | 3 | 5 | [Artifact](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14744920691) |
| 29 | 34 | 34 | 32 | 167 | 4 | 4 | [Artifact](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14744032141) |
| 31 | 33 | 33 | 33 | 130 | 5 | 0 | [Artifact](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14744596670) |
| 35 | 34 | 34 | 31 | 155 | 2 | 3 | [Artifact](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14745163228) |

## Retained real triggers

### ObsID 26655 · acisf26655N002_evt2.fits.gz · real/527

ICRS 350.8694519°, +58.8471507° (r90 1.27″). Peak MJD 59943.23930625 UTC. Photon t95−t05 50.3 s in a 60 s trial. Peak aperture flux 5.38e-12 ± 1.3e-12 erg cm⁻² s⁻¹ (0.5–7 keV, counting uncertainty only).

Frozen pointing target `CXOUJ232327.8+584842`, proposal 24500364, “Monitoring the thermal evolution of the Cas A CCO” (GO).

The window includes 31 aperture photons and 58 annular photons. Before and after the window the aperture holds 834 photons in 7982 s and 2403 in 22056 s. The trigger spans 27 frames and 28 detector pixels. Frozen classification **ambiguous**.

The extraction aperture encloses a persistent bright knot embedded in structured shell emission at Galactic latitude −2.10°. The source and scaled annulus show large repeated excursions together throughout the visit, including near-zero intervals. The 60 s selection falls on one of those excursions. The detector positions move with the sky, which is consistent with the knot. In this modulated field the trigger supplies little evidence for an isolated flare. No catalog counterpart was found.

Input `acisf26655N002_evt2.fits.gz`, 97,667,315 bytes, SHA-256 `fa1553d8f20a7ab2bd79a30c7596b3a8623f447f76841b73c1618781f8c17852`. Detector plot SHA-256 `2a2fdb22f13d1f7bb0cf867630d52ea006aa5f357b12149a3373ea4d901590d2`, artifact member `photometry/26655/acisf26655N002_evt2.fits.gz/real/527/diagnostics.png`. [Worker archive](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14743171356) · [Archive event product](https://cxc.cfa.harvard.edu/cdaftp/byobsid/5/26655/primary/acisf26655N002_evt2.fits.gz)

### ObsID 26655 · acisf26655N002_evt2.fits.gz · real/558

ICRS 350.8628782°, +58.8515520° (r90 1.41″). Peak MJD 59943.21742514 UTC. Photon t95−t05 54.8 s in a 60 s trial. Peak aperture flux 1.95e-12 ± 6.7e-13 erg cm⁻² s⁻¹ (0.5–7 keV, counting uncertainty only).

Frozen pointing target `CXOUJ232327.8+584842`, proposal 24500364, “Monitoring the thermal evolution of the Cas A CCO” (GO).

The window includes 13 aperture photons and 21 annular photons. Before and after the window the aperture holds 134 photons in 6092 s and 582 in 23946 s. The trigger spans 11 frames and 13 detector pixels. Frozen classification **ambiguous**.

The source and scaled annulus rise and fall together repeatedly near this part of the resolved shell. The 60 s peak falls on one of the joint maxima. The independent readout-strip statistic has a small local probability (9.4 × 10⁻⁹) even though the combined streak criteria were not met. These background and readout patterns merit caution. No catalog counterpart was identified at Galactic latitude −2.10°.

Input `acisf26655N002_evt2.fits.gz`, 97,667,315 bytes, SHA-256 `fa1553d8f20a7ab2bd79a30c7596b3a8623f447f76841b73c1618781f8c17852`. Detector plot SHA-256 `371c955af91e4a5a7376ddd7596610c0f381c7f7778068df420c161d2e2ad71b`, artifact member `photometry/26655/acisf26655N002_evt2.fits.gz/real/558/diagnostics.png`. [Worker archive](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14743171356) · [Archive event product](https://cxc.cfa.harvard.edu/cdaftp/byobsid/5/26655/primary/acisf26655N002_evt2.fits.gz)

### ObsID 26655 · acisf26655N002_evt2.fits.gz · real/908

ICRS 350.8508896°, +58.7910354° (r90 1.27″). Peak MJD 59943.41540074 UTC. Photon t95−t05 455.9 s in a 600 s trial. Peak aperture flux 1.56e-12 ± 7.3e-13 erg cm⁻² s⁻¹ (0.5–7 keV, counting uncertainty only).

Frozen pointing target `CXOUJ232327.8+584842`, proposal 24500364, “Monitoring the thermal evolution of the Cas A CCO” (GO).

The window includes 38 aperture photons and 119 annular photons. Before and after the window the aperture holds 345 photons in 22746 s and 101 in 6752 s. The trigger spans 36 frames and 36 detector pixels. Frozen classification **ambiguous**.

The local image shows a curved filament running through the source aperture and into the annulus. The trigger photons follow this structure, and their detector positions track the sky dither over the 600 s window. Comparable earlier source-count peaks appear during the same visit. Diffuse shell emission and varying detector coverage are plausible explanations for this low-latitude trigger. No catalog counterpart was returned.

Input `acisf26655N002_evt2.fits.gz`, 97,667,315 bytes, SHA-256 `fa1553d8f20a7ab2bd79a30c7596b3a8623f447f76841b73c1618781f8c17852`. Detector plot SHA-256 `b70d83bf45bd992216caf56d079ed16fc65b6b435d06bb533f21e5a108ea78e8`, artifact member `photometry/26655/acisf26655N002_evt2.fits.gz/real/908/diagnostics.png`. [Worker archive](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14743171356) · [Archive event product](https://cxc.cfa.harvard.edu/cdaftp/byobsid/5/26655/primary/acisf26655N002_evt2.fits.gz)

### ObsID 26655 · acisf26655N002_evt2.fits.gz · real/922

ICRS 350.9289323°, +58.7891337° (r90 1.29″). Peak MJD 59943.28294867 UTC. Photon t95−t05 501.6 s in a 600 s trial. Peak aperture flux 8.47e-13 ± 4.5e-13 erg cm⁻² s⁻¹ (0.5–7 keV, counting uncertainty only).

Frozen pointing target `CXOUJ232327.8+584842`, proposal 24500364, “Monitoring the thermal evolution of the Cas A CCO” (GO).

The window includes 21 aperture photons and 71 annular photons. Before and after the window the aperture holds 84 photons in 11374 s and 145 in 18124 s. The trigger spans 21 frames and 21 detector pixels. Frozen classification **ambiguous**.

The aperture lies near the shell's projected edge, with a substantial brightness gradient across its annulus. The CCD-wide rate increases across the event, and the window's CCD-background probability is 2.1 × 10⁻⁴¹. Both detector-coordinate traces follow sky motion, and the full light curve shows other excursions of similar size. The IRAS source 77″ away has a chance-projection upper estimate of 1.0. The existing evidence favors a background or exposure-related excursion over an isolated extragalactic flare.

Input `acisf26655N002_evt2.fits.gz`, 97,667,315 bytes, SHA-256 `fa1553d8f20a7ab2bd79a30c7596b3a8623f447f76841b73c1618781f8c17852`. Detector plot SHA-256 `fd672e840c928a6cab32086810725c342dd0e711a846c0b0d95917dac573b2f5`, artifact member `photometry/26655/acisf26655N002_evt2.fits.gz/real/922/diagnostics.png`. [Worker archive](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14743171356) · [Archive event product](https://cxc.cfa.harvard.edu/cdaftp/byobsid/5/26655/primary/acisf26655N002_evt2.fits.gz)

### ObsID 30306 · acisf30306N001_evt2.fits.gz · real/13

ICRS 49.8788265°, +41.4622643° (r90 2.33″). Peak MJD 61088.04219648 UTC. Photon t95−t05 162.9 s in a 180 s trial. Peak aperture flux 8.2e-13 ± 4.5e-13 erg cm⁻² s⁻¹ (0.5–7 keV, counting uncertainty only).

Frozen pointing target `PERSEUSCLUSTER`, proposal 25808977, “The sounds of feedback: deep and wide imaging of the cool core of the Perseus cluster” (CLP).

The window includes 10 aperture photons and 15 annular photons. Before and after the window the aperture holds 47 photons in 6236 s and 4 in 312 s. The trigger spans 10 frames and 10 detector pixels. Frozen classification **ambiguous**.

The CCD image shows a strong diffuse brightness gradient. The trigger photons are scattered across the inner aperture, and there are recurrent aperture excursions elsewhere in the exposure. The changing detector row is consistent with the field-derived dither. Only 312 s of the good-time interval follow this window. NED lists WISEA J031930.69+412753.5 at 9.96″, beyond the 3.38″ combined match radius, with no catalog redshift and a chance-projection upper estimate of 0.42. Diffuse emission and Poisson clustering remain plausible. The galaxy cannot support an extragalactic host identification.

Input `acisf30306N001_evt2.fits.gz`, 63,414,760 bytes, SHA-256 `238fa5522b85c49adf0b8340642052c46812d04bad564f298530fe36c163ca4a`. Detector plot SHA-256 `33bd4e516978edb92b3b9f4fbcd183381f3761e07e6e714229e78714e6fa898c`, artifact member `photometry/30306/acisf30306N001_evt2.fits.gz/real/13/diagnostics.png`. [Worker archive](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14740169910) · [Archive event product](https://cxc.cfa.harvard.edu/cdaftp/byobsid/6/30306/primary/acisf30306N001_evt2.fits.gz)

### ObsID 30306 · acisf30306N001_evt2.fits.gz · real/842

ICRS 50.0042349°, +41.5492216° (r90 2.76″). Peak MJD 61088.52941398 UTC. Photon t95−t05 46.6 s in a 60 s trial. Peak aperture flux 2.36e-12 ± 9.1e-13 erg cm⁻² s⁻¹ (0.5–7 keV, counting uncertainty only).

Frozen pointing target `PERSEUSCLUSTER`, proposal 25808977, “The sounds of feedback: deep and wide imaging of the cool core of the Perseus cluster” (CLP).

The window includes 8 aperture photons and 7 annular photons. Before and after the window the aperture holds 611 photons in 48394 s and 19 in 1632 s. The trigger spans 7 frames and 8 detector pixels. Frozen classification **no_counterpart**.

The aperture lies on a pronounced gradient in bright extended emission. Its trigger photons are distributed across the aperture. The full-exposure source and scaled-annulus light curves track similar count rates, with other source excursions comparable in the displayed binning. The detector positions span several rows and columns. No catalog object meets the matching criteria. A brief fluctuation in diffuse emission is a competing explanation, and this record supplies no host identification.

Input `acisf30306N001_evt2.fits.gz`, 63,414,760 bytes, SHA-256 `238fa5522b85c49adf0b8340642052c46812d04bad564f298530fe36c163ca4a`. Detector plot SHA-256 `e130165402706567313c2579063e4f624de278b93719cd4c0c9ea9787efdc086`, artifact member `photometry/30306/acisf30306N001_evt2.fits.gz/real/842/diagnostics.png`. [Worker archive](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14740169910) · [Archive event product](https://cxc.cfa.harvard.edu/cdaftp/byobsid/6/30306/primary/acisf30306N001_evt2.fits.gz)

### ObsID 30316 · acisf30316N001_evt2.fits.gz · real/48

ICRS 49.9133943°, +41.4813794° (r90 2.33″). Peak MJD 61076.57899357 UTC. Photon t95−t05 54.8 s in a 60 s trial. Peak aperture flux unavailable (0.5–7 keV, counting uncertainty only).

Frozen pointing target `PERSEUSCLUSTER`, proposal 25808977, “The sounds of feedback: deep and wide imaging of the cool core of the Perseus cluster” (CLP).

The window includes 12 aperture photons and 18 annular photons. Before and after the window the aperture holds 536 photons in 22378 s and 503 in 22678 s. The trigger spans 10 frames and 12 detector pixels. Frozen classification **ambiguous**.

The full-field image places the aperture in bright, spatially varying extended emission. Twelve photons are dispersed across the aperture in 60 s. The full-exposure aperture curve fluctuates around the scaled annulus rate. NED lists unclassified SSTSL2 J031939.46+412854.2 at 3.12″ with a chance-projection upper estimate of 0.24. The much less localized 2MAXI object is 76″ away with a 936″ catalog major-axis uncertainty. Neither object identifies a host. A fluctuation in the extended emission remains a plausible explanation.

Input `acisf30316N001_evt2.fits.gz`, 53,917,648 bytes, SHA-256 `fc0055d9d94a93e4b6a6b137b88958fc369cbdf665f3a62dd0c6cae44dc58694`. Detector plot SHA-256 `33af41023215b38b4890aad64b5705b1eccc2c93e459a8f367aa8c189591253e`, artifact member `photometry/30316/acisf30316N001_evt2.fits.gz/real/48/diagnostics.png`. [Worker archive](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14741374743) · [Archive event product](https://cxc.cfa.harvard.edu/cdaftp/byobsid/6/30316/primary/acisf30316N001_evt2.fits.gz)

### ObsID 30316 · acisf30316N001_evt2.fits.gz · real/207

ICRS 49.8943154°, +41.5374291° (r90 1.40″). Peak MJD 61076.69532658 UTC. Photon t95−t05 480.6 s in a 600 s trial. Peak aperture flux 6.9e-13 ± 4.8e-13 erg cm⁻² s⁻¹ (0.5–7 keV, counting uncertainty only).

Frozen pointing target `PERSEUSCLUSTER`, proposal 25808977, “The sounds of feedback: deep and wide imaging of the cool core of the Perseus cluster” (CLP).

The window includes 11 aperture photons and 10 annular photons. Before and after the window the aperture holds 56 photons in 32429 s and 34 in 12087 s. The trigger spans 10 frames and 11 detector pixels. Frozen classification **no_counterpart**.

The light curve shows a concentrated peak at the start of the 600 s window, followed by sparse photons. Both detector-coordinate traces closely track the independently estimated sky dither. The CCD-wide background shows no simultaneous increase. The aperture nevertheless sits on a bright diffuse gradient, with trigger photons extending nearly to its boundary. A sky-localized brightening and a fluctuation of the local emission are both viable. The catalog checks returned no positional counterpart or host.

Input `acisf30316N001_evt2.fits.gz`, 53,917,648 bytes, SHA-256 `fc0055d9d94a93e4b6a6b137b88958fc369cbdf665f3a62dd0c6cae44dc58694`. Detector plot SHA-256 `e434ee6e591579401960eb26891c7202a6eee97c6dedad7b669c1d7cd7390ab3`, artifact member `photometry/30316/acisf30316N001_evt2.fits.gz/real/207/diagnostics.png`. [Worker archive](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14741374743) · [Archive event product](https://cxc.cfa.harvard.edu/cdaftp/byobsid/6/30316/primary/acisf30316N001_evt2.fits.gz)

### ObsID 30298 · acisf30298N001_evt2.fits.gz · real/89

ICRS 49.9330583°, +41.5070342° (r90 1.30″). Peak MJD 61083.77470930 UTC. Photon t95−t05 442.9 s in a 600 s trial. Peak aperture flux 8e-13 ± 5.5e-13 erg cm⁻² s⁻¹ (0.5–7 keV, counting uncertainty only).

Frozen pointing target `PERSEUSCLUSTER`, proposal 25808977, “The sounds of feedback: deep and wide imaging of the cool core of the Perseus cluster” (CLP).

The window includes 21 aperture photons and 58 annular photons. Before and after the window the aperture holds 129 photons in 26144 s and 64 in 9850 s. The trigger spans 21 frames and 21 detector pixels. Frozen classification **ambiguous**.

The event lies near a bright corner of the CCD and within extended emission. The aperture photons follow the sky dither, while the CCD background also rises across the event, from roughly 6–7 to 9–10 counts/s in the diagnostic. The ledger gives a CCD-background probability of 3.6 × 10⁻²². The instrument pass therefore coexists with measurable chip-wide variability. SN 1968A and 2MAXI J0319+415 have offsets of 58″ and 93″ and chance upper estimates of 1.0. Neither supports identification of this event. Spatially varying and time-dependent background compromise the interpretation.

Input `acisf30298N001_evt2.fits.gz`, 43,680,977 bytes, SHA-256 `54999d8f290441c3833a9e488a85b459e8c1e8ed3cc30adc923fd33e86ffe344`. Detector plot SHA-256 `28df87aa36e38f8ad09dafcdda5cdb029008671979d9555a7bb91c21f01e27cf`, artifact member `photometry/30298/acisf30298N001_evt2.fits.gz/real/89/diagnostics.png`. [Worker archive](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14742945368) · [Archive event product](https://cxc.cfa.harvard.edu/cdaftp/byobsid/8/30298/primary/acisf30298N001_evt2.fits.gz)

### ObsID 30298 · acisf30298N001_evt2.fits.gz · real/159

ICRS 49.8948464°, +41.5248949° (r90 2.37″). Peak MJD 61083.88475366 UTC. Photon t95−t05 50.3 s in a 60 s trial. Peak aperture flux 2.34e-12 ± 8.4e-13 erg cm⁻² s⁻¹ (0.5–7 keV, counting uncertainty only).

Frozen pointing target `PERSEUSCLUSTER`, proposal 25808977, “The sounds of feedback: deep and wide imaging of the cool core of the Perseus cluster” (CLP).

The window includes 11 aperture photons and 17 annular photons. Before and after the window the aperture holds 824 photons in 35775 s and 19 in 760 s. The trigger spans 7 frames and 11 detector pixels. Frozen classification **ambiguous**.

The sky aperture samples a steep brightness gradient, with photons dispersed across its face. The long light curve has several comparable excursions at the displayed binning. A large fraction of the 60 s increase falls near its end, and only 760 s of this good-time interval follow it. NED lists galaxy PCC 4604 at 8.15″ with no redshift and chance-projection upper estimate 0.36. The broad localization of 2MAXI J0319+415 also overlaps the field. No secure host supports the flare interpretation, and clustering in bright diffuse emission remains plausible.

Input `acisf30298N001_evt2.fits.gz`, 43,680,977 bytes, SHA-256 `54999d8f290441c3833a9e488a85b459e8c1e8ed3cc30adc923fd33e86ffe344`. Detector plot SHA-256 `97c98029eb633eacf70197db07dd1162aaa27208a8e3d3e19c0e7f8d08e90003`, artifact member `photometry/30298/acisf30298N001_evt2.fits.gz/real/159/diagnostics.png`. [Worker archive](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14742945368) · [Archive event product](https://cxc.cfa.harvard.edu/cdaftp/byobsid/8/30298/primary/acisf30298N001_evt2.fits.gz)

### ObsID 24264 · acisf24264N002_evt2.fits.gz · real/0

ICRS 80.1201904°, -69.0468915° (r90 2.31″). Peak MJD 59369.49482369 UTC. Photon t95−t05 1648.4 s in a 1800 s trial. Peak aperture flux 2.67e-12 ± 1.9e-12 erg cm⁻² s⁻¹ (0.5–7 keV, counting uncertainty only).

Frozen pointing target `SNR0519-69.0`, proposal 22500449, “Deep Observations of the LMC SNR 0519-69.0” (GO).

The window includes 9 aperture photons and 6 annular photons. Before and after the window the aperture holds 1 photons in 2928 s and 4 in 20339 s. The trigger spans 8 frames and 9 detector pixels. Frozen classification **ambiguous**.

The sky image shows a loose concentration of photons near the aperture center, with a small pair toward its edge. The trigger occupies several phases of the detector dither. The distant bright CCD source lies well off the plotted readout line, and the CCD rate stays near its surrounding level. The aperture brightening is apparent against sparse flanks. Five Gaia entries lie within 2.58″. The nearest has a 1.07″ offset and a second entry 1.13″ away has stellar astrometry. Their local chance-projection upper estimates are approximately 0.95 in this crowded field. The flare is a plausible sky event whose counterpart and distance remain undetermined.

Input `acisf24264N002_evt2.fits.gz`, 7,282,341 bytes, SHA-256 `610e20209eda2a0342f19db41588f5a3c5d95e1bdea433325aaab5d4811c13fd`. Detector plot SHA-256 `c97e98b75c07e54a31fbf0d024d5138540f6ad5e81c36cba381b56297c5e8aea`, artifact member `photometry/24264/acisf24264N002_evt2.fits.gz/real/0/diagnostics.png`. [Worker archive](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14742179009) · [Archive event product](https://cxc.cfa.harvard.edu/cdaftp/byobsid/4/24264/primary/acisf24264N002_evt2.fits.gz)

### ObsID 28436 · acisf28436N002_evt2.fits.gz · real/0

ICRS 16.0000230°, -72.0297859° (r90 1.34″). Peak MJD 60432.20518324 UTC. Photon t95−t05 437.2 s in a 600 s trial. Peak aperture flux 1.17e-11 ± 9.6e-12 erg cm⁻² s⁻¹ (0.5–7 keV, counting uncertainty only).

Frozen pointing target `E0102-72.3`, proposal 25500475, “AO-25 Calibration Observations of E0102-72” (CAL).

The window includes 18 aperture photons and 43 annular photons. Before and after the window the aperture holds 62 photons in 10751 s and 137 in 23714 s. The trigger spans 18 frames and 18 detector pixels. Frozen classification **ambiguous**.

The aperture and its annulus intersect a bright, curved rim. SIMBAD places supernova remnant SNR B0102−72.3 at 7.51″, consistent with the resolved shell morphology. The brief source-aperture rise and photons along the field-derived dither trace are visible. Spatially nonuniform remnant emission limits background subtraction across the rim. NED's typed star at 1.34″ and an unclassified Gaia source at 1.25″ add alternatives, with local chance-projection upper estimates of 0.58 and 0.36. There is no secure flare counterpart or galaxy host.

Input `acisf28436N002_evt2.fits.gz`, 591,672 bytes, SHA-256 `c1f05f214c3e517e21fe0ad3f73c569a2e37ed1323b8c62b0b6a75ce07d30138`. Detector plot SHA-256 `10df88380c5bdef1ec9f9dafedabfe3b28a5604bd5cf6ed9c92aa548ca7f79c4`, artifact member `photometry/28436/acisf28436N002_evt2.fits.gz/real/0/diagnostics.png`. [Worker archive](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14742179009) · [Archive event product](https://cxc.cfa.harvard.edu/cdaftp/byobsid/6/28436/primary/acisf28436N002_evt2.fits.gz)

### ObsID 23538 · acisf23538N002_evt2.fits.gz · real/14

ICRS 6.2261908°, +64.1691555° (r90 1.41″). Peak MJD 59484.64451377 UTC. Photon t95−t05 153.0 s in a 180 s trial. Peak aperture flux 1.41e-12 ± 6.1e-13 erg cm⁻² s⁻¹ (0.5–7 keV, counting uncertainty only).

Frozen pointing target `TYCHO'SSNR`, proposal 22500080, “Tycho's SNR: A New Era” (GO).

The window includes 12 aperture photons and 22 annular photons. Before and after the window the aperture holds 180 photons in 26875 s and 11 in 2026 s. The trigger spans 10 frames and 12 detector pixels. Frozen classification **ambiguous**.

The 180 s pulse stands above the aperture baseline and its detector coordinates broadly follow the sky dither. The full-field image places the aperture inside a bright, resolved shell at Galactic latitude +1.45°. The local image contains extended emission through the extraction annulus, and the CCD-background probability is 3.2 × 10⁻⁵. Its amplitude did not trigger the pipeline veto. No catalog entry supports a host. A fluctuation of shell emission or a Galactic source remains possible.

Input `acisf23538N002_evt2.fits.gz`, 40,469,177 bytes, SHA-256 `e7d507093cb515d1114c3e1e4c1a487e6682a03df89714c41aec97e35d54e03d`. Detector plot SHA-256 `62422ebc7f1cf96844bfce1bb161569362bda9e0ed9b998e4d7842152dd0d761`, artifact member `photometry/23538/acisf23538N002_evt2.fits.gz/real/14/diagnostics.png`. [Worker archive](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14743130536) · [Archive event product](https://cxc.cfa.harvard.edu/cdaftp/byobsid/8/23538/primary/acisf23538N002_evt2.fits.gz)

### ObsID 23538 · acisf23538N002_evt2.fits.gz · real/19

ICRS 6.2398146°, +64.1701956° (r90 2.00″). Peak MJD 59484.53817675 UTC. Photon t95−t05 59.7 s in a 60 s trial. Peak aperture flux 3.91e-12 ± 1.2e-12 erg cm⁻² s⁻¹ (0.5–7 keV, counting uncertainty only).

Frozen pointing target `TYCHO'SSNR`, proposal 22500080, “Tycho's SNR: A New Era” (GO).

The window includes 23 aperture photons and 66 annular photons. Before and after the window the aperture holds 1684 photons in 17688 s and 1075 in 11333 s. The trigger spans 12 frames and 22 detector pixels. Frozen classification **ambiguous**.

At a second position in the same shell, the aperture straddles bright knots and gradients. The scaled annulus remains close to the source-aperture count rate across the visit, and the source curve shows several comparable excursions at its displayed resolution. The 23 photons in the 60 s window have a similar soft energy distribution to the surrounding annulus photons. Baseline emission is substantial on both sides. The catalog cone yields no positional counterpart at this low Galactic latitude. Shell emission and counting fluctuations provide a plausible explanation.

Input `acisf23538N002_evt2.fits.gz`, 40,469,177 bytes, SHA-256 `e7d507093cb515d1114c3e1e4c1a487e6682a03df89714c41aec97e35d54e03d`. Detector plot SHA-256 `53e23393a763146f25b30cc72d798017ab24fc616cf387445db509963969db33`, artifact member `photometry/23538/acisf23538N002_evt2.fits.gz/real/19/diagnostics.png`. [Worker archive](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14743130536) · [Archive event product](https://cxc.cfa.harvard.edu/cdaftp/byobsid/8/23538/primary/acisf23538N002_evt2.fits.gz)

### ObsID 24391 · acisf24391N002_evt2.fits.gz · real/19

ICRS 6.1985094°, +64.1311561° (r90 2.35″). Peak MJD 59363.34072088 UTC. Photon t95−t05 56.5 s in a 60 s trial. Peak aperture flux -2.39e-12 ± 6e-12 erg cm⁻² s⁻¹ (0.5–7 keV, counting uncertainty only).

Frozen pointing target `TYCHOSNR`, proposal 22500103, “Measuring Ion Temperatures in the Ejecta of the Tycho Supernova Remnant” (GO).

The window includes 18 aperture photons and 48 annular photons. Before and after the window the aperture holds 1499 photons in 22746 s and 400 in 5940 s. The trigger spans 12 frames and 18 detector pixels. Frozen classification **ambiguous**.

The TYCHOSNR pointing places this aperture across extended shell structure. The full-visit aperture and annulus have comparable rates and multiple excursions. An unclassified Gaia entry 3.08″ away has a 0.20 chance-projection upper estimate. The readout-strip probability is small (3.0 × 10⁻⁷), although the combined streak veto did not fire. Shell emission, the local brightness gradient and readout geometry limit interpretation of this low-latitude trigger.

Input `acisf24391N002_evt2.fits.gz`, 39,702,027 bytes, SHA-256 `408e305c62c15cccebe361f17452546b5907273eb7c3527155321fce7ba31910`. Detector plot SHA-256 `2fba6dc7d793d14b5c63fa121b1eada8434794665598b50f671368cf0fb0f89a`, artifact member `photometry/24391/acisf24391N002_evt2.fits.gz/real/19/diagnostics.png`. [Worker archive](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14743042773) · [Archive event product](https://cxc.cfa.harvard.edu/cdaftp/byobsid/1/24391/primary/acisf24391N002_evt2.fits.gz)

### ObsID 24853 · acisf24853N002_evt2.fits.gz · real/0

ICRS 194.1389757°, +27.4461816° (r90 2.07″). Peak MJD 59156.73566147 UTC. Photon t95−t05 1097.3 s in a 1800 s trial. Peak aperture flux 5.6e-13 ± 3.4e-13 erg cm⁻² s⁻¹ (0.5–7 keV, counting uncertainty only).

Frozen pointing target `NGC4839`, proposal 21800392, “A detailed view of the group NGC 4839 as it merges with the Coma cluster” (GO).

The window includes 10 aperture photons and 8 annular photons. Before and after the window the aperture holds 1 photons in 3553 s and 1 in 5592 s. The trigger spans 10 frames and 10 detector pixels. Frozen classification **ambiguous**.

A compact 1,800 s sky brightening at Galactic latitude +88.82° is flanked by one aperture photon on each side, with no corresponding CCD-wide increase. Gaia DR3 1463673349229204736 is 1.08″ from the epoch-propagated position, with parallax 5.006 ± 0.150 mas and RA proper motion −62.189 ± 0.149 mas yr⁻¹. Its local chance-projection upper estimate is 0.01074, just above the automated veto threshold. NED also identifies a star 1.33″ away. A stellar flare is the leading explanation. A galaxy 7.88″ away, without a redshift and with chance upper estimate 0.34, offers weaker host evidence.

Input `acisf24853N002_evt2.fits.gz`, 3,946,486 bytes, SHA-256 `f2b7ed4bfb4253ae4276023c7b1544f97d96d03f92928919856e33524d3b977c`. Detector plot SHA-256 `b46e0614839b5a20b0b4875c86ddf181e638c6ede9cb9d4398f4860cb98f3ed9`, artifact member `photometry/24853/acisf24853N002_evt2.fits.gz/real/0/diagnostics.png`. [Worker archive](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14743042773) · [Archive event product](https://cxc.cfa.harvard.edu/cdaftp/byobsid/3/24853/primary/acisf24853N002_evt2.fits.gz)

### ObsID 32380 · acisf32380N001_evt2.fits.gz · real/0

ICRS 243.5485946°, -14.3596881° (r90 1.24″). Peak MJD 61199.79665103 UTC. Photon t95−t05 162.4 s in a 180 s trial. Peak aperture flux 4.04e-12 ± 1.1e-12 erg cm⁻² s⁻¹ (0.5–7 keV, counting uncertainty only).

Frozen pointing target `2RXSJ161409.7-142143`, proposal 20201077, “Stellar activity with TESS and Chandra” (CCT).

The window includes 37 aperture photons and 1 annular photons. Before and after the window the aperture holds 15 photons in 4179 s and 122 in 5724 s. The trigger spans 28 frames and 32 detector pixels. Frozen classification **ambiguous**.

**Priority measured flare in a targeted pointing.** A compact group of 37 photons follows the independent field dither on both detector axes. Only one annular photon is present, and the CCD-wide rate shows no comparable rise. The source climbs abruptly and remains above its earlier baseline after the peak. The post-peak aperture collected 122 photons in 5,724 s, compared with 15 in 4,179 s before the window. That sustained emission limits the interpretation as an isolated outburst.

There is no Gaia or SIMBAD match. NED's only entry, 1RXS J161409.6-142147, is at RA 243.540000°, Dec −14.363056°, 32.3″ away with a 35″ catalog uncertainty. The frozen inventory lists proposal 20201077, PI Günther, “Stellar activity with TESS and Chandra,” category STARS AND WD, type CCT. Chandra observed target 2RXSJ161409.7-142143 on 2026-06-08 and released the data on 2026-06-11. Its finer FITS pointing target at RA 243.547897°, Dec −14.358500° is within about 5″ of the Chandra centroid. The trigger therefore provides evidence for rapid variation in the direction of a known, deliberately observed ROSAT source. The program specifically targeted stellar activity, making a stellar flare a strong alternative. Its distance and physical class remain undetermined. Independent visits and contemporary event notices bear on that interpretation. [Frozen pointing inventory](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14739228660)

Input `acisf32380N001_evt2.fits.gz`, 1,081,909 bytes, SHA-256 `3c5f46567187114d7deb741d2a709bbc4e3566801a90a9cdb3f23a7c66a944cf`. Detector plot SHA-256 `56acdc764c60f2d692835b7f74d9be6190a63d8cbc990a4fa71e72992df68909`, artifact member `photometry/32380/acisf32380N001_evt2.fits.gz/real/0/diagnostics.png`. [Worker archive](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14743042773) · [Archive event product](https://cxc.cfa.harvard.edu/cdaftp/byobsid/0/32380/primary/acisf32380N001_evt2.fits.gz)

### ObsID 24363 · acisf24363N002_evt2.fits.gz · real/0

ICRS 266.4170986°, -29.0078101° (r90 1.35″). Peak MJD 59756.31725780 UTC. Photon t95−t05 167.4 s in a 180 s trial. Peak aperture flux 1.65e-12 ± 8.8e-13 erg cm⁻² s⁻¹ (0.5–7 keV, counting uncertainty only).

Frozen pointing target `CMZMOLECULARCLOUD`, proposal 22910118, “Unique X-ray Window into Internal Structure of the Cold Gas in Dense Molecular Clouds” (GO).

The window includes 10 aperture photons and 7 annular photons. Before and after the window the aperture holds 12 photons in 1825 s and 123 in 26078 s. The trigger spans 10 frames and 9 detector pixels. Frozen classification **ambiguous**.

The source photons move with the expected dither in both axes. They lie within the central Galactic field, where the ledger returns 198 positional matches. SIMBAD lists the Sgr A West Northern Arm 0.13″ away and the variable star [GKM98] S1-1 at 0.20″, alongside numerous other stars. The aperture continues to receive photons outside the trigger. The spatial crowding and Galactic foreground prevent a unique physical interpretation.

Input `acisf24363N002_evt2.fits.gz`, 9,220,713 bytes, SHA-256 `06472e758aa92b55dea924d4fd7d692b86e9ca1860455bdd2065754a52ff0410`. Detector plot SHA-256 `235895110ba62fe4df59661d6a13cbfb4038f08db4ad8abd9704fb6e3db4dc99`, artifact member `photometry/24363/acisf24363N002_evt2.fits.gz/real/0/diagnostics.png`. [Worker archive](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14742989228) · [Archive event product](https://cxc.cfa.harvard.edu/cdaftp/byobsid/3/24363/primary/acisf24363N002_evt2.fits.gz)

### ObsID 24363 · acisf24363N002_evt2.fits.gz · real/1

ICRS 266.3826872°, -28.9728212° (r90 1.89″). Peak MJD 59756.52069790 UTC. Photon t95−t05 1626.6 s in a 1800 s trial. Peak aperture flux 6.12e-13 ± 3.7e-13 erg cm⁻² s⁻¹ (0.5–7 keV, counting uncertainty only).

Frozen pointing target `CMZMOLECULARCLOUD`, proposal 22910118, “Unique X-ray Window into Internal Structure of the Cold Gas in Dense Molecular Clouds” (GO).

The window includes 14 aperture photons and 25 annular photons. Before and after the window the aperture holds 26 photons in 18884 s and 12 in 7399 s. The trigger spans 14 frames and 14 detector pixels. Frozen classification **ambiguous**.

The photons cluster around a persistent source, with multiple separated peaks during the 1,800 s window and steady CCD background. Its position lies 1.48″ from the cataloged Chandra source 2CXO J174531.7-285822 and 1.60″ from SIMBAD X-ray source CXOU J174531.7-285823. An unclassified Gaia source lies 1.48″ away. Variability of a known X-ray source in this Galactic-center field provides a plausible explanation.

Input `acisf24363N002_evt2.fits.gz`, 9,220,713 bytes, SHA-256 `06472e758aa92b55dea924d4fd7d692b86e9ca1860455bdd2065754a52ff0410`. Detector plot SHA-256 `0deef01c60c9edf8a2073c0d574b253bd783eb6063b6eeadbaa41e6cf233ba15`, artifact member `photometry/24363/acisf24363N002_evt2.fits.gz/real/1/diagnostics.png`. [Worker archive](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14742989228) · [Archive event product](https://cxc.cfa.harvard.edu/cdaftp/byobsid/3/24363/primary/acisf24363N002_evt2.fits.gz)

### ObsID 24363 · acisf24363N002_evt2.fits.gz · real/2

ICRS 266.4330972°, -28.8136268° (r90 1.26″). Peak MJD 59756.52727764 UTC. Photon t95−t05 79.3 s in a 180 s trial. Peak aperture flux 1.72e-12 ± 8e-13 erg cm⁻² s⁻¹ (0.5–7 keV, counting uncertainty only).

Frozen pointing target `CMZMOLECULARCLOUD`, proposal 22910118, “Unique X-ray Window into Internal Structure of the Cold Gas in Dense Molecular Clouds” (GO).

The window includes 8 aperture photons and 1 annular photons. Before and after the window the aperture holds 41 photons in 19987 s and 41 in 7916 s. The trigger spans 7 frames and 7 detector pixels. Frozen classification **ambiguous**.

The source forms a compact image, with an early spike and photons tracing the independent dither. It lies 1.77″ from the previously cataloged Galactic-center X-ray source CXOGCS J174543.90-284849.2 and 1.70″ from a 2MASS object. The post-trigger aperture activity is greater than in the early baseline. A flare of the cataloged source is plausible. Source distance and host remain undetermined.

Input `acisf24363N002_evt2.fits.gz`, 9,220,713 bytes, SHA-256 `06472e758aa92b55dea924d4fd7d692b86e9ca1860455bdd2065754a52ff0410`. Detector plot SHA-256 `3484c3e39ab4e8aa66e36048e502a8e5860a64258a03cb5947062ed14ced2733`, artifact member `photometry/24363/acisf24363N002_evt2.fits.gz/real/2/diagnostics.png`. [Worker archive](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14742989228) · [Archive event product](https://cxc.cfa.harvard.edu/cdaftp/byobsid/3/24363/primary/acisf24363N002_evt2.fits.gz)

### ObsID 28581 · acisf28581N001_evt2.fits.gz · real/0

ICRS 181.5022854°, +52.6421366° (r90 2.17″). Peak MJD 60806.97792463 UTC. Photon t95−t05 47.0 s in a 60 s trial. Peak aperture flux 2.09e-12 ± 7.8e-13 erg cm⁻² s⁻¹ (0.5–7 keV, counting uncertainty only).

Frozen pointing target `NGC4102`, proposal 25700108, “Probing the AGN-Host Relationship in the LINER /CT-AGN NGC 4102” (GO).

The window includes 8 aperture photons and 0 annular photons. Before and after the window the aperture holds 48 photons in 5451 s and 64 in 7545 s. The trigger spans 8 frames and 8 detector pixels. Frozen classification **ambiguous**.

The eight photons are localized on an existing compact X-ray source. SIMBAD lists confirmed quasar SDSS J120600.49+523831.5 at 0.54″, z = 0.62 and chance-projection upper estimate 0.0033. NED lists WISEA J120600.51+523831.8 at 0.33″ with z = 0.6194. These positional and redshift measurements support a quasar association. The flanks contain 112 aperture photons, and other excursions occur in the full-visit light curve. The 60 s clustering warrants comparison with the survey controls and independent visits. Variability or counting noise from this known active nucleus remains a stronger interpretation than an isolated new FXT.

Input `acisf28581N001_evt2.fits.gz`, 1,809,429 bytes, SHA-256 `f5042f6cf4264987387c22edb905387f7fec2e59df80096106c92f490b38a758`. Detector plot SHA-256 `b3d03c9fa83ee8e077e8dfa93950994174bd95e4eae132ed9a9b3db9a0b19a3c`, artifact member `photometry/28581/acisf28581N001_evt2.fits.gz/real/0/diagnostics.png`. [Worker archive](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14742989228) · [Archive event product](https://cxc.cfa.harvard.edu/cdaftp/byobsid/1/28581/primary/acisf28581N001_evt2.fits.gz)

### ObsID 30303 · acisf30303N001_evt2.fits.gz · real/379

ICRS 50.0401488°, +41.5515525° (r90 2.74″). Peak MJD 61097.38557916 UTC. Photon t95−t05 157.1 s in a 180 s trial. Peak aperture flux 1.62e-12 ± 1e-12 erg cm⁻² s⁻¹ (0.5–7 keV, counting uncertainty only).

Frozen pointing target `PERSEUSCLUSTER`, proposal 25808977, “The sounds of feedback: deep and wide imaging of the cool core of the Perseus cluster” (CLP).

The window includes 11 aperture photons and 14 annular photons. Before and after the window the aperture holds 12 photons in 1629 s and 175 in 22736 s. The trigger spans 9 frames and 11 detector pixels. Frozen classification **no_counterpart**.

The CCD is dominated by extended Perseus-cluster emission. The 11 photons within the selected sky aperture have substantial spatial spread and arrive in two principal groups near the start and end of the 180 s interval. Detector coordinates scatter around the dither curve. The scaled-annulus rate tracks the long-term aperture baseline and other source excursions recur. No catalog object met the matching criteria. Local fluctuations of cluster emission are a plausible alternative to a hostless transient.

Input `acisf30303N001_evt2.fits.gz`, 36,192,851 bytes, SHA-256 `0ac82bc7d330e6e03f75ad0711259a9df2d1d4085c1299c3d502af289f192313`. Detector plot SHA-256 `d09e54e9ff0814267379566e6d3d9a6d38c202673ba97f1eaebd1fb70fd15335`, artifact member `photometry/30303/acisf30303N001_evt2.fits.gz/real/379/diagnostics.png`. [Worker archive](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14744645818) · [Archive event product](https://cxc.cfa.harvard.edu/cdaftp/byobsid/3/30303/primary/acisf30303N001_evt2.fits.gz)

### ObsID 30079 · acisf30079N002_evt2.fits.gz · real/0

ICRS 251.4794119°, -46.9617207° (r90 1.33″). Peak MJD 60575.82270147 UTC. Photon t95−t05 1302.9 s in a 1800 s trial. Peak aperture flux 7.3e-13 ± 3.7e-13 erg cm⁻² s⁻¹ (0.5–7 keV, counting uncertainty only).

Frozen pointing target `NGC6204`, proposal 24200285, “ZAMS X-rays and Primordial Planetary Atmospheres” (GO).

The window includes 18 aperture photons and 2 annular photons. Before and after the window the aperture holds 2 photons in 11554 s and 4 in 3741 s. The trigger spans 18 frames and 18 detector pixels. Frozen classification **catalog_incomplete**.

The photons are concentrated on the sky and closely follow two-axis dither. A sharp rise and weaker post-peak emission establish a convincing flare in the NGC 6204 stellar cluster. SIMBAD star UCAC4 216-133818 lies 0.63 arcseconds away, with chance upper bound 0.007. Although NED follow-up was incomplete, the identified star is a strong foreground explanation.

Input `acisf30079N002_evt2.fits.gz`, 1,154,912 bytes, SHA-256 `44787c96384558326c248fb613503360f9fd81721f324c5950f2611e65eeda9f`. Detector plot SHA-256 `d71c66cb515509d48b79d078bfd7338e7b7e1c4e4c67f922c07b2561edfea016`, artifact member `photometry/30079/acisf30079N002_evt2.fits.gz/real/0/diagnostics.png`. [Worker archive](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14744920691) · [Archive event product](https://cxc.cfa.harvard.edu/cdaftp/byobsid/9/30079/primary/acisf30079N002_evt2.fits.gz)

### ObsID 30367 · acisf30367N001_evt2.fits.gz · real/5

ICRS 49.9571804°, +41.5096091° (r90 1.41″). Peak MJD 61292.55711771 UTC. Photon t95−t05 162.4 s in a 180 s trial. Peak aperture flux 1.58e-12 ± 7e-13 erg cm⁻² s⁻¹ (0.5–7 keV, counting uncertainty only).

Frozen pointing target `PERSEUSCLUSTERE`, proposal 25808977, “The sounds of feedback: deep and wide imaging of the cool core of the Perseus cluster” (CLP).

The window includes 12 aperture photons and 17 annular photons. Before and after the window the aperture holds 50 photons in 5850 s and 256 in 26060 s. The trigger spans 10 frames and 11 detector pixels. Frozen classification **ambiguous**.

The detector track coherently follows a sky position. It lies immediately next to the diffuse, bright Perseus core, with persistent counts in both flanks and the annulus. The NGC 1275 Cluster entry 12.81 arcseconds away and SN 1968A at 24.01 arcseconds both have nearly unit chance upper bounds. Neither identifies this local fluctuation.

Input `acisf30367N001_evt2.fits.gz`, 35,851,003 bytes, SHA-256 `8791a324b9eb50fc49aea3fcf44d738cdb8546afff7e3fdd69ba4a909b0f97f0`. Detector plot SHA-256 `49e9927330799d036466ec53fc52494f2021e1586fb95fd7463a7a255e8a5509`, artifact member `photometry/30367/acisf30367N001_evt2.fits.gz/real/5/diagnostics.png`. [Worker archive](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14744920691) · [Archive event product](https://cxc.cfa.harvard.edu/cdaftp/byobsid/7/30367/primary/acisf30367N001_evt2.fits.gz)

### ObsID 30367 · acisf30367N001_evt2.fits.gz · real/18

ICRS 49.8832364°, +41.5581019° (r90 1.47″). Peak MJD 61292.73790907 UTC. Photon t95−t05 479.3 s in a 600 s trial. Peak aperture flux 1.59e-12 ± 1e-12 erg cm⁻² s⁻¹ (0.5–7 keV, counting uncertainty only).

Frozen pointing target `PERSEUSCLUSTERE`, proposal 25808977, “The sounds of feedback: deep and wide imaging of the cool core of the Perseus cluster” (CLP).

The window includes 8 aperture photons and 6 annular photons. Before and after the window the aperture holds 22 photons in 21304 s and 10 in 10186 s. The trigger spans 7 frames and 8 detector pixels. Frozen classification **ambiguous**.

The CCD rate is stable. Eight photons on a convincing two-axis dither track produce a short peak in the diffuse outskirts of Perseus, with comparable earlier aperture peaks. Galaxy PCC 4524 at 9.01 arcseconds has a chance upper bound of 0.45. Sparse diffuse-field fluctuations remain a plausible explanation.

Input `acisf30367N001_evt2.fits.gz`, 35,851,003 bytes, SHA-256 `8791a324b9eb50fc49aea3fcf44d738cdb8546afff7e3fdd69ba4a909b0f97f0`. Detector plot SHA-256 `7f822ff242cc761f99d86059e3a06e511d60b10990f78f556b60f920986422f7`, artifact member `photometry/30367/acisf30367N001_evt2.fits.gz/real/18/diagnostics.png`. [Worker archive](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14744920691) · [Archive event product](https://cxc.cfa.harvard.edu/cdaftp/byobsid/7/30367/primary/acisf30367N001_evt2.fits.gz)

### ObsID 22689 · acisf22689N002_evt2.fits.gz · real/0

ICRS 131.3121092°, -41.2322893° (r90 1.24″). Peak MJD 59216.26796986 UTC. Photon t95−t05 2978.6 s in a 3600 s trial. Peak aperture flux 8.03e-13 ± 4.6e-13 erg cm⁻² s⁻¹ (0.5–7 keV, counting uncertainty only).

Frozen pointing target `COLLINDER197`, proposal 21200606, “Evolution of Stellar Magnetic Activity in 'Missing Link' Clusters: Collinder 197” (GTO).

The window includes 59 aperture photons and 3 annular photons. Before and after the window the aperture holds 7 photons in 11306 s and 20 in 9635 s. The trigger spans 59 frames and 58 detector pixels. Frozen classification **catalog_incomplete**.

The full-visit light curve shows a well-isolated rise, structured hour-long emission and slower decline, spatially concentrated on a dithered sky position. The Collinder 197 field has a SIMBAD candidate young stellar object 2MASS J08451489-4113565 0.43 arcseconds away and a Gaia position at 0.45 arcseconds. A young-star flare is plausible despite incomplete NED follow-up.

Input `acisf22689N002_evt2.fits.gz`, 4,381,536 bytes, SHA-256 `05c76061b101f668e1cb621450bc1396d5f30cbbfb5b484b7c70534e2aea1247`. Detector plot SHA-256 `1980e9892c9d4353cffd7e2a92252d726d6c16325c6ff72cdb9ff493b868e3a9`, artifact member `photometry/22689/acisf22689N002_evt2.fits.gz/real/0/diagnostics.png`. [Worker archive](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14744032141) · [Archive event product](https://cxc.cfa.harvard.edu/cdaftp/byobsid/9/22689/primary/acisf22689N002_evt2.fits.gz)

### ObsID 30354 · acisf30354N001_evt2.fits.gz · real/50

ICRS 49.9694385°, +41.5555544° (r90 2.36″). Peak MJD 61266.91717271 UTC. Photon t95−t05 47.0 s in a 60 s trial. Peak aperture flux 2.65e-12 ± 1.2e-12 erg cm⁻² s⁻¹ (0.5–7 keV, counting uncertainty only).

Frozen pointing target `PERSEUSCLUSTERS`, proposal 25808977, “The sounds of feedback: deep and wide imaging of the cool core of the Perseus cluster” (CLP).

The window includes 8 aperture photons and 7 annular photons. Before and after the window the aperture holds 280 photons in 23970 s and 67 in 7559 s. The trigger spans 8 frames and 8 detector pixels. Frozen classification **no_counterpart**.

Eight photons occur on a persistent Perseus surface-brightness gradient. They scatter around the short dither track, the before/after aperture rates resemble the scaled annulus, and the ledger reports an anomalous CCD background statistic. No host association is recorded. A diffuse/background fluctuation remains plausible.

Input `acisf30354N001_evt2.fits.gz`, 34,391,123 bytes, SHA-256 `26f3da51434b02bc13422d1f0cdc6eb26b0eb8429681c805874d862080535d59`. Detector plot SHA-256 `6f46c873de27f7615976f802a12aeea7fa70f25224d7e4ecb00ef292949d119d`, artifact member `photometry/30354/acisf30354N001_evt2.fits.gz/real/50/diagnostics.png`. [Worker archive](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14744032141) · [Archive event product](https://cxc.cfa.harvard.edu/cdaftp/byobsid/4/30354/primary/acisf30354N001_evt2.fits.gz)

### ObsID 30354 · acisf30354N001_evt2.fits.gz · real/145

ICRS 49.9089063°, +41.4579524° (r90 2.03″). Peak MJD 61266.82181477 UTC. Photon t95−t05 153.3 s in a 180 s trial. Peak aperture flux 1.3e-12 ± 9.2e-13 erg cm⁻² s⁻¹ (0.5–7 keV, counting uncertainty only).

Frozen pointing target `PERSEUSCLUSTERS`, proposal 25808977, “The sounds of feedback: deep and wide imaging of the cool core of the Perseus cluster” (CLP).

The window includes 13 aperture photons and 18 annular photons. Before and after the window the aperture holds 177 photons in 15627 s and 166 in 15782 s. The trigger spans 11 frames and 13 detector pixels. Frozen classification **ambiguous**.

The sky events occupy a broad patch in bright cluster emission, with scatter around the detector dither trajectory. The full visit shows multiple aperture peaks and persistent diffuse counts. Galaxy PCC 4669 at 8.65 arcseconds has chance upper bound 0.38. This positional coincidence cannot securely identify a host.

Input `acisf30354N001_evt2.fits.gz`, 34,391,123 bytes, SHA-256 `26f3da51434b02bc13422d1f0cdc6eb26b0eb8429681c805874d862080535d59`. Detector plot SHA-256 `6f46d32ad996198864ea935d61e3a9b66864f4f4ac7cbb7a04043a597ef17554`, artifact member `photometry/30354/acisf30354N001_evt2.fits.gz/real/145/diagnostics.png`. [Worker archive](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14744032141) · [Archive event product](https://cxc.cfa.harvard.edu/cdaftp/byobsid/4/30354/primary/acisf30354N001_evt2.fits.gz)

### ObsID 30354 · acisf30354N001_evt2.fits.gz · real/155

ICRS 49.9388021°, +41.5029195° (r90 1.36″). Peak MJD 61266.96788583 UTC. Photon t95−t05 146.1 s in a 180 s trial. Peak aperture flux 1.16e-11 ± 1e-11 erg cm⁻² s⁻¹ (0.5–7 keV, counting uncertainty only).

Frozen pointing target `PERSEUSCLUSTERS`, proposal 25808977, “The sounds of feedback: deep and wide imaging of the cool core of the Perseus cluster” (CLP).

The window includes 11 aperture photons and 19 annular photons. Before and after the window the aperture holds 219 photons in 28267 s and 18 in 3142 s. The trigger spans 10 frames and 11 detector pixels. Frozen classification **ambiguous**.

The two detector coordinates track dither in a bright part of the Perseus core. The aperture has persistent photons and other temporal peaks. Potential hosts include PCC 4937 at 5.21 arcseconds (chance upper bound 0.18) and the catalog entries for SA 0426-042 at 9.49 arcseconds. These chance bounds and the structured cluster field leave the source identity unresolved.

Input `acisf30354N001_evt2.fits.gz`, 34,391,123 bytes, SHA-256 `26f3da51434b02bc13422d1f0cdc6eb26b0eb8429681c805874d862080535d59`. Detector plot SHA-256 `10b7b927c27d8854e567705ef54fff1b60a17f933fd5d934c992e5b4735ac10f`, artifact member `photometry/30354/acisf30354N001_evt2.fits.gz/real/155/diagnostics.png`. [Worker archive](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14744032141) · [Archive event product](https://cxc.cfa.harvard.edu/cdaftp/byobsid/4/30354/primary/acisf30354N001_evt2.fits.gz)

### ObsID 23324 · acisf23324N003_evt2.fits.gz · real/0

ICRS 127.7690303°, +9.7053843° (r90 1.43″). Peak MJD 59569.55084868 UTC. Photon t95−t05 2880.8 s in a 3600 s trial. Peak aperture flux 1.48e-12 ± 6.6e-13 erg cm⁻² s⁻¹ (0.5–7 keV, counting uncertainty only).

Frozen pointing target `LHS2021`, proposal 22200629, “Chromosopheric and Coronal Activity in the Lowest-Mass Stars” (GO).

The window includes 78 aperture photons and 61 annular photons. Before and after the window the aperture holds 6 photons in 10093 s and 19 in 6383 s. The trigger spans 75 frames and 78 detector pixels. Frozen classification **ambiguous**.

The rise, long decay and 75 separate frames form a convincing sky flare. It lies roughly 9 arcminutes from this stellar program's target. The independent flare real/1 is at LHS2021. This event's likely optical counterpart, Gaia 599875677054390784 (G=20.64), has a two-parameter position solution with neither proper motion nor parallax. Its stellar or extragalactic nature is open. ROSAT 1WGA J0831.0+0941 is 56.35 arcseconds away. A bright neighbor lies at the outer background-annulus edge. Priority serendipitous flare for independent visits and counterpart classification; the elevated decay extends into the post-window flank.

Input `acisf23324N003_evt2.fits.gz`, 5,538,954 bytes, SHA-256 `507a5aaf26b492ef2b5120242d9d33e281803667c8c1edc125fa905e6a09bbd0`. Detector plot SHA-256 `c541303e0eafad11aa98c8a76e2b2de9fea93c5e8ae112c6611ca7704aa80c0c`, artifact member `photometry/23324/acisf23324N003_evt2.fits.gz/real/0/diagnostics.png`. [Worker archive](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14744596670) · [Archive event product](https://cxc.cfa.harvard.edu/cdaftp/byobsid/4/23324/primary/acisf23324N003_evt2.fits.gz)

### ObsID 23324 · acisf23324N003_evt2.fits.gz · real/1

ICRS 127.6326220°, +9.7846623° (r90 1.23″). Peak MJD 59569.65675684 UTC. Photon t95−t05 162.1 s in a 180 s trial. Peak aperture flux 1.74e-12 ± 7.9e-13 erg cm⁻² s⁻¹ (0.5–7 keV, counting uncertainty only).

Frozen pointing target `LHS2021`, proposal 22200629, “Chromosopheric and Coronal Activity in the Lowest-Mass Stars” (GO).

The window includes 12 aperture photons and 0 annular photons. Before and after the window the aperture holds 25 photons in 19469 s and 3 in 427 s. The trigger spans 12 frames and 11 detector pixels. Frozen classification **ambiguous**.

A short soft flare follows the sky dither in both detector coordinates. This position coincides with the pointing target LHS2021. Gaia 599891555546067072, 0.60 arcseconds away at the event epoch, has parallax 59.806 ± 0.075 mas and proper motion about 670 mas/yr. The pointed nearby star is a strong foreground explanation. This is a different location from real/0.

Input `acisf23324N003_evt2.fits.gz`, 5,538,954 bytes, SHA-256 `507a5aaf26b492ef2b5120242d9d33e281803667c8c1edc125fa905e6a09bbd0`. Detector plot SHA-256 `de8e3b8415aa7af06dbcb2664a0337fca6201434a6b7d2a3063888970ecb2839`, artifact member `photometry/23324/acisf23324N003_evt2.fits.gz/real/1/diagnostics.png`. [Worker archive](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14744596670) · [Archive event product](https://cxc.cfa.harvard.edu/cdaftp/byobsid/4/23324/primary/acisf23324N003_evt2.fits.gz)

### ObsID 27801 · acisf27801N002_evt2.fits.gz · real/0

ICRS 290.3157305°, +43.9663192° (r90 2.36″). Peak MJD 60057.71255976 UTC. Photon t95−t05 169.5 s in a 180 s trial. Peak aperture flux 9.08e-13 ± 5.6e-13 erg cm⁻² s⁻¹ (0.5–7 keV, counting uncertainty only).

Frozen pointing target `A2319`, proposal 22800264, “How do Cosmic Plasmas Work? Insights from the Hottest Brightest Cold Front and Gas Perturbations in A2319” (GO).

The window includes 10 aperture photons and 14 annular photons. Before and after the window the aperture holds 85 photons in 11922 s and 33 in 5862 s. The trigger spans 10 frames and 10 detector pixels. Frozen classification **ambiguous**.

The photons partly follow the detector dither and the scaled annular rate stays low within the selected window. The sky aperture spans persistent bright, diffuse A2319 emission on a gradient, with additional short aperture peaks elsewhere in the exposure. The distant cluster and radio matches have unit local chance upper bounds. Diffuse-field statistics compete with source variability; no galaxy host can be assigned.

Input `acisf27801N002_evt2.fits.gz`, 8,226,114 bytes, SHA-256 `8a0bd9c096e515bf170ef6f0f15e343879fb2bff5fbed2654ad318e055f084df`. Detector plot SHA-256 `dcd75f1b9b26bc9f7a4184b474304df3b490bf395596e288a5132e30d97975b6`, artifact member `photometry/27801/acisf27801N002_evt2.fits.gz/real/0/diagnostics.png`. [Worker archive](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14744596670) · [Archive event product](https://cxc.cfa.harvard.edu/cdaftp/byobsid/1/27801/primary/acisf27801N002_evt2.fits.gz)

### ObsID 30332 · acisf30332N001_evt2.fits.gz · real/119

ICRS 49.9046797°, +41.5139585° (r90 1.41″). Peak MJD 61263.89848756 UTC. Photon t95−t05 141.8 s in a 180 s trial. Peak aperture flux 1.01e-12 ± 5.8e-13 erg cm⁻² s⁻¹ (0.5–7 keV, counting uncertainty only).

Frozen pointing target `PERSEUSCLUSTERN`, proposal 25808977, “The sounds of feedback: deep and wide imaging of the cool core of the Perseus cluster” (CLP).

The window includes 8 aperture photons and 6 annular photons. Before and after the window the aperture holds 69 photons in 16878 s and 44 in 12033 s. The trigger spans 7 frames and 8 detector pixels. Frozen classification **ambiguous**.

Eight events follow the two-axis dither and define an aperture spike. Both long flanks contain photons from the bright diffuse Perseus field and the peak lies in a window with elevated overall CCD statistics. The broad 2MAXI position is 51.81 arcseconds distant. Tentative local excess; the diffuse field and limited event count constrain interpretation.

Input `acisf30332N001_evt2.fits.gz`, 33,439,149 bytes, SHA-256 `3b7d751e06a906077c4701578407b17b16bd12083d40c6c52f2e7b594730fe45`. Detector plot SHA-256 `ae4dbf475bdb9dbee506382d9f63102dd8d2949ac1aa7355f80b80b66a2c801a`, artifact member `photometry/30332/acisf30332N001_evt2.fits.gz/real/119/diagnostics.png`. [Worker archive](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14744596670) · [Archive event product](https://cxc.cfa.harvard.edu/cdaftp/byobsid/2/30332/primary/acisf30332N001_evt2.fits.gz)

### ObsID 30332 · acisf30332N001_evt2.fits.gz · real/121

ICRS 49.8933815°, +41.5477224° (r90 3.12″). Peak MJD 61263.82730607 UTC. Photon t95−t05 49.9 s in a 60 s trial. Peak aperture flux 2.4e-12 ± 9.5e-13 erg cm⁻² s⁻¹ (0.5–7 keV, counting uncertainty only).

Frozen pointing target `PERSEUSCLUSTERN`, proposal 25808977, “The sounds of feedback: deep and wide imaging of the cool core of the Perseus cluster” (CLP).

The window includes 8 aperture photons and 4 annular photons. Before and after the window the aperture holds 125 photons in 10793 s and 249 in 18237 s. The trigger spans 5 frames and 8 detector pixels. Frozen classification **ambiguous**.

The sky aperture samples the bright gradient outside the Perseus core. Its eight photons in five frames span a broad patch relative to the short dither track, and another aperture spike reaches a comparable height elsewhere in the exposure. The foreground Gaia source 2.11 arcseconds away has stellar astrometry; a WISEA entry also overlaps with appreciable chance probability. Neither supplies a secure extragalactic interpretation.

Input `acisf30332N001_evt2.fits.gz`, 33,439,149 bytes, SHA-256 `3b7d751e06a906077c4701578407b17b16bd12083d40c6c52f2e7b594730fe45`. Detector plot SHA-256 `e8107cd89c09819c82695ddfb23ca691e45127242dfa2064fe2bd4410dabaee4`, artifact member `photometry/30332/acisf30332N001_evt2.fits.gz/real/121/diagnostics.png`. [Worker archive](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14744596670) · [Archive event product](https://cxc.cfa.harvard.edu/cdaftp/byobsid/2/30332/primary/acisf30332N001_evt2.fits.gz)

### ObsID 23554 · acisf23554N002_evt2.fits.gz · real/0

ICRS 56.7519829°, +43.1698841° (r90 1.42″). Peak MJD 59313.71116205 UTC. Photon t95−t05 1448.6 s in a 1800 s trial. Peak aperture flux 8.16e-13 ± 3.9e-13 erg cm⁻² s⁻¹ (0.5–7 keV, counting uncertainty only).

Frozen pointing target `AT2020WNT`, proposal 22500226, “Source of Energy and pre-explosion progenitor evolution of Super-luminous Supernovae” (TOO).

The window includes 40 aperture photons and 4 annular photons. Before and after the window the aperture holds 0 photons in 6240 s and 15 in 6431 s. The trigger spans 40 frames and 40 detector pixels. Frozen classification **catalog_incomplete**.

Forty photons follow a compact source through repeated dither cycles after an empty 6.2 ks baseline. The nearest Gaia source at 0.06 arcseconds has parallax 1.558 ± 0.129 mas and RUWE 0.97; its neighbor at 1.82 arcseconds has 1.509 ± 0.113 mas and similar proper motion. The field targeted a supernova, but these two foreground stars, possibly associated with each other, strongly favor stellar flaring. The readout-strip statistic is small; the sky image and two-axis dither strengthen the astrophysical interpretation. An all-untyped NED cone left the frozen catalog state incomplete.

Input `acisf23554N002_evt2.fits.gz`, 10,966,421 bytes, SHA-256 `16707b3d950574391e0dc106ef5fad6914509de1d805a8aa6c42035f28a96d0e`. Detector plot SHA-256 `d73834465702b79e5d4817c5622bec4a70ba50a2ded21cec1f5d4a2ecb2781ea`, artifact member `photometry/23554/acisf23554N002_evt2.fits.gz/real/0/diagnostics.png`. [Worker archive](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14745163228) · [Archive event product](https://cxc.cfa.harvard.edu/cdaftp/byobsid/4/23554/primary/acisf23554N002_evt2.fits.gz)

### ObsID 29634 · acisf29634N001_evt2.fits.gz · real/106

ICRS 49.9027179°, +41.4831035° (r90 1.33″). Peak MJD 61265.43655351 UTC. Photon t95−t05 163.3 s in a 180 s trial. Peak aperture flux 5.6e-13 ± 4e-13 erg cm⁻² s⁻¹ (0.5–7 keV, counting uncertainty only).

Frozen pointing target `PERSEUSCLUSTERS`, proposal 25808977, “The sounds of feedback: deep and wide imaging of the cool core of the Perseus cluster” (CLP).

The window includes 8 aperture photons and 7 annular photons. Before and after the window the aperture holds 70 photons in 23379 s and 8 in 3523 s. The trigger spans 6 frames and 8 detector pixels. Frozen classification **ambiguous**.

The eight photons in six frames align well with both sky-dither coordinates and show a short peak against the full-visit aperture curve. They fall in persistent, bright Perseus gas. The CCD-wide statistic is anomalous despite modest variations over the short plotted interval. The two galaxy offsets of 5.91 and 9.33 arcseconds have chance bounds 0.23 and 0.47. The dither supports sky emission but cannot distinguish a point flare from a local fluctuation of diffuse emission.

Input `acisf29634N001_evt2.fits.gz`, 30,887,316 bytes, SHA-256 `3eb9c742e366cf37f86eccac25aaac490acc85d238cfc88d60a0e758b01c56ad`. Detector plot SHA-256 `36878d852b59e5928a394bcbbe187189546f2a05d8c8cf8c3ca5f7b5a4fb15c1`, artifact member `photometry/29634/acisf29634N001_evt2.fits.gz/real/106/diagnostics.png`. [Worker archive](https://github.com/dovvnloading/Open-Space/actions/runs/44985879443/artifacts/14745163228) · [Archive event product](https://cxc.cfa.harvard.edu/cdaftp/byobsid/4/29634/primary/acisf29634N001_evt2.fits.gz)

## Instrument passes excluded by the primary foreground veto

Each entry records the counterpart that met the frozen veto. A SIMBAD stellar class, supernova, stellar binary, remnant or stellar cluster can occur in another galaxy. A proper-motion-only match can also belong to a nearby external stellar population. The catalog's low chance-projection estimate describes alignment, not physical distance.

The indexed odd shards contain 1 real instrument-pass foreground vetoes.

- **ObsID 27145 real/0**, shard 29, 124.9231956°, +8.3694744°, 18 photons. Positive Gaia parallax meets the five-sigma criterion. Input `0913c57f6aaea85db97833b2d567e1c223d71f03ac1929153646e8cfd4c85b24`. gaia 3098103328152814976 (stellar_astrometry), 0.49″, P≤0.00613, ϖ=1.9824165305023749±0.20879191 mas, μ=(-7.9761995392939635, -2.2303810095002383) mas/yr, RUWE=0.97065.
