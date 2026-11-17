# JWST repeat-imaging transient search

## November 2026 continuation

Two previously unmatched **nuclear brightenings** pass individual-exposure checks. The existing catalog identifications refer to their host galaxies. The cause of each flux increase remains unclassified.

| Position (J2000, degrees) | Change in aperture flux | Host distance and interpretation |
|---|---|---|
| **COSMOS-022**, 150.122033, +2.323386 | F162M +236.2 ± 18.2 nJy, F182M +212.1 ± 18.1 nJy. Both late-May observing dates are elevated relative to December–January. | COSMOS2020 953994, MUSE z = 0.5067. The source is coincident with the galaxy center. The matched rest-band F162M absolute difference magnitude is −16.45. [Images, exposures and host evidence](repeat-catalog-followup/README.md) |
| **GOODS-S 185**, 53.188672, −27.812554 | F277W +22.25 ± 1.15 nJy from October 2025 to October 2026. The earlier 2022–25 measurements were near 49–52 nJy. Five filters brighten. | JADES DR2 103607, host photo-z 2.47 (68% 1.91–2.51). The matched rest-band F277W absolute difference magnitude is −17.18. [Flux history, stamps and host evidence](goods-s-oct2026/README.md) |

The absolute magnitudes refer to the measured change in a fixed aperture under the host association. Neither object has a transient spectrum. AGN variability and a transient projected onto the host nucleus remain possible. TNS mirrors, published transient/AGN lists and the available program target lists constrain novelty. The direct TNS guest snapshot was inaccessible.

| Selection | Measured outcome |
|---|---|
| NEXUS spectroscopic z > 5 galaxies | No validated joint event in 1,577 temporal contrasts at 100 host positions across eight visits. Positive one-visit, flat-fν injections gave joint 90% and 50% recovery at AB 26.59 and 27.13. A separate September catalog check covers fewer hosts. [Search and rejected z = 5.612 deficit](high-z-hosts/README.md) |
| Blind GOODS-S difference images | About 2.33 arcmin² per filter. Of 31 two-filter coincidences, 28 followed flagged detector pixels. The remaining three have known-source or resolved-residual explanations. At AB 27.5 single-filter recovery was 72.1% (F277W) and 75.9% (F356W). [Blind search and injections](goods-s-oct2026/README.md) |
| MACS0416 Warhol and Spock arcs | Zero same-sign two-band peaks in 0.190136 arcmin². Equal-flux two-band injections gave 50% recovery at AB 26.65. A repeatable F444W excess lies 0.072″ from published D31-W5. [Split exposures and critical-curve distance](macs0416/README.md) |
| Cataloged sources in repeated fields | 69 catalogs, 37 matched filter pairs and 47,584 eligible source/filter comparisons in COSMOS, CEERS and three cluster fields. The independent peak and exposure checks include known quasar L1 of SMM J02399−0136. [All 42 catalog trigger groups and dispositions](repeat-catalog-followup/README.md) |
| Abell 2744 and MACS0647 | No factor-of-two, ≥7σ change in 700 compact-source/filter entries detected in both epochs. [Dates, filters and depth diagnostics](cluster-catalog-variability/README.md) |
| El Gordo | Earlier branch result: zero joint F150W/F444W detections in 0.879293 arcmin² in the July 2022 and August 2026 cluster-core comparison. [Measured area](el-gordo/screen.json) · [Injection bins](el-gordo/injections.csv) |

These are specific footprints, cadences and selection functions. Direct coverage at every epoch across the named surveys has not been established.

## Earlier NEXUS pairwise imaging

No new extragalactic transient has been established at this checkpoint. The search has acquired actual single-epoch NEXUS science and error maps. The first full-footprint subtraction compares NIRCam F444W from 2026 May 28 and August 7. A coverage-selected central square also has F200W in both epochs. An August–September F444W comparison has also been screened; its detector-level and source-history checks remain incomplete.

[qdr-inputs.csv](qdr-inputs.csv) records the URLs and SHA-256 hashes of the first eight compressed archive inputs and the derived central cutouts. [screen.py](screen.py) fits a relative astrometric translation and empirical compact-source PSFs, cross-convolves the epochs, searches both signs of the difference and measures conditional point-source recovery by injection on the observed background. The saved outputs include empirical calibration, usable area and conditional image-injection measurements. Candidate vetting, independent-exposure checks and survey-wide false-alarm calibration remain incomplete.

Data credit: the NEXUS collaboration, [QDR archive and release documentation](https://nexus-jwst.github.io/updates/). [Survey overview (Shen et al.)](https://arxiv.org/abs/2408.12713) · [Release notes (Zhuang et al.)](https://arxiv.org/abs/2603.04586) · [Published year-one transients](https://arxiv.org/abs/2609.06985).

The blind May–August run measured 18.719 arcmin² in F444W and 3.791 arcmin² in F200W after coverage and boundary masks. Its two-sided seven-sigma threshold and point-source shape screen retain 25 and 12 peaks respectively. These are detector outputs. The strongest F200W brightening coincides with the [previously reported AT 2026abcf](https://www.wis-tns.org/object/2026abcf/). [F444W measurements](may-aug-f444w/) · [F200W measurements](may-aug-f200w/). Host, exposure-level and additional-epoch checks are in progress.

The [four-field matched-catalog screen](repeat-catalogs/README.md) contains multi-filter rankings, not confirmed transients; its selection requires a cataloged positive-flux source in both epochs. A separate [MACS J0416 analysis](macs0416/README.md) reports a Warhol-arc residual near previously published, unresolved microlensed stars and a bounded two-band image search. It does not establish a new source.
