# Literature checks on the 16 short-cadence trigger targets

The deep-dip screen produced 142 isolated box detections on 16 white dwarfs and zero periodic detections. Several boxes cover the same epoch at different trial durations or exposure times. Two targets are established pulsators. One hosts a previously withdrawn astrometric planet claim. Individual event attribution uses the image and cadence checks separately.

All 16 positions have a unique Gaia-ID and position match in the survey TIC crossmatch. An exact Gaia DR3 identifier join to SIMBAD returned all 16 and supplied the aliases and spectral classifications below. The references identify stars with relevant published measurements. An atmospheric type alone leaves close companions and the origin of a single dip unresolved. [Survey crossmatch](tic_match.csv) · [SIMBAD identities and aliases](short_cadence/literature_identities.csv) · [Exact ADQL, response hash and retrieval time](short_cadence/literature_query.json)

| TIC | Catalogue or common name | Isolated boxes | Stellar or prior transit evidence |
|---:|---|---:|---|
| 8988748 | GD 1212 | 33 | Multiperiodic DAV. Kepler measured modes from 828.2 to 1,220.8 s. [Hermes et al.](https://doi.org/10.1088/0004-637X/789/1/85) |
| 14802607 | WDJ073607.90+144451.07 | 18 | Optical spectroscopy identifies a hot DOZ white dwarf. [Reindl et al.](https://doi.org/10.1051/0004-6361/202346865) |
| 24705587 | WD 0141−675 | 4 | The Gaia astrometric planet candidate was withdrawn. ESPRESSO and infrared data constrained the proposed companion. Its photospheric metal pollution traces accretion. [Rogers et al.](https://doi.org/10.1093/mnras/staf652) |
| 66409818 | GD 659 | 6 | DA1.4 white dwarf. |
| 199666369 | DN Dra / G 226−29 | 2 | Classical DAV with a rotationally split pulsation near 109 s. [Kepler et al.](https://arxiv.org/abs/astro-ph/0003475) |
| 219863212 | WD 1713+695 | 16 | Flux-standard study reported no detected TESS variability in its analysis. [Elms et al.](https://doi.org/10.1093/mnras/stae2265) |
| 233212451 | G 240−72 / WD 1748+708 | 12 | Strongly magnetic white dwarf identified by polarized light. [Angel et al.](https://doi.org/10.1086/181508) |
| 255908093 | Ross 640 | 17 | Metal-polluted helium-atmosphere WD. Ultraviolet spectra show Mg, Si and Fe. [Cottrell & Greenstein](https://ntrs.nasa.gov/citations/19800061028) |
| 298900716 | WD 1911+536 | 6 | An earlier TESS flux-standard study found no detectable variability. [Elms et al.](https://doi.org/10.1093/mnras/stae2265) |
| 300970893 | LAWD 26 / WD 0752−676 | 8 | DA9.0 white dwarf. |
| 323485429 | LAWD 37 / GJ 440 | 2 | DQ white dwarf. |
| 357686299 | WD 1819+580 | 4 | DA1.0 white dwarf. |
| 396487431 | GD 419 | 2 | DA2.6 white dwarf. |
| 402027326 | LB 1526 / WD 2331−475 | 2 | Hot DA white dwarf. |
| 408000326 | HZ 2 / WD 0410+117 | 6 | DA2.3 white dwarf. |
| 446165080 | GD 140 / WD 1134+300 | 4 | In the six-star CHEOPS search that found no significant transit-like event in its approximately day-long observations per target. [Morris et al.](https://doi.org/10.1051/0004-6361/202140913) |

## Pulsation frequencies outside the trigger windows

The 20-s light curves recover the characteristic stellar timescales even after masking the centers and shoulders of all isolated boxes in each product. A Lomb–Scargle scan of the undisturbed, quality-zero SAP light curves used the survey's six-hour flattening and an eight-MAD outlier cut. The ranges were 700–1,400 s for GD 1212 and 100–120 s for DN Dra. Each amplitude below is the fitted semi-amplitude of the strongest sinusoid in that range. GD 1212 has multiple variable modes. A single fitted amplitude gives no assignment for an individual 5–14% box. [Measurements, MAST URLs and FITS hashes](short_cadence/pulsation_checks.json) · [Reproduction script](short_pulsations.py)

| TIC | Sector | Fit cadences | Strongest period (s) | Fractional semi-amplitude |
|---:|---:|---:|---:|---:|
| 8988748 | 42 | 66,156 | 1,177.49 | 0.247% |
| 8988748 | 70 | 83,677 | 1,011.59 | 0.106% |
| 199666369 | 54 | 103,974 | 109.47 | 0.169% |
