# CD−27 2812 and a candidate second binary

CD−27 2812 (TIC 37606218) is a benchmark eclipsing F9 V + M-dwarf binary. Its known wide Gaia companion has separately been classified as an unresolved white dwarf plus main-sequence star. Joining these observations gives a **candidate 2+2 quadruple** with a white dwarf in one branch. The proposed four-star architecture and the orbital motion of the faint branch await confirmation. The two resolved sources and their wide association were already cataloged. A search through the references below and object aliases did not locate a publication explicitly identifying the combined 2+2 architecture. Publication priority remains open.

## Measurements

| Observable | Result |
|---|---:|
| Bright Gaia DR3 source, eclipsing binary | 2897976170198264064 |
| Faint Gaia DR3 source, WD+MS candidate | 2897976062821487360 |
| Angular separation at Gaia epoch 2016.0 | 53.7153 arcsec |
| Bright and faint parallaxes | 5.1457 ± 0.0146 and 5.0831 ± 0.0738 mas |
| Projected separation at the brighter source's distance | 10,439 au |
| Relative sky-plane motion and speed at that distance | 0.284 mas/yr and 0.262 km/s |
| Gaia spectroscopic period of the bright binary | 7.83571 ± 0.00024 d |
| Faint source Gaia BP−RP and absolute RP at bright source distance | 1.726 and 10.169 mag |

The brighter system has two measured dynamical masses, 1.3597 and 0.5624 solar masses. The observed sky-plane speed falls below 0.572 km/s, the escape-speed scale from those masses alone at the projected separation. Both the line-of-sight separation and relative radial velocity are unknown. The bright binary's center-of-mass velocity in Gaia is 4.79 km/s. High-resolution spectroscopy and TESS photometry independently establish its components. [Adshead et al. 2026](https://arxiv.org/abs/2604.13631)

El-Badry, Rix & Heintz listed this exact pair in 2021 with chance-alignment density ratio 0.00198. In our independent positional selection on the same mission's DR3 astrometry, only these two objects, among 3,365 sources with parallax signal-to-noise above five within 0.5° of the bright system, fall within 0.5 mas in parallax and 1 mas/yr in vector proper motion of CD−27 2812. Doubling both tolerances selects the same two. These diagnostics characterize the association. Neither count is a false-alarm probability. [Wide-pair catalog](https://doi.org/10.1093/mnras/stab323)

Rebassa-Mansergas et al. included the faint source in their catalog of 1,312 unresolved WD+MS systems. Its row has no reported temperatures, masses or orbital period, and its companion-contamination flag is zero. Gaia reports zero blended BP and RP transits among 48 and 46 observations respectively. The separate source 2MASS J06130097−2753401 has J = 14.903 ± 0.038 and Kₛ = 14.071 ± 0.054, with AAA quality and 000 confusion flags. At similar absolute RP magnitude, 867 of 868 nearby Gaia comparison stars have redder BP−RP colors. The comparison uses the published Gaia photometry again and was restricted to high-quality, nearby sources. Its color rank helps show why to obtain a spectrum. Selecting a source from a catalog partly based on color precludes interpreting this rank as a discovery significance. [WD+MS catalog](https://doi.org/10.1051/0004-6361/202554700) · [2MASS Point Source Catalog](https://doi.org/10.26093/cds/vizier.2246)

## Selection and prior coverage

The search began with every source in the published WD+MS sample. The frozen cones produced 81 spatial and kinematic matches to its 1,312 targets. Only the CD−27 2812 match has a resolved neighbor with a nonzero Gaia DR3 non-single-star flag. This exploratory selection preferentially finds well-measured, nearby systems and allows a finite orbital velocity. The count is a screening yield. [`analyze.py`](analyze.py) implements the actual cuts. [`data/sources.md`](data/sources.md) gives the remote queries, snapshot provenance and control selection.

On 2026-11-14, 90-arcsec coordinate searches around the bright system returned no entries in the VizieR December 2023 release of Tokovinin's Multiple Star Catalog or the Washington Double Star Catalog query. SIMBAD linked the bright object to its eclipsing-binary and spectroscopic studies. Searches on CD-27 2812, TIC 37606218, and the two Gaia identifiers located the established bright-binary literature. Catalog omissions cannot establish publication priority. A projected coincidence remains possible, as does a different explanation for the faint source's composite colors. [Multiple Star Catalog](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/ApJS/235/6) · [Washington Double Star Catalog](https://crf.usno.navy.mil/wds)

## Falsification and replay

A spatially isolated optical spectrum of Gaia DR3 2897976062821487360 can test for a hot white-dwarf continuum or broad Balmer lines together with cool-star molecular bands. A two-component fit should outperform a single-star fit after accounting for flux calibration, metallicity and extinction. Multi-epoch spectroscopy of the cool star can test for an orbit and recover a systemic radial velocity for comparison to the benchmark binary. A single-star spectrum that fits both optical and near-infrared data, or discordant systemic velocities after accounting for orbital motion, would weaken the proposed architecture. An established independent companion at this distance would tie the white-dwarf evolution to a binary with precisely known masses.

Run `python analyze.py` in this directory with Python 3.10 or newer. The standard-library replay reads the frozen search and control data and prints the contents of [`result.json`](result.json). `sha256sum -c SHA256SUMS` verifies the files.
