# Galactic nebulae misidentified as extragalactic H I in FASHI

A frequency-and-position audit finds **19 published FASHI extragalactic H I entries with plausible Galactic hydrogen recombination-line identifications** toward four known H II regions. Six rows are from DR1 and 13 from DR2. The strongest case is Sh 2-235. Four consecutive hydrogen recombination lines in DR2 give the same Galactic velocity within 1.7 km/s. As 21-cm H I emitters, the catalog assigns them four distances from 56.2 to 298.1 Mpc and four logarithmic H I masses from 9.108 to 10.351 in solar units. The entire quartet lies within 88.2 arcsec, about half the telescope's 174-arcsec beam FWHM.

These are **proposed source-level reidentifications** of emission from known nebulae, strongest for co-spatial multi-line systems. The FASHI team had **already described and illustrated this Galactic recombination-line contamination mechanism** in DR1 and acknowledged possible residual contamination in DR2. The new contribution is identification of the specific retained catalog rows, not discovery of the contamination mechanism or the nebulae. The two releases use overlapping observations. An entry is a row in one release and does not necessarily represent a spatially or statistically independent observation.

| Galactic H II region | DR1 catalog rows | DR2 catalog rows | Independent ionized-gas velocity (km/s, LSR) | Evidence |
|---|---:|---:|---:|---|
| Sh 2-235 | 2 | 6 | −23.07 | Four central DR2 α transitions, two further spatial components at α frequencies, two DR1 β transitions. The central lines have S/N 23.1–27.1 and widths 19.4–21.6 km/s. The external spectrum has width 20.52 km/s. |
| Sh 2-206 | 1 | 4 | −25.20 | Two DR2 α lines 30.0 arcsec apart, two DR2 β entries and one DR1 β entry. The α lines have S/N 34.2 and 48.7, with widths 23.9 and 25.2 km/s. |
| Sh 2-156 / IC 1470 | 3 | 1 | −49.55 | Three consecutive DR1 α transitions within 26 arcsec. DR2 retains H168α at almost exactly its original frequency and coordinates. Its broad line and −11.1 km/s velocity residual make this component individually less diagnostic. |
| Sh 2-297 | 0 | 2 | +1.94 | H169α and H170α are 6.0 arcsec apart. Their velocities agree within 1.9 km/s and lie 6.0 and 7.9 km/s above the external centroid, inside the external 22.6 km/s line width. |

## Search and checks

The observing unit is one FASHI catalog extraction. We used the published line centers and W50 linewidths, with no candidate list supplied to the search. Of 41,741 DR1 and 156,411 DR2 rows, 12,345 and 44,234 respectively have a positive width below 100 km/s. The independent comparison sample comprises 2,510 WISE Galactic H II regions with a published velocity. The WISE reference measurements for the four regions used higher-frequency transitions observed with other instruments. The source measurements are identified in [regions.csv](regions.csv). [WISE H II catalog](https://astro.phys.wvu.edu/wise/) · [Quireza et al. 2006](https://arxiv.org/abs/astro-ph/0603133) · [Lockman 1989](https://adsabs.harvard.edu/pdf/1989ApJS...71..469L) · [Balser et al. 2011](https://doi.org/10.1088/0004-637X/738/1/27)

The seed search considered H167α through H170α throughout each release. A region needed two different transitions, each within the published region radius and 180 arcsec of its center. The pair had to be within 180 arcsec of one another and agree within 12 km/s. Each line had to lie within 16 km/s of the external LSR velocity. We derive heliocentric radio velocities from the observed and hydrogen reduced-mass Rydberg frequencies and add the projected standard 20 km/s kinematic LSR solar motion. The two velocity gates allow for the 6.4 km/s FASHI spectral resolution and resolved nebular gas. These are exploratory selection choices.

Only after selection did the search use the full published angular extent of each of the four regions to look for additional α components and β lines in both releases. It examined H208β–H215β with the same width and external-velocity cuts. Five β matches survived. Outer components can have different spatial centroids in extended Galactic emission. The physical interpretation is strongest for the co-spatial multi-line seeds. Every accepted entry and the catalog distance and mass that would follow from the H I interpretation are in [identifications.csv](identifications.csv).

Forty comparison runs moved all WISE centers together along Galactic longitude by ±0.5° to ±10° in 0.5° steps while preserving latitude, radius and velocity. **Zero** comparison runs produced any seed region in either release. The true positions gave one DR1 and three DR2 seed regions. These controls test sensitivity to positional alignment. They share catalogs and correlated Galactic structure, and do not establish a calibrated false-discovery probability. Other isolated line coincidences outside the selected regions remain unclassified. [Control results](audit.json)

## Prior work and limits

The DR1 authors specifically describe removing Galactic recombination lines and illustrate an H167α example at J054023.71+355217.4. This is in the **same Sh 2-235 complex**, about seven arcminutes from the central quartet. DR2 explicitly acknowledges possible residual recombination contamination. The contribution here is the table of surviving released catalog IDs, the line-by-line assignments, validation against independent nebular velocities, and whole-catalog controls. A search on November 15, 2026 of the release papers, these regions with FASHI, and selected exact catalog IDs found no prior itemized list of these particular rows. This is a bounded literature check. [DR1 §4.5 and Fig. 12](https://zcp521.github.io/pub/fashi.pdf) · [DR2 §6.6](https://arxiv.org/html/2606.31539v1)

The catalog records cannot rule out a weaker superposed 21-cm galaxy blended with any individual Galactic line. The consistent velocities across independent atomic transitions and their localization on the nebulae favor a dominant foreground contribution. Source cubelets and recomputed line profiles would test blends and the spatial morphology directly. The survey paper makes the reduced cubelets available by author request. The row counts here do not estimate the survey-wide contamination fraction or change to its H I mass function.

## Reproduction

With Python, NumPy, pandas and SciPy installed, run the following command from this directory.

```sh
python replay.py --cache /tmp/fashi-inputs --output /tmp/fashi-results
```

The script retrieves and SHA-256 verifies the three official source tables linked in its FILES map. It reuses verified cached files on subsequent runs. Comparison of the generated CSV and JSON products with the copies here reproduces both the source IDs and negative controls. The distances are in Mpc, mass is log10(M_HI/M_sun), coordinates are J2000 degrees, frequencies and frequency errors are in MHz, and W50 and velocities are in km/s. [FASHI DR1 data](https://nadc.china-vo.org/res/r100875/?lang=en) · [FASHI DR2 data](https://nadc.china-vo.org/res/r103421/?lang=en) · [WISE Galactic H II region catalog v3.0](https://astro.phys.wvu.edu/wise/)

This work has used the data from the Five-hundred-meter Aperture Spherical radio Telescope (FAST). FAST is a Chinese national mega-science facility, operated by the National Astronomical Observatories of Chinese Academy of Sciences (NAOC).
