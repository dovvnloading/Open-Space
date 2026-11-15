# A candidate OH megamaser toward IRAS F15496+0740

FASHI DR2 source **20260153867** has a 1405.0059 ± 0.0060 MHz line at S/N 21.04. Hydrogen at 21 cm would place it at z = 0.01096. Interpreting it as 1667.35903 MHz OH puts it at **z = 0.1867274**, just 43.8 km/s from the spectroscopic redshift of a dusty galaxy inside the FAST beam. This is a **single-line OH megamaser candidate, not a confirmed OH detection**. The independent Parkes HIPASS spectrum shows **no significant positive detection of the predicted 1665-MHz companion** at 1403.35667 MHz (or of the FAST main line at adequate significance). A **foreground H I galaxy inside the FAST beam remains possible**, including one uncataloged by DESI or SDSS. [FASHI DR2](https://nadc.china-vo.org/res/r103421/?lang=en) · [DESI DR1](https://data.desi.lbl.gov/doc/releases/dr1/)

| Evidence | Measurement |
|---|---|
| FAST position and line | RA 238.01370°, Dec +7.53255°, position uncertainty 8.28 arcsec, peak 17.28 mJy, W50 = 29.23 ± 2.52 km/s, W20 = 53.73 km/s |
| Optical spectrum | DESI 39627969517324910 has z = 0.18690074. The position is 15.47 arcsec from the radio centroid, 1.87 times its quoted uncertainty. |
| Mid-infrared host | AllWISE J155203.48+073142.0 lies 0.18 arcsec from the optical position. W1−W2 = 1.508, W2−W3 = 3.881, W4 = 6.894 ± 0.086 mag. Photometric quality AAAA, contamination flags 0000. |
| Far-infrared source | IRAS F15496+0740 lies 9.9 arcsec from the optical host. It has 0.3763 Jy at 60 µm (quality 3), 0.6960 Jy at 100 µm (quality 2), and a 6-arcsec minor positional uncertainty. |
| Conditional OH luminosity | About 750 L_sun from the 0.618 Jy km/s direct sum, or about 650 L_sun from the 0.536 Jy km/s busy-function fit. This spread indicates the integrated-flux dependence on line measurement. |
| Falsifiable companion prediction | Same-velocity 1665.40184 MHz emission would be at **1403.35667 MHz**. No separate FASHI catalog entry establishes that line. |

The luminosity uses the H I radio-velocity channel scale to convert the published velocity integral to a frequency integral, followed by 4πD_L² in a flat cosmology with H0=70 km/s/Mpc and Ωm=0.3 (D_L=908 Mpc). WISE colors and infrared brightness support an OH host without identifying the radio transition. Five SDSS DR18 spectra within 3 arcmin of FAST have no galaxy at z_HI. An independent search of primary DESI DR1 galaxies with ZWARN=0 and DELTACHI2>30 also finds no galaxy within the beam at ±300 km/s of z_HI. The nearest qualifying target in the 0.01 ≤ z < 0.02 query is 156.2 arcmin away. These checks apply to the spectroscopic sample; uncataloged foreground gas could supply the line. A spatially resolved 1667-MHz detection or detection of the predicted companion would discriminate the two hypotheses. [AllWISE](https://irsa.ipac.caltech.edu/Missions/wise.html) · [IRAS Faint Source Catalog](https://irsa.ipac.caltech.edu/IRASdocs/surveys/fsc.html)

## Selection and controls

This follow-up used the full 156,411-row FAST, 2,454,349-row DESI and 173,044-row IRAS FSC catalogs. We first screened for a radio–optical offset under 20 arcsec and under two FAST position uncertainties (6-arcsec floor), and a redshift residual under 230 km/s in the 1667-MHz interpretation. The nearest IRAS source had to lie within 60 arcsec of FAST and 30 arcsec of DESI. Four candidates passed. All are in [infrared_screen.csv](infrared_screen.csv). WISE photometry, published prior searches and foreground tests were examined **after** selection. IC 1776, one of the other matches, is a foreground source with previously measured H I near its FAST velocity. This demonstrates an actual superposition failure mode in the selected set. [IC 1776 identifications](https://simbad.cds.unistra.fr/simbad/sim-id?Ident=IC+1776)

We repeated the infrared position cuts with the velocity window centered at 20 offsets, ±600 through ±6000 km/s, in 600 km/s steps. These windows yielded eight matches in total, with at most two at one offset. The physical window yielded four. The windows and infrared selection are exploratory. These counts cannot be converted into a calibrated individual-source probability. The accepted object has one DESI galaxy spectrum. [Control counts](controls.csv)

## Direct check in a separate radio survey

The Parkes HIPASS H422 spectrum covers both predicted frequencies. We measured fixed Gaussian templates with the FASHI W50 of 29.23 km/s, fitting a local linear baseline. The corresponding fitted peaks are +8.54 mJy/beam at 1405.0059 MHz and −21.50 mJy/beam at 1403.35667 MHz. Fixed-width measurements at 291 off-line channels have a median of +0.12 mJy and robust scatter of 8.54 mJy. This gives +0.99 and −2.53 empirical noise units. The FASHI main line peaks at 17.28 mJy. There is no positive Parkes companion detection, and the negative fluctuation does not establish absorption. The data and archive attribution are in [hipass_spectrum.txt](hipass_spectrum.txt). [ATNF extraction](https://www.atnf.csiro.au/cgi-bin/multi/release/plotspectrum.cgi?ra=15%3A52%3A03.29&dec=%2B07%3A31%3A57.2&datasource=hipass&xaxis=freq&lo=1395&hi=1419)

Run `python ../check_parkes.py --spectrum hipass_spectrum.txt --centers 1405.0059 1403.35667 --width-kms 29.23` to reproduce [both measurements](parkes_measurements.json).

## Previous OH observations and limits

We searched the full detection and upper-limit tables in Darling & Giovanelli's first two Arecibo reports and the third report's machine-readable tables of detections, non-detections, ambiguous cases and previously known megamasers. IRAS 15496+0740 does not appear among these targets. We also checked the 2006 Darling et al. optical follow-up list, Zhang et al.'s FASHI OH study and Roberts et al.'s ALFALFA follow-up using the IRAS identifier, nearby coordinates and WISE/FASHI aliases. We found no published OH detection or **source-specific prior upper limit** in those searches. The IRAS source is below the 0.6 Jy nominal PSCz 60-µm selection limit; a 5-arcmin cone around the WISE host in the PSCz main catalog and gap supplement is empty. The Arecibo and targeted FASHI searches drew heavily on PSCz. A fresh sensitivity bound would require the spectrum and noise at the relevant frequencies. [Darling & Giovanelli I](https://arxiv.org/abs/astro-ph/0005046) · [II](https://arxiv.org/abs/astro-ph/0102345) · [III tables](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/AJ/124/100) · [Zhang et al.](https://arxiv.org/abs/2401.15397) · [Roberts et al.](https://arxiv.org/abs/2506.06115)

The 2026 FASHI paper gives an aggregate total of 18 OH-consistent cases without individual IDs in its discussion. The actual DR2 table contains 20260153867. Overlap with that unpublished candidate list is unresolved. [FASHI DR2 §6.6](https://arxiv.org/html/2606.31539v1)

## Replay

Run the J044309−064913 candidate analysis first to cache its FASHI DR2 and 11 DESI DR1 redshift slices. Use the same cache with the following command. This script adds the checksum-verified HEASARC copy of IRAS FSC, reselects the four catalog matches and recalculates the 20 controls. The evidence files include the primary FASHI and DESI rows, the native AllWISE and IRAS catalog responses, the PSCz cone, and the SDSS DR18 3-arcmin spectral result. [HEASARC IRAS FSC](https://heasarc.gsfc.nasa.gov/W3Browse/iras/irasfsc.html) · [NOIRLab DESI DR1](https://datalab.noirlab.edu/data/desi)

```sh
python replay.py --cache /tmp/fashi-oh-sources --output /tmp/fashi-oh-infrared
```

This work has used the data from the Five-hundred-meter Aperture Spherical radio Telescope (FAST). FAST is a Chinese national mega-science facility, operated by the National Astronomical Observatories of Chinese Academy of Sciences (NAOC).
