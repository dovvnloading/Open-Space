# Far-infrared endpoints

The public IRAS Faint Source Catalog reproduces the first endpoint. The two AKARI single-scan fluxes and their astrometry remain author measurements. The public AKARI sky map shows no significant signal at the AKARI endpoint. A WISE-bright galaxy near the IRAS position provides a specific stationary alternative.

## Completed checks

| Source | J2000 RA and Dec (degrees) | Published flux densities | Independent check |
|---|---|---|---|
| IRAS F02211−4844 | 35.74075, −48.5125 | 60 μm 0.24 Jy, 100 μm 0.52 Jy | Native B1950 FSC row transformed to ICRS lies 0.139″ away. FSC fluxes are 0.2396 and 0.5164 Jy, qualities 3 and 2, with relative errors 16% and 23%. The 1σ astrometric ellipse is 18″ × 7″ at PA 79°. |
| AKARI-MUSL J0220440−491247 | 35.18379, −49.2135 | 65 μm 0.09 Jy, 90 μm 0.27 Jy | Table 2 transcription. Phan et al. flag the 65 μm flux as unreliable. Separation recomputed on the sphere is 47.458561′, versus 47.4586′ in the paper. |

Figure 6 shows three visits to the AKARI coordinate on June 26, 2006 at 13:35:09, 15:13:56 and 16:52:54 UT, with no convincing signal. The two reported detections are on December 26, 2006 at 20:13:46 and 21:52:29 UT. The paper gives June 26 to September 1, 1983 as the IRAS scan interval. These are observation times, whereas the coordinate system has a J2000 equinox. Phan et al. selected angular separations of 42–69.6′ under a 500–700 au distance scenario and a 7–17 Earth-mass flux model. Their two directions do not fix an orbital solution or an ensemble probability measure. [Paper, Table 2, Figures 5–6 and data statement](https://doi.org/10.1017/pasa.2025.10024)

| Product and fixed position | Forced flux (Jy) | Local robust scatter (Jy) | Contrast above control median |
|---|---:|---:|---:|
| IRIS 60 μm coadd at IRAS | 0.173 | 0.090 | 1.94 |
| IRIS 100 μm coadd at IRAS | 0.279 | 0.169 | 1.94 |
| ISSA 60 μm coadd at IRAS | 0.175 | 0.095 | 1.95 |
| ISSA 100 μm coadd at IRAS | 0.364 | 0.172 | 2.34 |
| AKARI 90 μm coadd at AKARI | −0.058 | 0.100 | −0.51 |
| AKARI 65 μm coadd at AKARI | 0.033 | 0.556 | −0.27 |

All 60/100 μm HCON products on ISSA/IRIS plate 45 and both AKARI positions are in `map_photometry.csv`. The IRAS products share observations. The stationary comparison source F02143−4941 / AKARI BSC J0216128−492720 gives 1.10 Jy and 1.31 Jy in IRIS 60 and AKARI 90 μm, with contrasts 12.0 and 12.9. Its catalogued fluxes are 0.918 and 0.885 Jy. Lower-S/N map fluxes have beam/sky systematics in addition to the measured field scatter. The public AKARI 90 μm nscan image has 9 at the endpoint. Its all-scan accumulation and the five probability panels in the paper describe different processing products. Neither endpoint has a bright-source catalogue v1 entry inside 60″.

The IRIS/ISSA measurements fit a fixed-width circular Gaussian and tilted planar background within two beam FWHM, with 4.0′ and 4.3′ beams at 60 and 100 μm. The AKARI measurements use a 90″ radius aperture, a 120–300″ median-sky annulus, and the published response factors 0.627 (65 μm) and 0.696 (90 μm). A fixed 10′-spaced control grid spans separations of 12–45′; the scatter is 1.4826 × MAD. The tabulated contrast is the flux minus the local control median divided by the scatter. The exact position was selected from the survey itself. Field contrasts are descriptive and cannot be translated into independent Gaussian false-alarm probabilities. Forced-flux values include signed noise. [IRIS beam measurements](https://doi.org/10.1086/427938) · [AKARI photometric calibration](https://doi.org/10.1093/pasj/psv033)

## Inputs and reproduction

Run `python measure.py --cache /path/to/cache` with NumPy and Astropy. `inputs.json` lists the source URL, acquisition time, byte count and SHA-256 of every input. Reproduction fails on a checksum change. The downloaded products include the arXiv paper, full CDS II/156A and II/298 catalogues, native IRIS FITS plates from IRSA, ISSA field 45 from NASA LAMBDA and AKARI all-sky v2.1 bundles from JAXA DARTS. IRIS and ISSA native coordinates use B1950 equatorial coordinates and the AKARI files use J2000 ecliptic coordinates. ICRS is used for comparing published positions. The unpublished monthly-unconfirmed-source data set is available by request according to Phan et al. The public coadd averages scans on which a mover was absent. These data cannot remeasure the per-scan signal or establish a 23-year linkage. [IRIS archive](https://doi.org/10.26131/IRSA94) · [ISSA](https://doi.org/10.26131/IRSA10) · [AKARI maps](https://darts.isas.jaxa.jp/en/datasets/darts:akari-fis-image-allsky-map-2.1) · [AKARI archive listing](https://www.ir.isas.jaxa.jp/AKARI/Archive/)

This research is based on observations with AKARI, a JAXA project with the participation of ESA.

## Possible stationary identification

NED lists the galaxy WISEA J022301.05−483040.5 at (35.75440°, −48.51127°) with photometric redshift 0.066250. Its AllWISE W4 (22 μm) measurement is 6.899 Vega mag at S/N 15.6, with cc_flags 0000. Its separation from the published IRAS position is 32.856″, PA 82.245°. Projecting onto the FSC error ellipse (18″ × 7″, PA 79°) gives an elliptical radius of 1.84, near its major axis. The optical/mid-infrared object can contribute to the far-infrared flux inside the IRAS beam. It remains a proposed counterpart until higher-resolution far-infrared observations or time-resolved astrometry can distinguish it.

The one-degree AllWISE control cone contains 22 sources with W4 S/N ≥10 and W4 magnitude ≤6.899. Its surface density is 7.00 per square degree. For a uniform independent field the expected count inside an ellipse of this size is 0.00073. The threshold and candidate were selected after inspecting the field, and galaxies cluster. This count describes local crowding. It is neither a posterior association probability nor a false-alarm rate for the original sky search. The raw queries and returned rows are in `counterpart_queries.json` and the adjacent text files. [NED redshift flags](https://ned.ipac.caltech.edu/Documents/Guides/Database)
