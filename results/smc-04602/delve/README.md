# DECam resolves both eclipses of OGLE-SMC-ECL-04602

Archival DECam photometry in the DELVE-MC measurement catalog puts both eclipses at the OGLE binary's position. On October 29, 2020, 49 exposures followed the secondary eclipse through ingress and egress. A simultaneous fit to four filters measured its center at **HJD 2459151.5696** (UTC time scale), orbital phase **0.6121**. A brighter Gaia source lies 2.11 arcseconds from the OGLE position, making the spatial check important. The DECam photometry resolves the two sources.

## Location and timing

At J2000 01:02:55.78, −72:52:23.8, the closest DELVE target in the selected catalog tile is 0154m727.1965, 0.11 arcsec from OGLE. Neighbor 0154m727.1896 coincides with the brighter Gaia source. DELVE also lists detections of the same stars on an overlapping tile. The analysis uses only one tile. It retains ordinary point-source measurements with positive errors ≤0.10 mag, selects the lowest quoted error for a repeated measurement of the same star and exposure, and pairs the target and neighbor by exposure identifier. This gives 248 matched g, r, i and z exposures on 32 nights.

The OGLE catalog epoch is HJD 2460001.4679 and period is 1.7295874 days. We estimated eclipse phases in three date windows from OGLE-IV I alone using two wrapped Gaussian dips, a constant, and a second-harmonic term. The eclipse core was within 0.03 orbital cycle of a predicted center. Baseline exposures lay more than 0.12 cycle from both centers. Each filter's baseline is the median of its nightly medians. There were 12, 20, 6 and 3 baseline nights in g, r, i and z. Positive differences denote dimming.

| Eclipse | HJD observing night (floor) | Matched core exposures | Target dimming (mag) | Neighbor change (mag) |
| --- | ---: | ---: | ---: | ---: |
| Primary | 2458484 | 2 | 0.527 | +0.006 |
| Primary | 2458674 | 1 | 0.486 | +0.009 |
| Primary | 2459150 | 12 | 0.454 | −0.061 |
| Secondary | 2458447 | 2 | 0.231 | +0.002 |
| Secondary | 2459151 | 30 | 0.245 | −0.033 |

For the last secondary night we fitted all 49 paired measurements. A Gaussian with a common center and width and free amplitudes and baselines in g, r, i and z gave depths of 0.341, 0.352, 0.327 and 0.357 mag. Ordinary and differential photometry, three eclipse shapes and five channel selections gave 30 center estimates spanning phase 0.6110–0.6140. These fits share exposures and the range measures model sensitivity. The same date range of OGLE data independently places the secondary near phase 0.611. The previous primary was observed near minimum and during egress. The DECam sequence provides an independent timestamp for the secondary and spatial evidence for both eclipses. OGLE supplies the long time baseline for any historical trend. Sources below the DECam resolution could still contribute.

## Prior reports and catalogs checked on November 14, 2026

The binary has been cataloged since OGLE-II. Aliases include SMC_SC9_141913, SMC110.7.5219 and SMC726.21.25362. It is absent from both the OGLE-III and OGLE-IV SMC remarks lists. Searches for 4602 and 141913 found no corresponding entry in the 27 targets of Hong et al. (2015), 90 targets of Hong et al. (2016), 18 targets of Zasche et al. (2014), or 21 targets of Zasche & Wolf (2019). SIMBAD links the two OGLE catalogs and the 2021 morphology paper. None of those SIMBAD titles describes an eclipse-timing analysis dedicated to this star. The VizieR VSX cone search returned no rows within 6 arcsec. These searches cannot establish publication priority. Gaia DR3 has no released epoch-photometry series among the four sources in that radius. The closest Gaia source is 4687420963524480640 (0.13 arcsec, G = 18.345) and the brighter neighbor is 4687420967788454272 (G = 18.067).

## Replay from the archived measurements

This directory includes the raw SQL/TAP query responses used for the DECam measurement, a frozen OGLE-IV I light curve, and machine-readable results. From this directory, with Python 3.12, NumPy and SciPy installed, run the following.

```sh
sha256sum -c SHA256SUMS
python delve04602.py data data/iv-04602-I.dat
python secondary04602.py data
```

A rerun of `delve04602.py` can differ from the committed `delve-analysis.json` and `delve-analyzed.csv` in late digits. The script uses SciPy's least-squares solver with its default stopping tolerances, so the stopping point depends on rounding in the installed NumPy and SciPy. On 16 November 2026 a rerun on Windows 11 with Python 3.14.0, NumPy 2.5.4 and SciPy 1.18.1 changed fitted values by at most 8.6 parts in a million. No number quoted in this folder changes, and the committed files are unchanged.

The first script converts UTC exposure start MJD to an approximate heliocentric Julian date at exposure midpoint, using exposure lengths from the DELVE exposure table. It writes matched exposures to data/delve-analyzed.csv and nightly measurements to data/delve-analysis.json. The second writes the event fits to data/delve-secondary.json. sources.json gives the GET URL, byte count and SHA-256 for every archived input, including the exact TAP queries. A data service may change row order or response headers on a fresh query. The archived bytes permit offline replay.

Original observations and catalog credit belong to OGLE, DELVE and the original DECam observing programs, Gaia/DPAC, CDS SIMBAD and VizieR. [NOIRLab DELVE data and citations](https://datalab.noirlab.edu/data/delve) · [DECam observations catalog](https://datalab.noirlab.edu/tap) · [OGLE-IV measurements](https://www.astrouw.edu.pl/ogle/ogle4/OCVS/smc/ecl/phot/I/OGLE-SMC-ECL-04602.dat) · [OGLE-II discovery catalog](https://acta.astrouw.edu.pl/Vol48/n4/pap_48_4_1_7.pdf) · [SIMBAD source](https://simbad.cfa.harvard.edu/simbad/sim-id?Ident=OGLE+SMC110.7+5219) · [Hong 2015 catalog](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/AJ/150/1) · [Hong 2016](https://arxiv.org/abs/1604.02200) · [Zasche 2014](https://arxiv.org/abs/1412.0905) · [Zasche & Wolf 2019](https://arxiv.org/abs/1902.04961)
