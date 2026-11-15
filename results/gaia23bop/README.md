# A roughly 600-day dimming timescale in NGC 1893 V22

Gaia23bop (AT 2023gdq) is the previously catalogued UX Orionis-type variable [LPC2012] V22. Three observing programs independently record a broad recurrence of optical faint states. The best periods in individual light curves span **595–619 days**. An equal-weight average of their periodograms peaks at **605.25 days**. The spread between surveys and the changing depths make approximately 600 days the appropriate precision. At the faint phase defined from Gaia alone, ZTF is typically **0.58 mag fainter** in both g and r.

## Identity and prior work

The alert position is RA 80.71265°, Dec +33.42975°. Gaia DR3 181174059990205184 lies 0.029 arcsec away (G = 15.143, parallax 0.3350 ± 0.0338 mas). SIMBAD identifies the object as Cl* NGC 1893 CUF 2116. VSX places V22 within 0.05 arcsec and lists UXOR with a blank period. The similarly named V21 lies 3.50 arcsec from the alert coordinates and has a published 9.87-day period.

Lata et al. (2012) reported its variability. Xue et al. (2019) identify V22 as star 922, report no period, and label it field/irregular. The Gaia team reported a previous dimming when it issued the 2023 alert. The addition here is the measured long recurrence timescale and the agreement of its faint phase across independently acquired optical time series. Searches on the alert name, TNS name, coordinates, NGC 1893 V22 and CUF 2116 did not turn up an earlier measurement of this long period. The apparent position in the NGC 1893 field does not settle membership.

[Lata et al. (2012)](https://doi.org/10.1111/j.1365-2966.2012.22070.x) · [Xue et al. (2019), Tables A1 and D1](https://arxiv.org/abs/1810.04923) · [Gaia alert](https://gsaweb.ast.cam.ac.uk/alerts/alert/Gaia23bop/) · [VSX positional check](https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=B/vsx/vsx&-c=80.71265%2033.42975&-c.rs=10&-out=Name,RAJ2000,DEJ2000,Period,Type)

## Independent measurements

The search spans 50–2,000 days using the same frequency grid in each band. Nightly medians suppress same-night pseudo-replication. ASAS-SN uses three-day bins after subtracting a separate median for each camera. A sinusoid evaluated at the Gaia frequency and phase scores the other surveys without changing either parameter. Positive correlations mean fainter measurements occur where Gaia predicts a faint phase. Formal period errors would omit the source's correlated variability.

| Data | Accepted detections | Time bins | Best period (d) | Power | Strongest power at 210–240 d | Correlation at Gaia phase |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Gaia G | 80 | 52 | 619.2 | 0.426 | 0.162 | 0.652 |
| ZTF g | 68 | 50 | 596.4 | 0.633 | 0.526 | 0.722 |
| ZTF r | 154 | 88 | 607.0 | 0.642 | 0.549 | 0.671 |
| ASAS-SN aperture | 578 | 426 | 595.3 | 0.218 | 0.058 | 0.398 |
| Pan-STARRS1, filter-centered | 59 | 30 | 275.5 | 0.334 | 0.071 | 0.185 |

The two ZTF filters share a survey and observing seasons. ASAS-SN is strongly blended in this cluster field. Its aperture includes nearby variable V21 and its absolute flux cannot measure V22's amplitude. Three ASAS-SN cameras give positive correlations of 0.36–0.50 separately. Gaia resolves the target. The ZTF query restricts source positions to 1.98 arcsec around V22, and its g and r centroid medians can be checked in the supplied input. Pan-STARRS1 extends coverage to MJD 55174 and records individual fades, but its preferred period differs. It contributes no independent period confirmation.

## Alternatives and limits

The 229.74-day yearly sampling alias of the Gaia period gives fixed-phase correlations of 0.010 in ZTF g, −0.171 in ZTF r and 0.086 in ASAS-SN. Both ZTF bands and the camera-centered ASAS-SN series favor the long timescale. Omitting one recurrence cycle at a time leaves positive Gaia-phase correlations in each survey. The weakest such result is r = 0.15 in ZTF r, indicating that some cycles carry much more information than others. The Gaia phase-window median magnitude contrast is 0.47 mag. The phase windows extend ±0.17 cycles about the faint and bright phases.

**The archival evidence establishes a preferred variability timescale; its predictive stability remains open.** A temporal split at MJD 58700 gives the first 24 Gaia nights a 670.5-day optimum whose extrapolation to the last 28 nights has r = −0.20. Sparse sampling and variable depths limit cycle timing. A rank-matched Ornstein–Uhlenbeck red-noise check at the actual Gaia times shows that processes with correlation times of 30, 90, 180, 360 and 600 days equal or exceed the Gaia maximum anywhere in the search band in 0.15%, 1.45%, 6.10%, 12.54% and 16.59% of 2,000 trials (plus-one estimates). These are conditional single-target simulation results. The Gaia source was chosen in an exploratory alert screen, and related measurements of one physical star are correlated. An orbital period, a disk geometry and future dated eclipses would require additional evidence.

## Reproduce

Run `python verify.py` from this folder with NumPy, SciPy and pandas. Output includes full-precision measurements, phase contrasts, camera checks, leave-one-cycle results, red-noise tests and the temporal split. The included CSVs are column subsets of the original archive responses. Date conventions retain the survey time systems, which differ by amounts negligible for a year-scale signal.

- `gaia.csv` contains the native Gaia alerts time and G magnitude. Numeric magnitudes are retained for analysis. [Gaia light curve](https://gsaweb.ast.cam.ac.uk/alerts/alert/Gaia23bop/lightcurve.csv)
- `ztf.csv` contains the IRSA light-curve service response around 80.71265°, +33.42975° with radius 0.00055°. The query returns OIDs 1701102400032249 and 1701202400039197. Use catflags = 0 and magnitude error below 0.15. [IRSA query](https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?POS=CIRCLE%2080.71265%2033.42975%200.00055&BANDNAME=g,r,i&FORMAT=CSV)
- `asassn.csv` comes from the Sky Patrol V2 master list, ID 661428833884, centered 1.37 arcsec from the alert at 80.71233°, +33.43002°. Use quality G, positive flux and flux error below 0.3 mJy. The 682 raw rows retain the original image IDs and camera names. [Sky Patrol access](https://asas-sn.github.io/skypatrol/queries.html)
- `ps1.csv` contains Pan-STARRS1 DR2 detections for objID 148110807126486302, at 80.71264185°, +33.42975356°. The check retains positive PSF flux measured above S/N 10, converts Jy to AB magnitudes and removes each filter median. [MAST detection API](https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/detection.csv?objID=148110807126486302&columns=%5BdetectID,obsTime,filterID,psfFlux,psfFluxErr,infoFlag,infoFlag2,infoFlag3,ra,dec%5D)
