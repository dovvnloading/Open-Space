# A candidate OH megamaser at WISEA J044309.47−064913.0

FASHI DR2 source **20260037137** has a 1407.2941 ± 0.0162 MHz emission line at catalog S/N 13.66. Two DESI spectra within 2.8 arcsec of its centroid give galaxy redshifts near 0.1845. If the radio line is the 1667.35903 MHz hydroxyl transition, its redshift is **0.18479785**, near both optical velocities, and its isotropic luminosity is approximately **470 solar luminosities**. The published 21-cm H I interpretation instead puts the emitter at z = 0.00932, with a catalog distance of 36.1 Mpc. These data establish a **single-line OH megamaser candidate, not a confirmed OH detection**. The independent Parkes HIPASS spectrum shows **no significant positive detection of the predicted 1665-MHz companion** at 1405.64218 MHz (or of the FAST main line at adequate significance). A **foreground H I galaxy unresolved within the 174-arcsec FAST beam remains possible**. A second OH line or interferometric localization is needed for a secure physical identification. [FASHI DR2](https://nadc.china-vo.org/res/r103421/?lang=en) · [DESI DR1](https://data.desi.lbl.gov/doc/releases/dr1/)

| Measurement | Value |
|---|---|
| FASHI position (J2000) and nominal centroid uncertainty | 70.79013°, −6.82007°, 12.74 arcsec |
| Line peak, W50, and W20 | 6.82 mJy, 44.91 ± 6.84 km/s, 120.74 km/s |
| Integrated flux | 0.39328 Jy km/s by channel sum, 0.38761 Jy km/s by busy-function fit |
| DESI 39627623118146959 | z = 0.18425289, 2.14 arcsec from FAST, optical minus OH velocity −137.93 km/s |
| DESI 39627623118146957 | z = 0.1845830 (rounded), 2.71 arcsec from FAST, optical minus OH velocity −54.37 km/s |
| AllWISE counterpart | J044309.47−064913.0, 2.51 arcsec from FAST, W1−W2 = 0.833 mag, W2−W3 = 3.941 mag, W4 = 6.444 ± 0.062 mag, four A-quality bands, contamination flags 0000 |
| NVSS counterpart | J044309−064915, about 5.6 arcsec from FAST, 1.4 GHz continuum flux density 5.1 ± 0.5 mJy |
| Prediction | A 1665.40184 MHz OH component at the main-line redshift would be centered at **1405.64218 MHz**. |

The DESI targets are only 1.18 arcsec apart and their redshifts differ by about 84 km/s. The spectra may sample overlapping light, gas motions, or distinct components. They do not independently demonstrate a merger. The mid-infrared colors and nearby continuum emission strengthen the case for a dusty radio host. Neither measurement identifies the line species. The 174-arcsec FAST beam leaves room for a faint superposed foreground H I galaxy. [DESI spectra and identifiers](https://data.desi.lbl.gov/doc/access/) · [AllWISE](https://irsa.ipac.caltech.edu/Missions/wise.html) · [NVSS](https://www.cv.nrao.edu/nvss/)

A separate search using IRAS infrared sources identified FASHI 20260153867 toward IRAS F15496+0740. Its 1405.0059 MHz line would place it just 44 km/s from its optical host under the OH interpretation. [Second candidate, evidence and controls](iras-f15496/README.md)

## Search, foreground checks, and controls

The discovery screen compared all 156,411 DR2 radio positions with 2,454,349 DESI DR1 primary galaxy redshifts from 0.17 through less than 0.28, requiring ZWARN = 0 and DELTACHI2 > 40. An initial 90-arcsec, ±350 km/s crossmatch returned 851 optical associations involving 780 radio rows. Ranking by positional residual exposed the two DESI targets here. Within 3 arcsec and ±150 km/s they are the only two optical entries passing in the whole sample and both correspond to this radio source. Velocity offsets use c ln[(1+z_opt)/(1+z_OH)].

Keeping positions fixed, 40 comparison windows shifted the required optical velocity by ±300 through ±6000 km/s in 300 km/s steps. Thirteen windows returned at least one association, with 16 radio/optical associations in total. None returned two optical targets at one radio source. These narrower cuts were examined after seeing the candidate. The trials and the two spectra have correlations, so the controls are descriptive and cannot be used as a discovery significance. [Output by shift](controls.csv)

The DESI field extraction includes all primary ZWARN = 0 targets regardless of redshift or spectral type. Of 32 in the rectangular extraction, 23 are within 174 arcsec and six within the 87-arcsec half-power radius. None is within 300 km/s of the 21-cm velocity. There is no second cataloged FASHI line at the predicted companion frequency within 9 arcmin. The public tables cannot supply a channel-by-channel upper limit there. FASHI reports W50 on a velocity axis referenced to H I. To estimate the conditional luminosity, we convert the catalog channel-summed flux on that radio axis into integrated frequency flux using 1420.40575 MHz / c. A flat cosmology with H0 = 70 km/s/Mpc, Ωm = 0.3 gives D_L = 897.5 Mpc at z_OH. Multiplying the integrated flux by 4πD_L² gives 469 L_sun. The 1665-MHz line may be weaker than the original detection. Spatially resolving the 1667-MHz signal against the redshifted optical host would also discriminate foreground H I. [Catalog definitions and data availability](https://arxiv.org/html/2606.31539v1)

## Independent Parkes spectrum

The HIPASS H295 cube covers this position and both frequencies. At J2000 04:43:09.63 −06:49:12.3, the 384-channel Parkes spectrum spans 1395.00–1418.94 MHz at 62.5 kHz spacing. We fit a Gaussian of fixed 0.211 MHz FWHM and a local linear baseline within ±1 MHz of each line. The same fit at 292 off-line channel centers gives a robust peak-amplitude scatter of 7.76 mJy/beam. Neighboring channel fits are correlated. After subtracting the off-line median, the signal-to-noise ratios are +0.29 at 1407.2941 MHz and −0.89 at 1405.6422 MHz. The FASHI peak is only 6.82 mJy. Even the detected main line is below one Parkes noise unit in this test, so these measurements leave both interpretations open.

Run `python check_parkes.py` to replay both forced fits on [hipass_spectrum.txt](hipass_spectrum.txt). The archive attribution is retained with the spectrum. [ATNF spectrum extraction](https://www.atnf.csiro.au/cgi-bin/multi/release/plotspectrum.cgi?ra=04%3A43%3A09.63&dec=-06%3A49%3A12.3&datasource=hipass&xaxis=freq&lo=1395&hi=1419) · [Barnes et al., HIPASS observations and reduction](https://doi.org/10.1046/j.1365-8711.2001.04102.x)

## Previous work and priority

Zhang et al. reported 27 OH emitters in an IRAS PSCz-targeted FASHI search. Roberts et al. used optical spectroscopy to identify additional OH emitters among apparent H I detections in ALFALFA. The 2026 FASHI DR2 release paper reports 18 OH-consistent cases in its optical cross-check and describes their exclusion when matching redshifts were available. It gives aggregate counts in §6.6. The official DR2 extragalactic H I table still contains this specific source. We found no itemized prior OH identification of this source in the checked targeted FASHI survey, Roberts et al., or searches for the FASHI, AllWISE and NVSS names. Whether it appeared in the FASHI team's internal list remains unresolved. [Zhang et al. 2024](https://arxiv.org/abs/2401.15397) · [Roberts et al. 2025](https://arxiv.org/abs/2506.06115) · [Zhang et al. 2026](https://doi.org/10.1007/s11433-026-3072-1)

## Reproduction and source products

With Python, NumPy, pandas and SciPy installed, run the following command. The full recalculation downloads the FAST catalog and eleven redshift slices of the NOIRLab DESI DR1 zpix table. It takes approximately 275 MB of source data, checks the FAST SHA-256, checks for capped DESI slices, verifies the optical sample size and unique IDs, and writes fashi.csv, optical_matches.csv and controls.csv.

```sh
python replay.py --cache /tmp/fashi-oh-sources --output /tmp/fashi-oh-replay
```

The independent field evidence is included in desi_field.csv, allwise.tbl and nvss.tsv. The DESI field query used the same TAP endpoint as replay.py and selected targetid, mean_fiber_ra, mean_fiber_dec, z, zerr, spectype, zwarn and deltachi2 from desi_dr1.zpix. It used the radio position ±0.05° in Dec and ±0.05°/cos(Dec) in RA with zwarn=0 and zcat_primary=true. AllWISE returned its native IPAC table in a 5-arcsec cone. NVSS returned its VizieR catalog VIII/65/nvss in a 0.25-arcmin cone. All three use the radio centroid above. [DESI Data Lab](https://datalab.noirlab.edu/data/desi) · [IRSA AllWISE catalog](https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-dd?catalog=allwise_p3as_psd) · [NVSS VizieR catalog](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=VIII/65)

This work has used the data from the Five-hundred-meter Aperture Spherical radio Telescope (FAST). FAST is a Chinese national mega-science facility, operated by the National Astronomical Observatories of Chinese Academy of Sciences (NAOC).
