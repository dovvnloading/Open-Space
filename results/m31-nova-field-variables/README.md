# Two red variables near historical M31 nova candidates

## A 441-day variable near M31N 2015-07c

ZTF measured a 2.18-magnitude peak-to-peak cycle in Gaia DR3 381349291519512064, 1.56 arcseconds from the published position of M31N 2015-07c. Across 97 good ZTF i-band nights from September 2018 to September 2024, a sinusoid near **441 days** fits better than the 220-day harmonic. Dropping each calendar year in turn gives 439–456 days. The July 2015 discovery occurred 0.08 cycles after the maximum extrapolated from ZTF alone. Gaia observed this particular star before and after the discovery. [Original discovery report](https://www.astronomerstelegram.org/?read=7795) · [M31 nova catalog](https://ashafter.sdsu.edu/m31_nova_table.html)

| Observation | Measured brightness |
| --- | --- |
| Gaia G, 2015 June 24 | 19.70 mag |
| Rozhen R, 2015 July 13 | 19.73–19.82 mag |
| Gaia G, 2015 July 25 | 19.88 mag |
| Rozhen R, 2015 August 14 | 20.26 ± 0.13 mag |
| Gaia G, 2015 August 19 | 20.12 mag |

The two instruments independently measured fading within their respective passbands. The discovery team already questioned a classical nova interpretation because the August Hα excess was weak. The new measurement is a six-year optical cycle with a contemporaneously bright, red counterpart at the historical event. Gaia reports BP−RP = 2.49 mag. The closest ZTF alert centroid is 0.07 arcseconds from Gaia. The original position, the existence of a faint neighboring Pan-STARRS detection, and the absence of original calibrated discovery frames leave the association provisional. A distinct eruption adjacent to the pulsating star remains possible. The distance and stellar subtype also remain open. [Original Hα photometry](https://www.astronomerstelegram.org/?read=7921)

Fourteen independent NEOWISE W1 visit medians from 2014 to 2023 favor 436.5 days, with a nearly tied 218.5-day solution. The respective weighted residual sums are 15.61 and 15.88. The optical data distinguish these much more clearly, with scores of 197.44 at 441 days and 296.51 at 220 days. NEOWISE has a broad beam and a roughly half-year cadence. Its signal supports long-period variability in this aperture. Detections near minimum light are incomplete in both surveys. The period spread under omission of an observing year measures sensitivity to sampling and is not a formal confidence interval.

## A second field with stronger timing ambiguity

Gaia DR3 369138458974960384 lies **1.76 arcseconds** from the position of M31N 2016-06a. The International Variable Star Index lists this object as an irregular variable with no period. Gaia subsequently reported the December 2019 brightening as Gaia20afl. Its Gaia alert light curve also records G ≈ 19.13 on 2016 July 5, three weeks after the MASTER discovery, and G ≈ 19.13 on 2023 July 7. Earlier Pan-STARRS measurements show repeated large changes before MASTER first reported the candidate. [MASTER discovery and finder chart](https://www.astronomerstelegram.org/?read=9164) · [Gaia alert](https://gsaweb.ast.cam.ac.uk/alerts/alert/Gaia20afl/) · [VSX entry](https://www.aavso.org/vsx/index.php?view=detail.top&oid=8703563)

A five-passband fit across Pan-STARRS z/y, Gaia G and ZTF i/r yields a characteristic **646-day** timescale. Omitting one passband at a time gives 639.5–685 days. The independent W1 visit medians favor 674.5 days, with a competing 144-day solution. Pan-STARRS alone favors 722 days. ZTF plus Gaia alone cannot resolve 321, 214 and 636 days. The survey-dependent periods and poor joint fit preclude a precise ephemeris. The neighboring H II region, the astrometric offset and MASTER's reported June 16 non-detection after its June 14 detection further limit attribution of the historical event. The variable's repeated optical brightenings are secure, while the nature of the MASTER transient is unresolved. [Discovery and upper limit](https://www.astronomerstelegram.org/?read=9164) · [MASTER preliminary photometry](https://www.astronomerstelegram.org/?read=9166)

## Selection and replay

An exploratory comparison of 1,435 coordinate-bearing entries from the Shafter and Pietsch M31 compilation with 8,854 ALeRCE ZTF objects having at least four detections inside a 6,500-arcsecond radius around (10.684°, 41.269°) selected the 2015 field through repeated r-band alerts. An independent position comparison with the Gaia Photometric Science Alerts catalog located the 2016 field. The numbers describe the input screen. These two fields were selected for follow-up and do not estimate a contamination fraction. SIMBAD's 2015 record lists the two discovery-team telegrams. The VSX and Gaia DR3 long-period tables did not supply a period for either of the two sources in the checks on 2026 November 14. The catalog and literature checks do not establish publication priority. [ALeRCE database](https://alerce.online/) · [Gaia alerts catalog](https://gsaweb.ast.cam.ac.uk/alerts/alerts.csv) · [M31 catalog](https://ashafter.sdsu.edu/m31_nova_table.html)

Run from this directory with Python 3 and NumPy.

```sh
python analyze.py > reproduced.json
diff measurements.json reproduced.json
```

The script verifies the byte digests of each bundled archive response against [sources.json](sources.json). That manifest gives the exact original query URLs. The measurements come from ZTF at NASA/IPAC IRSA, ESA Gaia DR3, the Gaia Photometric Science Alerts team, Pan-STARRS DR2 at MAST and NASA/IPAC NEOWISE. Pan-STARRS neighbor lists and SIMBAD query returns preserve field context.

Only unflagged direct ZTF measurements with errors smaller than 0.30 mag enter the fits. The code merges same-passband observations by integer MJD, imposes an optical uncertainty floor, and fits a constant plus sine and cosine in each band on a common 130–1,500 day grid spaced by 0.5 day. Passbands have independent baselines, amplitudes and phases. Calendar-year omissions and separately fitted survey subsets expose dependence on seasonal aliases. The W1 selection, within-visit scatter floor, Gaia reject flags, and Pan-STARRS signal-to-noise and quality cuts are explicit in [analyze.py](analyze.py). The scores rank fits within each dataset. Correlated variability and survey thresholds limit their interpretation as statistical evidence.
