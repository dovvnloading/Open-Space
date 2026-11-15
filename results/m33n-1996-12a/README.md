# A 544-day Mira at the site of M33N 1996-12a

Williams & Shafter measured a transient in M33 on 1996 December 9 in their 70 Å Hα filter. Their six-object survey gives one magnitude for this candidate, 19.1 after extinction correction, with no tabulated B magnitude or spectrum. The published position is **0.394 arcsec** from the oxygen-rich fundamental-mode Mira candidate **01342067+3034391**. The authors quote approximately 1 arcsec positional accuracy. The later CFHT catalog gives a period of **543.85 days**, assigns all six machine classifiers to the positive class and reports a high-confidence visual check. This association supplies a specific alternative to a thermonuclear eruption. [Original survey, Tables 2–4](https://doi.org/10.1086/422833) · [CFHT catalog, Table A7](https://doi.org/10.1093/mnras/stae1140)

## Independent observations

| Measurement | Observed property |
|---|---|
| CFHT fitted magnitudes | Mean i = 20.430 and r = 22.635. Published i and r amplitude parameters are 1.774 and 2.060 mag. |
| Pan-STARRS1 DR2 | Object 144690235861403627 lies 0.03 arcsec from the CFHT position. PS1 measures i ≈ 18.50 on four exposures on MJD 56989 and 21.22 on MJD 56619. |
| ZTF DR24 | The nearest i light curve has a median position 0.07 arcsec from the Mira. Its 24 unflagged measurements run from i = 18.346 to 20.666 across 1,560 days. |
| ZTF public alerts | Fifteen positive r difference detections of ZTF19aawtdor in two seasons. The first nine span 92 days in 2020–2021. Six more span 54 days in 2023, alongside one g-band alert. Difference magnitudes refer to the excess over the ZTF reference. |

These measurements independently establish a large-amplitude red variable at the cataloged eruption position. Neither broadband i magnitudes nor ZTF r difference magnitudes measure the 1996 Hα flux. Mira pulsation shocks can produce Balmer emission, and a narrowband image admits stellar continuum. The strength of either contribution on the discovery date is unknown. [Balmer emission in Miras](https://ned.ipac.caltech.edu/level5/Gray/Gray29.html) · [PS1 detection columns](https://outerspace.stsci.edu/spaces/PANSTARRS/pages/298812397/PS1+Detection+table+fields) · [IRSA ZTF archive](https://irsa.ipac.caltech.edu/Missions/ztf.html)

## A comparison from the same survey

M33N 1997-09a, the other faint single-epoch Hα source (19.3 mag), lies 1.54 arcsec from the 523.87-day Mira candidate 01342428+3042433. SIMBAD already assigns the nova and the Mira identifier M33SSS J013424.31+304243.3 to the same object. This is a positive comparison with a larger offset. SIMBAD lists only the 2004 and 2012 nova papers for M33N 1996-12a in the saved cone result. Searches by both nova names, Mira identifiers, paper titles and the combinations with “Mira”, “symbiotic” and “long-period variable” found no earlier explicit link for the 1996 source. Catalogs and unindexed literature can contain further work. [M33 nova inventory](https://www.mpe.mpg.de/~m31novae/opt/m33/M33_table.html) · [Later nova compilation](https://doi.org/10.1088/0004-637X/752/2/156) · [Yuan et al. Mira survey](https://ui.adsabs.harvard.edu/abs/2017AJ....153..170Y/abstract)

The other four historical objects are more than 2.9 arcsec from the nearest gold-sample Mira. Controls shift the original six positions together through the local stellar background, retaining crowding differences among fields. See the checked-in JSON for control frequencies using the published gold selection and the full catalog. Of 100,000 shifted six-field groups, 696 contain a gold match at least as close as the 1996 source. Fifty-eight also contain a second gold match at least as close as the 1997 comparison. The respective counts using all table rows are 3,521 and 1,244. These comparisons are conditional on selection of this survey and these candidate properties. The original exploratory work had screened 62 nova inventory positions for ZTF alerts, examined a separate set of M87 nova positions, and considered several M33 candidates and catalogs. The six-field control therefore does not give a global discovery significance.

## Historical phase and competing explanations

A 543.85-day constant-period sinusoid fitted to Pan-STARRS1 and ZTF nightly i data places the 1996 date about 0.12 cycle before maximum. A free long-period solution near 546 days shifts that estimate. Half-period aliases, irregular sampling, selection at the detection limit, cycle-to-cycle changes and filter-dependent phase limit this extrapolation over more than 16 cycles. The historical phase does not establish that the Mira produced the Hα brightening.

An unrelated nova within the historical position error, or a white dwarf eruption in a symbiotic system with a Mira donor, remains possible. The survey found its faint candidates by inspecting aligned CCD frames at different epochs. Its published observing log is not a list of source-level nondetections. The original discovery and comparison CCD images and a contemporaneous spectrum were not obtained in this study. Relative registration of the historical residual against modern stars would test the blend alternative. A nova-line spectrum or a contemporaneous broadband flux history would help distinguish an eruption from a Mira maximum. The 1996 and 1997 detections were both below the original survey's 18-mag magnitude-limited threshold. We do not infer a revised M33 nova rate.

## Reproduction

Run `python analyze.py` in this directory with NumPy and SciPy. The script reads frozen CSV measurements, checks that they regenerate the scientific summaries and prints machine-readable output. The provenance file identifies the archival queries, row selections and checksums. The coordinate list uses the original six positions in Table 3. In particular the 1995-2 position differs from a later compilation.
