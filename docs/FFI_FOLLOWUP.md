# Resolved-source checks for the LMC image search

Six of the 13 recurring TESS full-frame signals have a nearby Gaia star whose independently measured G-band period matches both TESS sectors. The published OGLE catalog comparison identifies period-compatible variables near 12 signals, allowing half- and double-period matches. **One signal in LMC field 4 remains a source-identification lead.** It recurs at roughly 1.11 days in a crowded field that contains several known variables. [TESS image measurements](../results/ffi-lmc-grid.json) · [OGLE crossmatch](../results/ffi-ogle.json) · [Gaia epoch-period measurements](../results/ffi-independent.json)

## Direct Gaia period matches

We queried every Gaia DR3 source within 30 arcsec of each discovery centroid and searched every available G-band epoch series. The 13 cones contain 2,990 catalog entries in total. Only 32 distinct sources have released epoch photometry in these cones. Each search uses a weighted floating-mean sinusoid, the G-band flux errors, the survey's rejection flags and periods from 0.15 to 8 days. The comparison uses the strongest peak over that entire interval. Qualifying matches require at least 20 accepted transits over 100 days, at least 40% fitted weighted variance, a Gaia position within 15 arcsec of the S1 center, and a frequency within one inverse-sector-baseline of both TESS frequencies. These are screening criteria. The analysis does not assign a survey-wide false-discovery probability or exclude an unresolved blend. [Gaia DR3 time-series data model](https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_photometry/ssec_dm_epoch_photometry.html) · [Period-search implementation](../examples/ffi_gaia_followup.py)

| TESS centroid (ICRS degrees) | TESS S1 / S2 period (d) | Gaia period (d) | Distance (arcsec) |
| --- | ---: | ---: | ---: |
| 76.99986, −69.40131 | 0.67234 / 0.67148 | 0.671175 | 1.32 |
| 77.41663, −69.34141 | 0.99317 / 0.99507 | 0.995202 | 1.29 |
| 77.81863, −69.36103 | 1.18311 / 1.18337 | 1.182511 | 4.06 |
| 76.90371, −69.21339 | 1.22207 / 1.21624 | 1.217908 | 2.44 |
| 77.08309, −69.17205 | 1.80715 / 1.80950 | 1.812336 | 4.89 |
| 77.46720, −69.20884 | 1.94953 / 1.93765 | 1.942869 | 2.87 |

Gaia samples these fields sparsely. Its strongest frequency may be an alias or harmonic, particularly for eclipsing systems with two minima. For example, the 1.55-day TESS signal at (77.79707°, −69.20872°) has a Gaia star 1.5 arcsec away with a 3.128-day photometric period. It falls outside this direct-match count. [Complete catalog cone receipts and Gaia photometry digests](../results/ffi-independent.json)

## Unresolved LMC-4 signal

The discovery centroid is (77.0550166°, −69.2370982°). The second-sector centroid is 6.30 arcsec away. The TESS periods are 1.112209 and 1.116938 days and their peak-pixel fitted variance fractions are 0.32 and 0.69. The discovery and validation brightest variable pixels both lie near an edge of their cutout. A larger cutout or a different roll angle can improve localization. The 30-arcsec Gaia cone has **238 catalogued stars**, five brighter than G=16. The 21-arcsec TESS sampling leaves multiple plausible emitters within one pixel. [Image source and per-sector diagnostics](../results/ffi-lmc-grid.json) · [Gaia catalog counts](../results/ffi-independent.json)

| Source | Distance to S1 / S2 centroids | Independent observations |
| --- | ---: | --- |
| OGLE-LMC-LPV-20399 | 2.22″ / 6.26″ | OGLE lists 19.963, 197.32 and 19.706 days. A blind 0.15–500 day search of 756 I-band epochs over 4,228 days recovered 19.962 days. The exact TESS frequencies fit only 0.44% and 0.10% of weighted variance. |
| OGLE J050812.41−691410.1, Gaia DR3 4661229089008619008 | 5.82″ / 10.90″ | This is the only Gaia source in the cone with released epochs. The best 0.15–8 day G-band fit on 30 measurements over 911 days is 4.355 days (62% of weighted variance). The exact TESS frequencies fit 9.65% and 12.48%. Sabogal et al. reported brightness-state changes in OGLE-II data. Gaia places it in the combined Be/GCAS/SDOR/WR classifier group. |
| UCAC2 1797921, MACHO 1.4532.22, OGLE J050813.83−691422.9 | 9.70″ / 4.14″ | Sabogal et al. reported a Be-star candidate. In 270 OGLE-II I observations, subtracting a running 24-day median makes a weak 1.123-day feature the strongest short-period peak. It explains about 8% of weighted variance and neither chronological half selects it as its strongest peak. |
| OGLE-LMC-RRLYR-05589 and OGLE-LMC-HB-0108 | 13.41″ / 12.30″ and 22.38″ / 18.00″ | OGLE reports 0.247605 and 220.437 days, respectively. Neither is period-compatible under the original OGLE match rule. |

The original MACHO red and blue light curves resolve UCAC2 1797921 separately. Both have measurements spanning 2,723 days, with 1,212 accepted red epochs and 1,214 blue epochs. After the same detrending, neither global maximum agrees with the TESS period. The strongest local peak in 1.07–1.16 days explains under 3% in each band. Neither band recovers the OGLE-II 1.123-day feature strongly. The two MACHO bands sample the same observing program, and the ground-based programs predate TESS. Evolving variability remains possible. The evidence supports a recurrent TESS image signal with an unresolved emitter and gives no basis for claiming a newly discovered star. [Ground-based measurements and archive receipts](../results/ffi-independent-ground.json) · [OGLE-III LPV catalog and photometry](https://ogle.astrouw.edu.pl/ogle/ogle3/OIII-CVS/lmc/lpv/) · [Published Be-star source table](https://cdsarc.cds.unistra.fr/viz-bin/ReadMe/J/MNRAS/361/1055?format=html) · [OGLE-II difference photometry](https://ogle.astrouw.edu.pl/ogle/ogle2/dia/) · [MACHO archive](https://macho.nci.org.au/)

## Reproduce

The following commands operate on the 13 published recurrent positions. They retain URL-keyed response caches and put the source SHA-256 values in the emitted JSON. The ground-based script also hashes each individual light curve extracted from the OGLE-II field archive. The scoped MACHO analysis uses the survey's per-band quality flag, reported magnitude errors below 0.08 mag and a six-MAD outlier cut. The short-period search reports both original magnitudes and magnitudes after subtracting the median of each centered 24-day interval. Local 1.07–1.16-day comparisons are post hoc. They do not constitute an independent blind detection.

```sh
python -m pip install -e '.[archive]'
python -m examples.ffi_gaia_followup results/ffi-lmc-grid.json --cache runs/ffi-gaia > gaia-check.json
python -m examples.ffi_ground_followup results/ffi-lmc-grid.json --cache runs/ffi-ground > ground-check.json
```
