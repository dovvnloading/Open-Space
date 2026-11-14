# Search for additional transits around 90 TESS candidate hosts

Across 393 two-minute TESS SPOC sector light curves, three of 630 periodogram peaks passed a retrospective repeatability and catalog-harmonic screen. Two track the evolving brightness modulation of V546 Car. The third reverses sign in simple-aperture photometry. **No additional transiting companion advanced in this search.** The result covers this sample, period range and selection procedure. It does not constrain the occurrence of additional planets.

## What the surviving peaks did

| Target | Fitted periods | Measurement at the fixed discovery ephemeris | Disposition |
| --- | --- | --- | --- |
| V546 Car, TIC 460950389, TOI-6715 | 1.10428 and 6.63196 d | Both are harmonics of the stronger 2.21068 d modulation. At 2.21068 d the flux deficits in sectors 63 and 64 have nominal signal-to-noise ratios of 48.4 and 26.7. The same phase in sector 90 has a ratio of −3.3. At half and three times the period the broader signal leaks into the held-out box. | Broad stellar variability, with changing amplitude and phase. The aperture shows percent-level modulation across much more of the cycle than the fitted box. |
| TIC 208297760, TOI-4417 | 7.12996 d | The sector 66 corrected curve shows a 649 ppm deficit at a nominal ratio of 9.84. Simple-aperture photometry at the same times gives −553 ppm, a brightening. The other held-out sector, 93, gives 24 ppm in the corrected curve and −14 ppm in the aperture. | Fails the aperture and replication checks. |

A separate high-ranking case illustrates a background false alarm. The TOI-5697 sector 44 dip at BTJD 2514.5623 drives a peak near 10.78571 d. During this event the aperture flux drops 0.82%, the background rises 3.96% and the first-moment column centroid moves −0.00780 pixel. Other predicted events lack a consistent deficit. The values use the median within 0.06 d of the event and flanking intervals 0.23–0.4 d away. The bundled control script repeats this measurement.

As a catalog control, TOI-6699.01 recurs at 5.32019 d with a fitted depth of 847 ppm. Its catalog period is 5.32200 d. It survives the epoch mask but the period comparison flags it. A known period can remain in a light curve when a catalog transit window drifts.

## Sample and search

On 2026-11-14 the NASA Exoplanet Archive TOI table contained 8,175 rows. The newest release date in that download was 2026-11-07. A row qualified with disposition PC, stellar radius below 0.85 solar radii, TESS magnitude between 0 and 12.5, period below 15 d and reported planet radius below 5 Earth radii. Sorting matching rows by decreasing numeric TOI and keeping the first row per TIC gave 277 distinct hosts. The first 150 formed the target list. A query to MAST CAOM supplied two-minute TESS SPOC light curves. Ninety of the 150 had consecutive sectors usable by this search. The manifest contains 7,756,957 cadences before quality, binning, clipping and transit masks. Selection by TOI number favors recently assigned numbers. It is not an unbiased planet survey.

For each host the initial selection took up to five sectors from the longest consecutive run (latest run in a tie), plus up to two latest remaining products. The exact mapping and products are frozen. On these selected sectors the search used the longest consecutive run, breaking ties by the latest sector. When all selected sectors formed one run with three or more sectors, the final sector was held out. Eight hosts have only the two search sectors. They could yield exploratory peaks but could not pass the held-out requirement.

The search removes nonzero quality flags, averages groups of at least three cadences into 10-minute bins, and divides each segment by a running median of up to 145 bins. Segments split across gaps longer than 0.12 d. Finite residuals between −15 and +5 robust scatter units are retained. All published TOI ephemerides for these TICs are masked with a half-width of one catalog duration, floored at 0.035 d and capped at 0.14 times the period. A missing duration defaults to 0.09 d.

A box least squares (BLS) period search uses a frequency grid corresponding to roughly 0.55–13 d. A local refinement can pass slightly beyond the grid endpoints. The trial durations span roughly one to 3.5 h. The seven strongest separate frequency peaks per host are retained, including aliases. The rankings reflect coarse search power before refinement. Peak depths are fractions of the median-normalized stellar flux. Epochs are BJD in the TDB time system minus 2,457,000.

Retrospective triage requires a BLS depth-to-error ratio of at least 8, two training sectors each with a local ratio of at least 3, and at least one held-out sector at 4 or higher at the unchanged period and epoch. Local depths compare the mean in the central box with shoulders at 0.8–2.5 times its full duration from center. The error assumes independent photometry. These ratios do not assign a false-alarm probability and correlated noise can strongly inflate them. A peak is flagged when its ratio to any known TOI period lies within 0.008 times the larger of one and that ratio of any of 1/4, 1/3, 1/2, 2/3, 1, 3/2, 2, 3 or 4. Sixty-four of 630 peaks carry this flag. Remaining triage survivors receive a separate simple-aperture depth check. An unchanged short-baseline ephemeris may miss a real transit many sectors later. No injection-completeness measurement, statistical validation or source localization was attempted for the rejected peaks.

## Reproduction

The snapshot of the original TOI CSV had SHA-256 `c8e6455e6694af2f8e38a969dbaa56f93ba97f5d5976606ce4c4e780dff62002`. The frozen ephemeris subset and selection table avoid a dependency on future archive edits. Each manifest row records the direct public MAST cloud product URL, original cadence count and SHA-256 of eight ordered measurement columns. The fetcher normalizes floating-point NaNs and serializes each column as little-endian 64-bit values for a digest that is independent of FITS container metadata. The raw FITS files are read on demand. [NASA TOI table](https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+%2A+from+toi&format=csv) · [TESS data products](https://archive.stsci.edu/missions-and-data/tess/data-products) · [V546 Car in SIMBAD](https://simbad.cds.unistra.fr/simbad/sim-id?Ident=V546+Car)

The measured environment used Python 3.12, NumPy 2.3.5, SciPy 1.17.0, pandas 2.2.3 and Astropy 8.0.2. Run in this directory with an empty cache or an existing matching cache.

```sh
python fetch.py /path/to/cache
python screen.py /path/to/cache /path/to/replay
python control.py /path/to/cache > /path/to/replay/background-control.csv
diff peaks.csv /path/to/replay/peaks.csv
diff sector_checks.csv /path/to/replay/sector_checks.csv
diff background-control.csv /path/to/replay/background-control.csv
```

`selection.csv` records all 150 queried TICs and the sectors available and used. `ephemerides.csv` contains the mask inputs. `products.csv` identifies the 393 source products. `peaks.csv` contains all seven peaks for each host, including unfavored outcomes. `sector_checks.csv` contains the per-sector depth at each peak, its white-noise ratio and the number of epochs with any in-box coverage. Negative depths are brightenings. NaN indicates inadequate in-box or shoulder coverage. `background-control.csv` preserves the raw instrumental diagnostic. The top-level manifest covers the files needed to replay this result.
