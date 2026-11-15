# Apollo eclipse and deep-moonquake screens

The final Apollo 17 Grade-A file contains 13,000 arrivals for 3,250 events. We use the first arrival on any geophone. The six uneclipsed full Moons in the archive provide control windows with the same length and offset from the tabulated full Moon. Four adjacent windows have the same length as each eclipse window.

| Global penumbral contacts in UT | Events | Six phase controls | Four adjacent controls |
|---|---:|---|---|
| 1976 Nov 6 20:48:17 to Nov 7 01:14:05 | 4 | 1, 2, 6, 0, 3, 0 | 2, 4, 2, 5 |
| 1977 Apr 4 02:06:52 to 06:29:44 | 6 | 1, 2, 6, 0, 3, 1 | 3, 3, 6, 3 |

For each eclipse, one of the six uneclipsed control counts equals or exceeds the eclipse count. All ten cataloged events fall outside the 80–120° lunar-module direction used by Civilini et al. The April global umbral interval from 03:30:58 to 05:05:42 UT contains zero starts. The archive lists waveform packages covering both dates. Catalog counts alone leave recording duty cycle, detection efficiency and cooling at Taurus–Littrow uncertain. The contacts apply to the lunar disk. Earlier LSPE work recorded another eclipse. [Data and waveform inventory](https://doi.org/10.22002/d7mk4-f8t44) · [Catalog method and lunar-module direction](https://doi.org/10.1029/2022JE007704) · [Earlier LSPE report](https://ntrs.nasa.gov/api/citations/19760026027/downloads/19760026027.pdf) · [November contacts](https://eclipse.gsfc.nasa.gov/LEplot/LEplot1951/LE1976Nov06N.pdf) · [April contacts](https://eclipse.gsfc.nasa.gov/LEplot/LEplot1951/LE1977Apr04P.pdf) · [Full-Moon ephemeris](https://eclipse.gsfc.nasa.gov/phase/phases1901.html)

## Deep nests

The separate passive seismic catalog contains 7,052 events whose 2004 class is A followed by a positive integer, with 316 such nest identifiers. We count distinct-event, distinct-nest pairs up to four lags, including across midnight. For each lag, the statistic is the highest count among all unordered pairs of nests. Each of 4,999 randomizations keeps every event time and every nest's exact UTC-day count, shuffling only the assignment of nests to observed times within the day. This conditions on aggregate intraday bursts, observed gaps and daily tidal modulation.

| Maximum lag | Event pairs | Maximum count for one pair of nests | Null median maximum | Upper-tail p |
|---|---:|---:|---:|---:|
| 5 min | 6 | 1 | 1 | 1.0000 |
| 20 min | 227 | 3 | 3 | 0.7542 |
| 60 min | 1,150 | 8 | 8 | 0.6986 |
| 360 min | 7,934 | 56 | 54 | 0.3312 |

A1/A14 and A8/A18 tie at eight one-hour pairings. A8/A10 lead at six hours. The maximum-count test accounts for searching all labels at each lag and favors common nests. Recovery of overlapping seismograms and relative risk among rare pairs would require a different study. [Nakamura catalog and Panning transcription](https://github.com/mpanning/lunar_catalog/tree/ad8fea79c840e0291d779a41a7a5654404fc89eb) · [Deep-nest identification](https://doi.org/10.1029/2004JE002332)

## Reproduction

The script uses NumPy and pandas, and verifies each input file by SHA-256.

```sh
curl -fLsS 'https://data.caltech.edu/records/d7mk4-f8t44/files/GradeA_thermal_mq_catalog_final.csv?download=1' -o grade.csv
curl -fLsS 'https://raw.githubusercontent.com/mpanning/lunar_catalog/ad8fea79c840e0291d779a41a7a5654404fc89eb/levent.csv' -o deep.csv
python analyze.py grade.csv deep.csv > reproduced.json
```
