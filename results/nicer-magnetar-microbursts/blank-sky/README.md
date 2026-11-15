# Microsecond coincidences in seven NICER blank fields

A scan of **52,543,114 cleaned events in 3,231,475.85 seconds** of NICER blank-sky data found a maximum of **five hits in a nominal 10 microseconds**. The two published magnetar-direction flashes have seven and 32 cleaned hits in that window. No blank-field interval reached either multiplicity.

Four of the 28 non-overlapping five-hit windows triggered only the slow chain. Three occur within 16.24 milliseconds of one another in BKGD_RXTE_8 on 25 May 2023. The fourth is in a different field, BKGD_RXTE_5. All apparent energies in those four windows are below 2 keV, where fast triggering is less efficient. The two published flash windows each include at least one event above 2 keV. The blank-field observations establish a lower-multiplicity comparison and preserve a negative result at the actual flash multiplicities. The exact physical cause of each flash remains open. [Event energies and detector IDs](windows.csv) · [NICER timing chains](https://heasarc.gsfc.nasa.gov/docs/nicer/mission_guide/)

| Maximum hits in 10 μs | Number of observations |
| ---: | ---: |
| 0 | 24 |
| 1 | 1,851 |
| 2 | 240 |
| 3 | 41 |
| 4 | 15 |
| 5 | 7 |
| 6 or more | 0 |

## Selection and counting

The parent population is the 12 August 2025 NICER master-table snapshot, SHA-256 `7140f95ae198e98e0c1e7aac6ebca90a6829d4f325ae7d91d867cbed83ba4f90`. Select the exact target names BKGD_RXTE_1 through _6 and BKGD_RXTE_8 with positive catalogue exposure. This selects 2,179 ObsIDs. Of those, 2,178 yielded cleaned event files. The archive returned HTTP 404 for 0012050104 (166 seconds of catalogue exposure, public-date field MJD 88068). The stored observing time is the sum of the actual good-time interval lengths for each available file. The 24 files with zero in-band counts contribute their recorded exposure. [HEASARC master table](https://heasarc.gsfc.nasa.gov/FTP/heasarc/dbase/dump/heasarc_nicermastr.tdat.gz) · [Frozen ObsIDs, exposure, event counts, maxima and input hashes](observations.csv)

For each cleaned instrument event file, retain pulse-invariant channels 50 through 800 inclusive (0.50–8.00 keV). Sort the calibrated event times. Starting at every retained event, count all hits through the endpoint obtained by adding 10 microseconds. The largest width excess due to binary64 rounding among the threshold windows is 0.014 μs. The per-observation maximum includes every overlapping window. For the window list, retain the earliest five-hit window and resume after its final hit. Every listed group falls inside a single good-time interval and hit five different detectors. Fast triggering is verified against both the event-flag bit and the recorded fast pulse height.

The control uses archive pipeline cleaning, including the high-rate episodes present in those products. Chu et al. used their own background cuts and a Bayesian-block search and reported 3,308.9 ks of blank-sky exposure. This fixed-duration count answers a different question. Under a homogeneous Poisson model, zero groups at seven or more hits in the observed exposure gives a one-sided 95% upper rate of 0.93 per million seconds. Conditions vary across pointings, so transferring that bound to magnetar observations requires an exposure and background model. [Chu et al., §§2–3](https://doi.org/10.3847/1538-4365/ae5655)

## Reproduce

The first command checks the frozen tables using the Python standard library. The second re-downloads the seven files with threshold crossings, checks their SHA-256 digests, evaluates every possible starting event and verifies the detector-level records. It requires NumPy and Astropy. Use `--all` to rerun the complete search against the archive or `--obsid` to check an individual observation. The cache is optional.

```sh
python verify.py
python verify.py --threshold-files --cache /tmp/nicer-blank-fields
```
