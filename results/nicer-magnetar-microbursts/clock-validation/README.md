# Published NICER burst starts have a 69-second timing discrepancy

A literal TDB interpretation of the calendar strings in Chu et al.'s machine-readable burst catalogue places eight ordinary bursts from eight different sources about **69.184 s early**. The strongest quarter-second interval in a blind three-minute search falls within 0.20 s of the shifted prediction for eight of eleven visits with event coverage. Six fall within 0.008 s. A twelfth selected visit has no good time in the test. The centre and width of the search were fixed before evaluating any of these sources. [Original machine-readable catalogue](https://content.cld.iop.org/journals/0067-0049/284/1/25/revision1/apjsae5655t4_mrt.txt) · [Catalogue paper](https://doi.org/10.3847/1538-4365/ae5655) · [All results and exposure times](results.csv)

The sign matters for time-coincidence searches. For these eight events the photons arrive later than predicted from the printed string read as TDB. Adding about 69.184 s to the catalogue time before converting it to spacecraft time brings the predicted position in the event list into agreement at this test's resolution. That offset equals TT minus UTC during the sampled years (2017–2024). The observations support a discrepancy in the representation of the published calendar times. They do not identify the software step that produced it, establish the absolute clock precision of these bursts, or independently check every catalogue entry. [NICER time-system and barycentering definitions](https://heasarc.gsfc.nasa.gov/docs/nicer/analysis_threads/time/)

## Independent sample and measurement

Start with the machine-readable table and exclude the two sources that prompted the hypothesis, 1E 1048.1−5937 and CXOU J010043.1−721134. For each remaining source select the highest reported blackbody flux among events with T90 between 20 and 300 ms, inclusive, and no other catalogue burst in the same observation within 150 s. This yields exactly one burst for each of twelve other sources. All twelve were kept, including a visit with zero current cleaned exposure. [Sample and measurement code](check.py)

Use the cleaned 0.5–8 keV event times, each observation's spacecraft orbit and its FITS RA_OBJ and DEC_OBJ. Treat the published string as TDB at the solar-system barycentre. Invert the Earth and spacecraft light-travel correction to find the predicted local arrival time. The calculation uses Astropy's built-in Earth ephemeris, linear spacecraft-position interpolation in metres, the event MJD reference, TT-to-TDB conversion and the FITS TIMEZERO. This approximation resolves a seconds-scale discrepancy. High-precision absolute timing also depends on HEASoft ephemeris and relativistic corrections and on the precise source direction. [Archive files, hashes and observation IDs](inputs.csv) · [HEASoft barycorr](https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/help/barycorr.html)

Search all possible event-anchored 250 ms intervals within ±90 s of the literal prediction. Record the start of the highest-count interval, taking the earliest in a tie. This start can precede the burst onset while still containing all its photons. The last two columns independently count events in 0.8 s windows centred on the two predictions. The count windows are fully exposed except where indicated. The rows are in source-name order.

| Source | Blind peak offset, s | Counts at printed time | Counts at +69.184 s |
| --- | ---: | ---: | ---: |
| 1E 1841−045 | +69.191 | 10 | 52 |
| 1E 2259+586 | +69.189 | 11 | 28 |
| 1RXS J170849.0−400910 | +50.189 | 17 | 119 |
| 4U 0142+61 | +69.179 | 62 | 109 |
| PSR J1846−0258 | +69.186 | 10 | 20 |
| SGR 1806−20 | +69.177 | 2 | 616 |
| SGR 1830−0645 | no exposure | — | — |
| SGR 1935+2154 | +69.192 | 1 | 1,337 |
| Swift J1555.2−5402 | +68.993 | no exposure | 86 |
| Swift J1818.0−1607 | +69.067 | 1 | 24 |
| Swift J1822.3−1606 | −66.136 | 2 | 39 |
| XTE J1810−197 | −36.945 | 187 | 127 |

Three covered observations put their strongest interval elsewhere. For 1RXS J170849.0−400910 and Swift J1822.3−1606, the shifted aperture still contains more photons than the literal aperture. For XTE J1810−197 the literal aperture contains more. These comparisons do not uniquely identify the intended burst for those three sources. Sky-coordinate errors, strong unrelated counts, the difference between a peak window and an onset, and evolving pipeline screening can affect individual checks.

## Two Fermi detections test the correction on a different spacecraft

Convert the literal and shifted predictions for all twelve selected bursts to **NICER-local UTC** and compare both lists to all 1,254 event-detection UTCs in Godwin et al.'s Fermi Gamma-ray Burst Monitor catalogue. With a 1 s coincidence window and no source prefilter, zero Fermi events match the literal times and two match the shifted times. Both GBM entries identify the same magnetar as NICER. SGR 1830−0645 is the source for which the current cleaned NICER visit has no good exposure around either prediction. Fermi supplies an independent check there. [GBM 17-year catalogue](https://magnetars.sabanciuniv.edu/sgr_catalog/) · [Godwin et al.](https://arxiv.org/abs/2609.02809) · [Full UTCs and signed residuals](gbm-matches.csv)

| GBM trigger | Source | GBM minus literal NICER prediction | GBM minus shifted NICER prediction | Expected geometric time difference | Residual after geometry |
| --- | --- | ---: | ---: | ---: | ---: |
| 210213417 | SGR 1830−0645 | +69.225 s | +0.041 s | +0.010 s | +0.031 s |
| 221110989 | SGR 1935+2154 | +69.192 s | +0.008 s | −0.005 s | +0.013 s |

The geometric term projects the difference between the geocentric spacecraft position vectors onto the FITS target direction. Positive values mean the wavefront reaches Fermi later. NICER positions are interpolated from its orbit file. Fermi positions are from the EVNTRATE record nearest the GBM catalogue event in the corresponding TRIGDAT file. The GBM catalogue reports the detection time, and its hardware trigger time can be different. Bandpasses, onset definitions, approximate ephemeris and source directions limit millisecond comparisons. The 13–31 ms residuals support the sign and scale of the large correction. They cannot measure an intrinsic interband delay. [GBM trigger bn210213417](https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/triggers/2021/bn210213417/current/) · [GBM trigger bn221110989](https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/triggers/2022/bn221110989/current/)

## All-catalogue GBM comparison

The adjacent [catalogue-wide NICER experiment](../catalog-clock/windows.csv) supplies literal spacecraft-time predictions for all 1,132 Chu et al. rows. Here those predictions are independently compared with all 1,254 GBM event UTCs at four fixed time offsets. NICER MET in the input already includes TIMEZERO. Use a TT reference epoch MJD 56658 + 0.000777592592592593 and compare local spacecraft clock times. Search the whole GBM list without filtering on source, then check its published most-probable source against the NICER target. All timing matches at the tabulated tolerances pass this source check. [Every matched pair within 1 s](catalog-gbm-matches.csv) · [Crossmatch code](catalog-gbm.py)

| Local time tolerance | −69.184 s | Literal TDB | +69.184 s | +138.368 s | +69.184 s outside April 2020 storm |
| --- | ---: | ---: | ---: | ---: | ---: |
| 0.10 s | 1 | 2 | 10 | 2 | 5 |
| 0.25 s | 4 | 4 | 20 | 5 | 7 |
| 1.00 s | 16 | 16 | 33 | 19 | 11 |

Entries count **distinct GBM detections**, keyed by GBM event MET. The 0.25 s shifted matches represent 20 separate NICER rows across eight visits, including 18 from SGR 1935+2154, one from SGR 1830−0645 and one from SGR 1806−20. Six have on-board GBM trigger numbers. The other 14 were detected in the archival GBM event search. At 1 s, two extra NICER rows match GBM detections already paired to other nearby NICER bursts, so 35 row pairs represent 33 GBM detections. This is a catalogue-level consistency check for the time convention. Different bands, instrumentation, event definitions, approximate barycentering and the light travel time between spacecraft limit the precision of individual pairings.

The single dense 28 April 2020 NICER burst storm (ObsID 3020560101) produces every coincidence at all three control offsets even when the tolerance is widened to 1 s. Outside that visit, with a 0.25 s tolerance the +69.184 s offset retains seven detections in seven different observations across three magnetars. The three control offsets each retain zero. Six of these seven have GBM triggers: 210213417, 210914259, 230228885, 221110989, 221114541 and 221014487. The seventh has GBM event MET 687289396.031. Their observed Fermi-minus-shifted-NICER differences span 0.0003–0.1185 s before the inter-spacecraft propagation correction. The two individually checked triggers above demonstrate its scale. GBM time comparisons do not require a currently valid NICER good interval. Bursts missing from the GBM catalogue supply no timing constraint.

## Independent published times

Borghese et al. separately tabulated seven NICER burst epochs at whole-second TDB precision for SGR J1935+2154. Search every one for the nearest Chu et al. entry from the same ObsID within ±90 s. Two match, with differences that agree with the photon test. The other five have no match under this rule. These are independent analyses of shared NICER observations and this comparison uses the reported times directly, without an orbit calculation. [Borghese et al. 2022, Table 4](https://academic.oup.com/view-large/370754906) · [Seven original epochs](published-times.csv)

| ObsID | Earlier paper's TDB epoch | Chu event | Earlier paper minus Chu | Residual after adding 69.184 s to Chu |
| --- | --- | ---: | ---: | ---: |
| 3655010201 | 2020-05-18 13:39:39 | 358 | +69.513 s | +0.329 s |
| 3020560106 | 2020-05-20 14:06:36 | 351 | +69.247 s | +0.063 s |

The scope of this finding is the printed absolute timestamps. A time-scale offset common to an event's photons would cancel in a duration measurement. On 14 November 2026, searches combining the catalogue title with “69.184”, “TDB UTC”, “correction”, “erratum” and “burst times” located the primary catalogue and mission documentation, without locating a published catalogue correction.

## Replay

Use Python 3.10 or later with NumPy, Astropy and lxml. From this directory, the following command retrieves approximately 55 MiB of primary data, verifies the SHA-256 checksums, repeats the full source selection and all three comparisons, and checks the resulting measurements against the frozen table. It accepts an existing cache. The offsets use spacecraft seconds relative to the literal TDB prediction. The tested alternative adds the constant locally. Differences between adding the constant there and adding it before astronomical conversion are below 0.008 s for these dates and well inside the 0.8 s counting aperture.

```sh
python check.py /path/to/nicer-clock-cache
python catalog-gbm.py ../catalog-clock/windows.csv /path/to/nicer-clock-cache/gbm-catalog.html
```
