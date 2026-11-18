# Search for detached emission in GRB 261115A

A screen of Fermi GBM CTIME for the November 15, 2026 trigger found **no qualifying detached count excess** in the selected pre-burst and late intervals. The two best-facing NaI detectors had no coincident threshold crossing in **2,381** windows. The known prompt emission was recovered by the two primary detectors and the third source-facing detector. This null concerns one burst, one pulse-height band, and the intervals and widths below.

## Observed result

Times are seconds from the Fermi trigger at 22:10:53.44 UTC. All complete native bins in the pre and late intervals passed the joint quality and exposure cuts. The windows slide by one native bin and include four nominal widths from 0.256 to 2.048 s. The first two rows form the primary test family. The other two are separate time controls scored using the primary trial count.

| Region (s) | Complete native bins | Windows | Joint threshold crossings | Lowest unadjusted joint p | Interval at the minimum (s) |
| --- | ---: | ---: | ---: | ---: | --- |
| −32 to −3 | 113 | 441 | 0 | 0.06246 | −6.336 to −4.288 |
| +8 to +40 | 499 | 1940 | 0 | 0.04045 | +17.536 to +19.584 |
| −106 to −77 | 112 | 437 | 0 | 0.01546 | −95.938 to −95.682 |
| +77 to +109 | 499 | 1940 | 0 | 0.13089 | +101.250 to +101.762 |

Both primary minima adjust to 1. The strongest late interval has 1,125 counts in n9 and 1,048 in n6. The separate source-facing n7 has p = 0.163 there and off-axis n4 has p = 0.043. Neither the third-detector rule nor the directional veto determined the disposition because no primary pair passed. The independently reported GECAM-B prompt emission sits within the excluded −3 to +8 s guard. [Fermi GBM localization](https://gcn.nasa.gov/circulars/46299) · [GECAM-B timing](https://gcn.nasa.gov/circulars/46305)

The fixed prompt control (−1 to +4 s) records **5,274** counts in n9 and **4,758** in n6, both over about 4.90 s of live exposure. Both rates are roughly twice their before-band rates. The third detector n7 records 4,770 counts. The n5 directional diagnostic also brightens in the known prompt control, so its viewing angle alone cannot establish a pure background reference. In the primary detectors the before-band rates are about 2–3% higher than the after-band rates. Each reported window uses the more conservative of its two sideband comparisons. Subinterval rates and all per-window results are deposited. [Search protocol](PROTOCOL.md) · [Summary and coverage](screen.json) · [All windows](windows.csv.gz) · [Recomputed claims](claims.json)

## Count-space calibration

In a separate calculation with the observed sideband counts fixed, we inverted each count test to estimate the mean signal needed to pass in n9, n6 and n7 at a particular covered window. The model adds the same mean count to each detector and models each detector's counts with a Poisson distribution whose background mean uses the lower measured sideband rate. A union bound with each detector passing at least 96.67% of the time gives at least 90% joint probability without a between-detector independence assumption. This requires about **86** added mean counts per detector for a 0.256 s window and **250–251** for a 2.048 s window across the two intervals. The calculation omits the off-axis veto and assumes the signal falls inside the window. It gives no flux or fluence limit, source completeness, or coverage of other times, durations and energies. Real signals produce unequal detector counts. [Sensitivity calculation](sensitivity.py) · [Sensitivity output](sensitivity.json)

## Inputs and reproduction

The five trigger CTIME v00 files have 8 pulse-height channels each. This analysis uses channels 2 through 4, with detector-specific bounds from about 26 to 296 keV. The TCAT v00 geometry was used for detector selection. All six files were retrieved from the official [HEASARC trigger directory](https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/triggers/2026/bn261115924/current/) and their exact identities are in [inputs.json](inputs.json). The full input set totals 2,986,560 bytes and is stored on D outside Git. The archive can revise products, so the scripts check hashes before parsing any data.

Protocol and primary code were frozen at `dcc7784c5683629297c2cdf2baef85bbd75fd0fe` before opening count columns. The first run rejected four zero-duration rows per file at +597.976 s before reading counts. Commit `150762184b8684306c51b482e20a877e2c93fab5` excludes that padding, which cannot intersect any search, control or background interval. Commit `6e441e417a207e1b8ded062e6e919a9b6d4c0617` fixes the JSON type of the native-bin count. Neither repair changed the protocol or thresholds. The independent window builder reran all **4,758** signal and control windows. Counts, exposures and conservative p-values matched exactly across all five detectors. Both implementations use the same exact binomial tail routine and the same underlying Fermi measurements. [Independent window builder](replay.py) · [Agreement receipt](replay.json)

Use Python 3.14.0, NumPy 2.5.5, SciPy 1.17.1 and Astropy 8.0.2. Download the six files named in the manifest into an input directory and run the following with a fresh output directory.

```text
python screen.py --inputs RAW_DIR --output RUN_DIR
python sensitivity.py RUN_DIR RUN_DIR/sensitivity.json
python replay.py RAW_DIR inputs.json PEER_DIR --compare RUN_DIR/windows.csv
```

Check the deposited [file digests](SHA256SUMS). The uncompressed `windows.csv` SHA-256 is `08067f7aa48f9364bf917868a703ebe85d9cdd1586547644aab92a5e58a9cc81`. The [GBM product documentation](https://fermi.gsfc.nasa.gov/ssc/data/access/data_products.html) and [instrument caveats](https://fermi.gsfc.nasa.gov/ssc/data/analysis/GBM_caveats.html) describe the time and exposure data. We excluded cadence transitions from all scanned windows. CTIME exposure estimates, background curvature, and same-spacecraft backgrounds limit the count model. This work did not analyze TTE events or independent GECAM-B counts.

An event-alias search on November 18 UTC checked `GRB 261115A`, `261115924`, precursor and extended-emission queries. The [event circulars](https://gcn.nasa.gov/circulars/events/grb-261115a?endDate=&query=&startDate=&view=index) document the known burst and optical follow-up. The generic report of multiple prompt pulses predates this work. [Published GBM precursor searches](https://arxiv.org/abs/2004.03246) establish prior use of these data for finding separated emission. The result here is a bounded negative check on this trigger.
