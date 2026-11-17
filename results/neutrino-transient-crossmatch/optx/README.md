# The independent inputs contain 73 optical flare episodes and 362 X-ray triggers

The two matching tables are `optical_events.csv` and `xray_events.csv`. Coordinates are degrees, and epochs are UTC. Catalogue construction used transient properties, mission classifications and observing metadata. No neutrino location or time was used to select a row. The GRB case directory contains post-match follow-up.

## Optical sources

| Input | Frozen coverage and selection | Contribution |
|---|---|---|
| TNS through the LINCC HATS mirror | All 12 partitions of the 2026-11-17 01:01 UTC collection. Checked 213,398 unique object IDs against the publisher row count. Exact TDE, TDE-H, TDE-He, TDE-H-He, TDE-featureless classifications with a 2024–26 discovery date. | 69 distinct TNS-classified TDEs. ZTF identifiers/group names appear in 63 reporting rows, ATLAS in 41. |
| Langis et al. TDECat | All 134 objects in the publisher CSV at Git revision bca013a266cd5d2a120aa4c6b34a2dc54bee1e5e. Projected 22 in the date range. | All 22 already occur in TNS. Discovery UTC available. |
| ZTF Bright Transient Survey | Full TDE sample query, unrestricted by magnitude in the query. The BTS selection and spectroscopic follow-up vary with time. Rows require a 2024–26 discovery. | 22 published peak timestamps. Epoch is JD−2458000, converted to UTC. |
| Z+LIONS | All 1,425 rows of the survey's magnitude-limited ZTF/Rubin sample served on November 17. Two are classified TDE. One supplied peak predates 2024, one falls in range. | One additional peak timestamp. Z+LIONS covers part of Rubin's footprint and excludes July 2026 from its published sample criteria. |
| Quintin et al. Fink/ZTF | Joined paper sky positions to fitted peak dates. Kept its TDE-candidate or ANT-candidate group with a 2024–26 peak and excluded ambiguous events and candidate SNe. | Four independent optically detected nuclear flares with classified candidate origins, unsampled peak errors 1.0–5.2 d. |

TNS names and source aliases are preserved in the selection projections. One row per episode, using a published survey peak before discovery time: 46 discovery epochs, 23 observed peaks, four fitted peaks. The four day-resolution peak dates are encoded at 12:00 UTC and preserve their uncertainties in `nuclear_events.csv`. No two retained rows are within 1 arcsec. The two different peaks of one nucleus are identified by episode when relevant. One LIONS host showed an archival optical maximum years before its TNS event, so that archival peak did not supply an in-range epoch. Nuclear flares with ambiguous classification, unreported flares, and unclassified or proprietary Rubin/ATLAS alerts do not enter the catalogue. The TNS snapshot precedes the cutoff by eleven hours. The BTS dynamic response was served about four minutes after cutoff and carries no independent snapshot-creation timestamp.

[LINCC catalogue](https://data.lsdb.io/hats/tns/2026-11-17/) · [TNS](https://www.wis-tns.org/content/tns-getting-started) · [Langis et al., TDECat, CC BY-NC 4.0](https://github.com/dlangis/TDECat) · [ZTF BTS](https://sites.astro.caltech.edu/ztf/bts/explorer_info.html) · [Z+LIONS](https://sites.astro.caltech.edu/LIONS/lions.html) · [Quintin et al.](https://arxiv.org/abs/2511.19016)

## X-ray sources

| Input | Selection | Retained |
|---|---|---:|
| Einstein Probe WXT | Autonomous WXT initial detections by the EP team account in NASA GCN. Search subjects starting with an EP designation and an instrument detection term. First WXT timestamp and J2000 position must parse, agree with designation date and precede follow-up. | 192 of 205 names |
| Swift XRT LSXPS | Discovery-date selection in frozen 11:00:01 UTC dump, classifications 1/3 (confirmed/outburst). Joined all 83 to first-detection times and discovery observation IDs on LSXPS detail pages, then to target names and pointings in the HEASARC Swift master release dated November 14. | 34 of 83 |
| Swift BAT | Confirmed Swift GRB-table row with GCN flight BAT position notice and explicit on-board rate/image trigger; independently supplement from the detailed mission BAT catalogue where trigger method is rate/image. | 140 |

LSXPS discovery processing time can differ from the time of photon acquisition. The XRT sample uses the public first-detection observation epoch. It excludes 47 fields with external transient targets or ambiguous target names and two further sources within 2 arcmin of the pointing centre. Some included sources are Galactic. These restrictions intentionally omit many targeted follow-ups, and Swift pointing history makes the coverage spatially and temporally uneven.

The EP inventory excludes seven names with corrections or retractions requiring review, four ambiguous initial/follow-up localizations and two nonstandard packets. Other subject formats, public alerts outside GCN or unreported triggers can escape this extraction. The terminal daily archive contains circulars through 46308. The contemporaneous GCN index and four newer individual records were checked through the cutoff. They contain no additional initial WXT trigger. GCN 46312 refined a known burst light curve and retained the independent trigger epoch. Historical revisions of public tables do not necessarily preserve their original issue dates.

For cross-mission aliases, four pairs within 0.25° and seven days have an explicit joint EP/GRB identification in their EP circular. Their separations are under 0.03° and trigger offsets under 40 s. Each is counted once using the BAT mission row. No other pair of the three input lists met that screen. `xray_aliases.csv` contains the GCN IDs. The last verified BAT burst occurred on October 17, 2026, the last included XRT detection on October 18, and the last EP trigger on November 15. Search sensitivity and sky coverage for the public inputs are heterogeneous.

[GCN archive](https://gcn.nasa.gov/circulars) · [LSXPS transients](https://www.swift.ac.uk/LSXPS/transients/) · [Hourly LSXPS tables](https://www.swift.ac.uk/LSXPS/oldTables.php) · [HEASARC Swift pointing archive](https://heasarc.gsfc.nasa.gov/W3Browse/catalog/swiftmastr.html) · [BAT mission catalogue](https://swift.gsfc.nasa.gov/results/batgrbcat/) · [Swift GRB table](https://swift.gsfc.nasa.gov/archive/grb_table/fullview/) · [GCN flight trigger archive](https://gcn.gsfc.nasa.gov/swift_grbs.html)

## Reproduction

The projected CSVs are the frozen scientific inputs. The raw upstream byte hashes, URLs, HTTP validators when provided and retrieval times are in `receipts.json`. `bat_inventory.csv` records SHA-256 for each flight-notice payload. `data/swift_join.csv` records the hash of each discovery page. The working GCN archive, Swift pointing master, all-sky parquet partitions and full mixed-type tables are represented by receipts plus the frozen relevant projections. Upstream rolling copies and mirrors can change or expire. A digest mismatch is reported on reacquisition.

Use Python 3.11+, lxml, pandas and fastparquet for full projection from raw data. In this directory, replace CACHE below with an absolute path. For Swift detail acquisition, run `swift_pages.py` in CACHE, then run `swift.py CACHE` from this directory. `swift_pages.py` expects the frozen file `swift11.gz`.

```sh
python fetch_receipts.py CACHE
python project_optical.py CACHE
python nuclear.py CACHE/nuclear
python optical.py
python ep.py CACHE/archive
python swift.py CACHE
python bat.py CACHE fetch
python combine_xray.py
sha256sum -c SHA256SUMS
```

To regenerate the analysis-facing tables using only the frozen projections, run `python optical.py`, `python swift.py`, then `python combine_xray.py`. The catalogue scripts preserve inclusion/exclusion reasons and do not run any neutrino association test.
