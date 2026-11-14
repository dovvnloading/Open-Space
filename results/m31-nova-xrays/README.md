# Two candidate X-ray counterparts of optical novae in M31

Chandra measured a point source 0.246 arcsec from the Fe II nova M31N 2017-06c more than four years after its optical discovery. A second point source appeared 0.178 arcsec from M31N 2021-06b about four months after discovery. Both were bright in two October 2021 visits. By the 2023 and 2024 visits, counts at both optical positions had fallen close to the local background. The X-ray sources themselves appear in the Chandra Source Catalog (CSC). This study proposes the two optical associations and measures the X-ray excesses in the discovery and later comparison visits.

| Optical eruption | CSC 2.2 source | Optical-to-X-ray delay in the two visits | Oct 2021 HRC counts / local background | 2023–24 HRC counts / local background |
|---|---|---|---:|---:|
| M31N 2017-06c | 2CXO J004232.8+411223 | 1,592.28 and 1,594.97 days | 120 / 27.50 | 32 / 24.58 |
| M31N 2021-06b | 2CXO J004258.9+411808 | 128.28 and 130.97 days | 140 / 25.50 | 33 / 30.67 |

![Exposure-normalized Chandra HRC photon images at both optical positions in 2021 and 2023–24](hrc-comparison.png)

The circles mark the optical apertures. All four panels share a count-rate scale and 0.5-arcsec Gaussian smoothing. The full lists and aperture calculations use unsmoothed photons. [Figure code](figure.py)

The two 2021 visits provided 52.70 ks of exposure. The nine later HRC visits provided 98.77 ks. Every count in the table comes from a circular 2-arcsec aperture centered on the **optical** coordinates. The 4–8 arcsec annulus has 12 times the source area. The two individual visits show positive excesses at each target. Conservative 95% lower bounds on the decline in HRC aperture count rate are 4.79-fold for 2017-06c and 7.59-fold for 2021-06b. These bounds use one-sided exact Poisson limits on four on/off means with a Bonferroni correction, a locally uniform background in each epoch, and comparable HRC response at these nearly identical pointings. Four ACIS pointings in 2024 independently gave 3 and 4 aperture counts, with expected backgrounds of 3.42 and 3.33 counts, respectively. The event lists cover the same sources in two older HRC control fields. [All measurements](measurements.txt) · [Official file identifiers and hashes](observations.csv)

A late X-ray phase can accompany nuclear burning on a white dwarf after a nova. A four-year delay makes 2017-06c particularly interesting. The observations establish the X-ray state in October 2021 and leave the intervening years unmeasured. HRC has limited energy resolution. A separate transient at the same position, or a later unobserved eruption, could contribute. The optical classifications differ. M31N 2017-06c has spectroscopy, while 2021-06b has optical and ultraviolet measurements. [2017 discovery](https://www.astronomerstelegram.org/?read=10484) · [2017 spectroscopy](https://www.astronomerstelegram.org/?read=10487) · [2021 optical discovery](https://www.astronomerstelegram.org/?read=14702) · [Published M31 nova X-ray monitoring](https://doi.org/10.1051/0004-6361/201322426)

## Search, competing sources and controls

An exploratory cross-match of optical novae, XMM-Newton and Chandra led to the two October pointings. For the defined October test, we took all 325 optical eruptions with JD 2456000 through 2459505.78 in the November 13, 2026 living catalog. We matched them against the 100 distinct CSC sources with HRC wide-band flux significance of at least 3 in **both** ObsIDs 23871 and 26163. This is a retrospective sample. The optical catalog includes M31 candidates outside the observed field. Each mission centroid was compared with every eligible optical position on the unit sphere. The matching radius was 2.5 arcsec. [Optical input](nova_positions.csv) · [X-ray inputs](epoch_detections.csv)

Six pairs qualified, covering five optical eruptions. Four of these pairs had significant X-rays in earlier observations before the corresponding optical eruption. This includes the apparently close 0.425-arcsec match at 2017-01a, which had significant X-rays in 103 pre-eruption visits. The two sources reported here had zero pre-eruption observations reaching a flux significance of 3 in any reported band, over 178 and 180 CSC extractions. The final two detections in each history are the October 2021 visits. Differences in instrument response and exposure affect the historical limits. [Source histories](source_history.csv.gz)

For a crowding control, we shifted the whole optical catalog through 36 directions at each of 5, 10, 20 and 40 arcsec. Three optical positions lie within 0.5 arcsec of a doubly detected source in the real field. The maximum was one in the 144 shifted fields. At 2.5 arcsec, one shifted field matched as many novae as the original. These correlated shifts describe local chance alignments and do not give a posterior probability for either association. The fixed-aperture binomial tails in the replay describe a local background test assuming locally uniform surface density.

A weaker case, M31N 2019-07c, shows photon excesses in the optical aperture but fails the two-visit catalog threshold. Its second wide-band flux significance is 2.77, and one pre-eruption CSC observation reaches 3.65. The released cutouts include this case so that its exclusion can be checked.

Searches through November 14, 2026 used both nova names, optical-coordinate aliases, both 2CXO names, and both Chandra observation IDs. The Pietsch table updated in 2022 records optical spectroscopy or UV follow-up at these positions. The later living catalog lists the optical eruptions. We found no paper or telegram reporting these two optical-to-X-ray associations in the checked sources. The CSC positions and flux extractions are public mission results. [Pietsch catalog](https://www.mpe.mpg.de/~m31novae/opt/m31/M31_table.html) · [Living catalog](https://ashafter.sdsu.edu/m31_nova_table.html) · [Chandra Source Catalog](https://cxc.harvard.edu/csc/)

## Replay

From this directory, run the offline calculation using Python, NumPy and SciPy.

~~~sh
python replay.py
~~~

The raw calibrated event files are publicly downloadable. To re-extract the cutouts, install Astropy and run the following command.

~~~sh
python extract.py /path/to/chandra-cache regenerated.csv.gz
~~~

The extractor verifies each input SHA-256, applies its sky WCS to the optical coordinates, and exports every level-2 event within 10 arcsec, including the spatial offsets, event time and mission pulse-invariant channel. Compare the uncompressed CSV bytes to the bundled cutout. The full inputs total about 0.87 GB. The later HRC and ACIS observations were selected from the Chandra observation catalog independently of the CSC source positions; their coverage uses all 2023–24 visits in the returned central-pointing set.

Catalog provenance and exact joins are recorded in [catalog_queries.md](catalog_queries.md). The optical positions are redistributed from Shafter, Hornoch & Zhao, *A Living Catalog of Novae in M31*, under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). The catalog credits W. Pietsch for the pre-May-2022 entries. This research uses data obtained from the Chandra Data Archive and the Chandra Source Catalog, provided by the Chandra X-ray Center. [Chandra catalog publications](https://cxc.harvard.edu/csc/publications/)
