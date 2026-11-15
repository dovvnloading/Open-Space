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

## The optical spectrum and X-ray energy bands

The three 600 s Liverpool Telescope SPRAT exposures from 2017 June 13 resolve Balmer emission from **M31N 2017-06c**. Williams & Darnley reported these same observations, identified the Fe II class and measured an Hα width of about 1,100 km/s. We retrieved the level-2 spectra from programme JL17A03 and measured the apparent widths per exposure. The following values include the instrumental line spread. [Observers' report](https://www.astronomerstelegram.org/?read=10487) · [SPRAT data products](https://telescope.livjm.ac.uk/TelInst/Inst/SPRAT/) · [FITS manifest](sprat-observations.csv)

| Exposure start (UTC) | Hα centroid (Å) | Hα FWHM (Å) | Hα FWHM (km/s) | Hβ FWHM (km/s) | Nearby sky [O I] 6300 FWHM (Å) |
|---|---:|---:|---:|---:|---:|
| 03:48:53 | 6558.52 | 29.36 | 1341 | 1354 | 16.33 |
| 03:59:02 | 6558.65 | 30.14 | 1377 | 1342 | 15.52 |
| 04:09:10 | 6557.95 | 29.24 | 1336 | 1519 | 16.79 |

![Three repeated SPRAT spectra of M31N 2017-06c and the sampled Hα line profiles](sprat-profiles.png)

The plotted spectra use a common continuum interval for normalization and the inset scales each continuum-subtracted Hα profile to its own peak. [Plot code](plot_sprat.py)

Each fit uses an unweighted Gaussian plus a sloped continuum over ±80 Å about the laboratory wavelength, on the archive's 9.2 Å grid. Moving the edges to ±60 or ±100 Å and independently extracting nine slit rows around the spatial peak with local sky subtraction gives Hα widths of 1,279–1,421 km/s across all fits. This range measures method and exposure sensitivity. The three share the same wavelength calibration. The nearby sky line is narrower in every exposure. Gaussian quadrature subtraction using its width as an approximate instrumental reference yields Hα widths near 1,090–1,180 km/s for the default fits, consistent in scale with the original report. Slit illumination differs for the sky and a point source. A coarse grid, blending, sky subtraction and the line-spread function limit the precision. The archive supplies relative detector counts for this extraction. [All fits](sprat-lines.csv) · [Sampled spectra and slit-sky profiles](sprat-spectra.csv)

For **M31N 2021-06b**, MASTER reported an unfiltered magnitude of 15.9 at this sightline in August 2018 and a limiting magnitude of 19.5 on 2021 June 8. In June 2021 Swift/UVOT measured uvw1 = 17.82 ± 0.08 AB on June 13, fading below the uvw1 detection limit of 20.0 AB by July 2. The earlier MASTER detection complicates identification of the eruptive object in this crowded field. The 2021 data support a UV-bright transient. Recurrence, blending and the optical classification remain open. A 30-arcsec Liverpool archive search returned 2019 spectra of two nearby, differently positioned targets, with no level-2 SPRAT spectrum of this eruption. [MASTER](https://www.astronomerstelegram.org/?read=14703) · [UVOT June](https://www.astronomerstelegram.org/?read=14716) · [UVOT July](https://www.astronomerstelegram.org/?read=14762)

Two Swift/XRT photon-counting visits bracket the 2021 October 18 and 20 Chandra detections. A six-pixel aperture (14.1 arcsec radius) at each optical position contained at most one 0.3–10 keV event. These short exposures yield upper limits. The single event at 2021-06b had PI = 41 (nominally 0.41 keV), while the October 13 annulus predicts 0.87 soft-band background events there. It carries no individual source attribution.

| Date (2021 UTC) | Sightline | Mean aperture exposure (s) | Aperture events / scaled annular background (0.3–10 keV) | One-sided 95% upper limit (counts/ks) |
|---|---|---:|---:|---:|
| Oct 13 | 2017-06c | 1515 | 0 / 2.85 | 1.98 |
| Oct 13 | 2021-06b | 1528 | 1 / 2.02 | 3.10 |
| Oct 27 | 2017-06c | 125 | 0 / 0.19 | 23.95 |
| Oct 27 | 2021-06b | 122 | 0 / 0.16 | 24.62 |

The Swift selection uses the two HEASARC `swiftmastr` pointings with positive PC exposure, start MJD 59492–59518, RA 10.30°–11.10° and Dec 40.96°–41.55°. Both cover both targets in their calibrated exposure maps. The cleaned event files' PC exposures are 1,533 and 131 s. We scale annular counts by the ratio of summed map exposures in the source and background regions (12–24 pixel background annulus). We keep grades 0–12 and PI 30–1000. The exposure in the table is the mean map exposure in the source aperture. Limits use the one-sided 95% Poisson bound on the total aperture counts and omit background subtraction to give conservative source limits. All rates refer to the aperture and have no PSF or energy-flux correction. Diffuse emission and neighbouring sources enter the large Swift apertures. Their response and coverage differ from the October 18–20 HRC observations. [Event and exposure-map manifest](swift-observations.csv) · [Swift/XRT documentation](https://www.swift.ac.uk/user_objects/docs.php)

In the four August 2024 ACIS pointings (117.803 ks total), both 2-arcsec optical apertures contained **zero events between 0.3 and 2 keV**, giving a marginal 95% upper limit of 0.0254 aperture counts/ks in this band per target. At 2–7 keV the apertures contained 1 and 2 events against expected backgrounds of 1.17 and 0.92, respectively. The ACIS bands use the PI-channel midpoint, (PI − 0.5) × 0.0146 keV. The background region is the original 4–8 arcsec annulus. ACIS rates here divide by the nominal observation exposures as in the original replay. These late counts and the shallow Swift visits leave the spectrum, temperature and luminosity of the 2021 HRC sources undetermined. [Band-by-band counts and upper limits](xray-bands.csv) · [ACIS PI convention](https://cxc.cfa.harvard.edu/ciao/dictionary/pi.html)

Run `python sprat.py` for the line-fit CSV and `python xray_bands.py --swift-cache /path/to/cache` for the X-ray CSV. The latter fetches four small public Swift products by checksum and requires Astropy. Omitting the cache option replays the ACIS rows offline. To regenerate the SPRAT data, select night of **20170612** (local observing night), instrument **SPRAT**, level **2**, programme **JL17A03** in the Liverpool archive. Download the three products in the manifest and run `python sprat.py --extract /path/to/fits --output regenerated.csv`. The extractor checks SHA-256 of the decompressed FITS. The output contains three detector-count spectra from SPEC_SS, a second extraction from LSS_NONSS for each, and the median slit sky for the wavelength and resolution check. The archive records an arc-fit status flag of 3 and the same comparison arc for all three exposures. The sky-line checks apply at the two tabulated lines. [Liverpool archive](https://telescope.livjm.ac.uk/cgi-bin/lt_search) · [Replay code](sprat.py) · [X-ray extraction](xray_bands.py)

The Liverpool Telescope on La Palma is operated by Liverpool John Moores University at the Instituto de Astrofísica de Canarias' Observatorio del Roque de los Muchachos, with support from the UK Science and Technology Facilities Council. We credit the JL17A03 observers and the UK Swift Science Data Centre, Swift team and HEASARC archive for the public data. [Telescope and instrument references](https://telescope.livjm.ac.uk/Ack/)

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
