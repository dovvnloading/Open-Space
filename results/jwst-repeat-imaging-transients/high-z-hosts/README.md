# NEXUS spectroscopic high-redshift host search

We tested the sky positions of galaxies with secure NIRSpec redshifts above 5 for compact, host-centered changes in NIRCam imaging. The eight NEXUS-Deep visits run from 2025 June 1 to 2026 August 7. This is a new host-selected search across an eight-epoch light curve. The blind pairwise subtraction elsewhere on this branch uses a different selection.

## Eight-epoch search and depth

| Search | Tested galaxy positions | Temporal contrasts | Positive one-visit injection placements per magnitude | AB at 90% recovery | AB at 50% recovery |
|---|---:|---:|---:|---:|---:|
| F200W | 120 | 1,867 | 815 | 26.11 | 26.70 |
| F444W | 122 | 1,936 | 839 | 26.14 | 26.78 |
| Same-sign F200W and F444W | 100 | 1,577 | 686 | 26.59 | 27.13 |

There are 127 distinct eligible spectroscopic positions. The counts specify where the cadence and masks allow an actual temporal test. The sample is not a continuous search of the approximately 50 arcmin² NEXUS-Deep survey. In two bands, a candidate requires the same object, test window and reference epochs, with signed scores in the same direction, magnitude ≥3 in each filter and joint score magnitude ≥7. The joint score is the sum of both single-band scores divided by the larger of √2 and the empirically measured paired-control scatter for that window length. No high-redshift host passes this test. Three of the 5,743 paired control contrasts pass, so candidate-level vetting remains essential. Joint recovery assumes equal injected flux density in both bands. All depth values use interpolation between the actual injections, and are conditional on the tested host-epochs.

### F200W and the single-band trigger

We binned each 30 mas F200W science mosaic onto the 60 mas sampling by averaging complete 2×2 pixel blocks. For independent-pixel propagation the binned errors are the quadrature sum divided by four. Blank-aperture and control distributions supply the empirical noise corrections. We obtained 854 covered measurements at the 127 unique parent positions. One position passed the single-band threshold in F200W.

MSAID **25213** (RA 268.433576°, Dec +65.254523°, secure galaxy z = 5.612) had a **negative** F200W contrast of −29.66 nJy, score −7.45, in January and March 2026 compared with its other six visits. F444W gave −1.24 and the joint score was −5.89 for the identical window. The galaxy spectrum measures the host redshift; the deficit has no independent redshift. We do not assign it a transient classification or absolute magnitude.

* The difference signal broadens as the aperture grows. With the alternative 0.48–0.66 arcsec sky annulus, the 0.12 arcsec and 0.36 arcsec apertures give deficits of 11.57 and 52.21 nJy. The ratio 0.22 differs from the median 0.846 for same-epoch stars. Subtracting the expected flux growth for a centered point source leaves a residual 3.50 times the robust scatter of 38 surrounding comparison positions. This diagnostic tests the centered-source hypothesis only. [Aperture and sky sensitivity](msaid25213/aperture_sensitivity.csv) · [Stellar and local-sky growth test](msaid25213/radial_contrast.csv) · [Two-band stamps](msaid25213/stamps.png)
* We downloaded all 39 public Level-2 NIRCam images returned by a position/epoch archive query at this host, including 15 F200W and 15 F444W frames, plus five additional filters. Only one F200W detector frame covers the March minimum in that inventory. Its unflagged 0.18 arcsec aperture gives **55.83 ± 6.90 nJy**, versus **24.14 nJy** in the QDR mosaic. That frame has only one usable registration anchor, so its archive WCS was retained. An independent February 2026 observation gives **64.55 ± 4.60 nJy**, and a clean September 2026 frame gives **65.49 ± 6.72 nJy**. Both January frames have unusable pixels in the aperture. The sole August F200W detector frame in the archive footprint has aperture DQ mask 32772 (JUMP_DET and TELEGRAPH). Differences in reduction, aperture centering and single-frame coverage limit interpretation. We have no clean, independent dither confirmation of the mosaic trough. [Detector photometry, DQ and verified file hashes](msaid25213/detector_photometry.csv)
* The nearest reported transient in the published NEXUS sample and the field's frozen TNS-broker snapshot is AT 2026rgp, 38.7 arcsec away. AT 2026abck lies 83.2 arcsec away. The nearest published variable is 11.5 arcsec away. These positional checks exclude those identifications. They do not establish that the deficit is an astrophysical event. Direct public TNS cone requests returned HTTP 403 here. [Nearest records by catalogue](msaid25213/nearest_reports.csv)

The single-band measurement remains in the output with its diagnostics. It does not survive the centered point-source and independent-image checks required for a discovery.

### September 2026 catalogue check

We separately compared the September 28 source catalogues with the mean May 28 and August 7 catalogue fluxes. We required a unique catalogue counterpart within 0.15 arcsec of the same galaxy in all three visits. The field controls provided per-visit relative scaling and a brightness-dependent residual scatter. Ten high-redshift positions could be tested in F200W, 24 in F444W and seven in both. The largest absolute single-band scores were 1.77 in F200W and 3.04 in F444W. The maximum absolute joint score was 1.95. Catalogue selection excludes some faint or shifted counterparts. The image injection depths in the table above apply only to the eight QDR visits. [September results and original input catalogue hashes](september/)

### F444W detail

Among 127 distinct eligible sky positions, 122 had adequate coverage for at least one temporal comparison, spanning z = 5.003–8.835. There were 865 measured host-epoch apertures and 1,936 usable contrasts. None exceeded the two-sided selection threshold of 7. The largest absolute score was 4.83, a negative single-epoch residual for MSAID 9800 at z = 8.149. No source receives a transient classification. Among 440 selected comparison objects, 427 contribute 7,430 F444W control contrasts. Four controls exceed threshold. Their difference images show positive and negative structure around the host or stellar centroid, illustrating why a numerical threshold alone would not establish a transient.

| Injected F444W AB magnitude | Recovered / 839 host-epochs |
|---|---:|
| 25.0 | 837 |
| 25.5 | 830 |
| 26.0 | 798 |
| 26.5 | 641 |
| 27.0 | 250 |
| 27.5 | 29 |

A linear interpolation between measured bins gives 90% and 50% recovery at AB 26.14 and 26.78. These magnitudes refer to the flux of the injected stellar template inside its 0.54 arcsec reference radius. The depth applies to an isolated positive, single-visit event centered on an eligible galaxy. The test is insensitive to sufficiently offset sources and steady emission over the whole time series. It cannot assign a volume or an all-field SN rate.

## Selection and measurements

* The public Epochs 1–8 MSA target tables define the parent sample. Keep MSAIDs whose visually inspected redshifts have confidence 2, all exceed 5, and agree within Δz/(1+z) < 0.015. Use the median reported secure redshift. Three different MSAIDs have exactly duplicated coordinates among the 130 eligible identifiers. Retain one position per pair and record the aliases. Redshift uncertainties are absent from the source table. An aperture light curve at a spectroscopic galaxy position measures combined host and potential transient flux.
* Compare with 440 lower-redshift objects with priority class ≥2 and positive catalogued F444W flux, sampled reproducibly in eight brightness quantiles. Zero-redshift controls above 0.45 µJy anchor the per-epoch registration, relative photometry and PSF. The F444W epoch fits use 11–14 stars. The native release carries an uncertainty map for each epoch.
* On the 60 mas F444W grid, measure a 0.18 arcsec radius aperture with a pixel-phase-aware fractional mask and a local median sky from 0.36–0.498 arcsec. Require complete aperture science and positive finite error coverage. Subtract the background from each aperture and correct the relative response using the stars' median per-epoch ratio. The light-curve pipeline errors omit correlated-noise and model terms. The hypothesis test separately estimates the contrast uncertainty.
* For each consecutive window of one, two or three visits, subtract the mean of at least two remaining visits with valid coverage. Subtract the median temporal difference measured at surrounding apertures 1.14 and 1.44 arcsec from the galaxy. At least 14 such apertures must have coverage throughout the identical epoch selection. The noise scale uses the larger of 1.4826 times the surrounding-aperture median absolute deviation and the propagated independent pipeline error. Inflate it, when needed, by the control-residual dispersion for the same window length and host-brightness quartile. The resulting score is an empirical selection statistic, not a Gaussian false-alarm probability.
* Build the empirical point-source template for each visit from normalized, aligned stellar stamps and truncate at 0.54 arcsec. For every eligible host-epoch and each test magnitude inject it into the actual science stamp at the host phase, remeasure the aperture including its sky, and require a positive seven-score detection at the injected epoch. The local comparison apertures are outside the injected region. The redshift-target catalogue fixes locations in advance. These tests measure recovery on host backgrounds, conditional on the real observing cadence and mask. They do not include image reduction or independent-exposure vetting efficiency.
* The branch's frozen field comparison list incorporates Zhuang et al. transient positions, TNS-based public broker records, nuclear variables and MSA target classes. The host result table gives the nearest comparison and separation. No high-redshift host requires a claim of novelty, cosmic-ray rejection or a host-transient assignment.

## Reproduction and files

Run from this directory with Python, NumPy, pandas, SciPy, Astropy and Matplotlib. Allow network access to the NEXUS release, sufficient runtime to decompress the mosaics, and space for compressed postage stamps. The cache contains the dated original target tables, compressed stamps and SHA-256 digests of every full input mosaic. The reproduction commands use separate caches for the two bands.

```sh
python acquire_nexus.py --band f444w --cache cache-f444w
python screen_hosts.py --cache cache-f444w --out f444w --comparisons ../catalogs/field_comparisons.csv
python acquire_nexus.py --band f200w --cache cache-f200w
python screen_hosts.py --cache cache-f200w --out f200w --comparisons ../catalogs/field_comparisons.csv
python join_bands.py .
python september_catalog.py cache-f200w/targets.csv cache-september NEW_EMPTY_OUTPUT_DIRECTORY
python inspect_trigger.py cache-f200w/targets.csv msaid25213
python check_morphology.py cache-f200w cache-f444w msaid25213
python radial_test.py cache-f200w msaid25213
```

The per-band result tables provide every targeted galaxy and best contrast (`host_screen.csv`), per-visit aperture photometry including missing coverage (`host_lightcurves.csv`), all tested contrasts for hosts and controls (`contrasts.csv.gz`), actual injections and binned recoveries, epoch astrometry and scaling, and the mosaic checksums. The reference epoch mask is a bit mask with bit 0 for epoch 1. The contrast field `start` is zero-based, while `epoch` and `start_epoch` are one-based. The tabulated host coordinates and redshifts in the per-band results apply throughout. `host_catalog_photometry.csv` contains sparse survey-published F200W and F444W photometry in µJy at epochs with a corresponding MSA target entry, including its original missing-value sentinels. These catalogue fluxes and the smaller-aperture forced fluxes have different units and aperture corrections.

Data and attribution: [NEXUS QDR](https://ariel.astro.illinois.edu/nexus/qdr/), [survey overview](https://arxiv.org/abs/2408.12713), [quick-release notes and coverage caveats](https://arxiv.org/html/2603.04586v3), [Year-1 transients](https://arxiv.org/abs/2609.06985), [reported z=5.017 AT 2026abck](https://arxiv.org/abs/2609.38373). AT 2026abck was reported in F200W with a nondetection in F444W, illustrating the band dependence of this sensitivity.

Acquisition verifies the deposited target-table and mosaic byte hashes before accepting data. Stamp caches bind the complete ordered target table, extraction configuration, upstream mosaic digest and compressed stamp digest; old caches without this identity, partial caches, and changed inputs are rejected. Use fresh separate per-band cache directories. The September replay verifies its six frozen catalogue inputs, regenerates and checks the target selection from the frozen target tables, and requires a new empty output directory; deposited receipts are not rewritten. Changed releases require a separately frozen campaign and updated code/receipts.
