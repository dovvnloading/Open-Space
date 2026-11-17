# DESI DR1 radial velocities and spectroscopic validation

**Two Gaia-matched stars have high line-of-sight speeds confirmed by absorption lines, each on one observing night.** Gaia DR3 1717410430709663616 has an uncorrected DESI RV of +503.49 ± 4.20 km/s. The published backup-program correction gives **+501.74 km/s** and a solar-motion-corrected line-of-sight component of approximately **+658 km/s**. Gaia DR3 4444095982221534720 has +473.53 ± 5.88 km/s and a line-of-sight component of +589 km/s. One night cannot constrain the centre-of-mass speed of an unresolved binary. The main README gives Galactic-frame posteriors, catalogue matches and origin tests.

## What survived the radial search

The Astro Data Lab table `desi_dr1.mws` contains 4,640,216 coadds satisfying Redrock class STAR, RVSpecFit warning zero and 0 < formal RV error < 20 km/s. The screen admits both signs. For this branch, $v_{\rm GSR,los}=v_r+(11.1\cos l+245.6\sin l)\cos b+7.25\sin b$, with velocities in km/s. This is one projected component. The catalogue query used raw barycentric RVs, including the backup program, and applied its published per-target correction during follow-up. A threshold in this component misses transverse and mixed-direction ejecta.

| Coadds | Difference between RVSpecFit and Redrock above 60 km/s | Fraction and 95% exact binomial interval |
|---|---:|---:|
| 57 with absolute projected speed above 580 km/s | 50 | 87.7% [76.3%, 94.9%] |
| 20 of these above 1,000 km/s | 20 | 100% [83.2%, 100%] |
| 4,641 catalogue control rows, random_id < 0.1 | 1 | 0.022% [0.00055%, 0.120%] |
| 57 moderate-speed controls matched by program and red-arm S/N | 0 | 0% [0%, 6.27%] |

Redrock and RVSpecFit fit the same exposure. The table measures a gross RV-fit discrepancy, using `abs(vrad − 299792.458 × rr_z) > 60` km/s. It does not estimate the false-discovery fraction of a multi-epoch unbound-star catalogue. Controls include survey matching when a matching subset is available. The finite intervals assume independent rows. Seven of the 57 tail rows have agreeing velocities. Four have a Gaia DR3 source ID, and three carry DESI's missing-ID sentinel. Four sky-matched sources were inspected for follow-up.

| Gaia DR3 | DESI TARGETID | Heliocentric RV, km/s | Uncorrected GSR LOS, km/s | Assessment |
|---|---|---:|---:|---|
| 1717410430709663616 | 2305843022006328308 | +503.49 ± 4.20 | +659.27 | One exposure, 2022 May 11. Hβ and Hα centroids are +491 ± 21 and +513 ± 12 km/s. The three Ca II lines give +503 ± 14, +483 ± 12 and +529 ± 18. The backup wavelength correction is +1.75466 km/s to subtract. |
| 4444095982221534720 | 39628017630184178 | +473.53 ± 5.88 | +589.15 | One exposure, 2022 April 13. Hβ, Hα and Ca II 8542 give +491 ± 25, +495 ± 17 and +443 ± 36 km/s. The negative parallax and disagreement between atmospheric-parameter pipelines reduce distance precision. |
| 1545579445011222272 | 39633208156884458 | +584.88 ± 4.95 | +647.84 | One exposure. RUWE 2.38 and negative parallax. Projected 7.63 arcmin from M106, which has a systemic RV around +507 km/s in SIMBAD. Milky Way membership remains unresolved. |
| 1525802834304705024 | 39633113709543745 | +567.84 ± 4.27 | +638.32 | One exposure. Gaia does not report proper motion or parallax. Projected 7.57 arcmin from M63, systemic RV around +500 km/s. Milky Way membership remains unresolved. |

In the extension to absolute projected speeds above 500 km/s, nine Gaia-matched coadds pass the same RV cuts and agree across the two fitters. Two have multi-night RVs. Gaia 1070549357629086976 has eight exposures over two nights and falls 3.41 arcmin from the centre of M82. Gaia 4411527554453783040 has two exposures separated by one day and lacks Gaia proper motion. Coordinates, surrounding galaxies, epoch data and input identifiers are in `lower/`. Neither enters the Galactic 6D candidate set.

The independent full 6D *nominal* parallax check found 1,880,343 joined coadds with positive raw parallax above 5σ, RUWE < 1.4, at least nine visibility periods, S/N ≥ 5 in at least one arm and the RV cuts above. We evaluated the 150 rows with nominal heliocentric 3D speed above 700 km/s. All four with nominal Galactic-frame speed above 1,000 km/s have grossly discrepant Redrock RVs. This inverse-uncorrected-parallax prefilter supplements the separate zero-point-corrected and spectrophotometric inference in the parent analysis. Exact inputs and cuts are in `control/`.

## Transverse candidates and exposure tests

Two transverse candidates in the parent analysis have multiple nights of DESI data. We checked true single-exposure RVSpecFit measurements, exposure IDs, dates, fibre status and fit-warning bits. For the constant-RV test, independent formal errors determine the within-night weighted mean. We added a 1 km/s calibration floor per night in quadrature and tested the means against one common velocity.

| Gaia DR3 | Exposures / nights | Baseline | Constant-RV test | Independent distance-model checks |
|---|---:|---:|---|---|
| 2752545725694798336 | 11 / 4 | 7.03 d | χ²(3) = 3.61, p = 0.306 | SPDist's quality and training flags pass. Its dereddened BP–RP = 0.596, T = 6,098 K, log g = 4.11, median absolute G = 4.12 and distance = 5.72 (+0.72/−0.52) kpc. SpecDis gives 8.10 ± 1.84 kpc. |
| 6340918677332110720 | 3 / 2 | 30.91 d | χ²(1) = 0.019, p = 0.889 | SPDist reports disagreement in the atmospheric parameters and a training-domain failure for both temperature–metallicity and metallicity–colour. Its 3.35 kpc median and 2.10–3.92 kpc interval must be treated as an extrapolation. SpecDis gives 6.94 ± 1.41 kpc. |

The cadence constrains detectable changes during those observations. Constant velocity across these baselines leaves longer periods, low amplitudes and low-inclination binaries viable. The coadd line plots in `tangential/` put the measured velocity on individual unmasked pixels. The strong Ca II 8542 line gives −96 ± 5 and −127 ± 8 km/s, respectively. The SPDist quantiles tabulated here are the formal published distribution. The authors recommend adding a 0.167 mag intrinsic scatter in quadrature when propagating distance errors. SPDist and SpecDis use shared DESI data and may share training information.

Two additional candidates from the moderate transverse-speed branch have clean calcium triplets. Gaia 4435527217525195264 has two observations on DESI observing night 20220606 (MJD 59737) within one night. Its three Ca II centroids are −68 ± 4, −70 ± 5 and −56 ± 11 km/s. Gaia 4550108286326235520 has one observation on MJD 59732. Its calcium lines support the coadd RV of −84.9 ± 1.9 km/s. Both are cool stars with weak Balmer lines. Their calcium diagnostics and exposure tables are in `moderate/`. The DESI observing dates in `epochs.csv` use the local observing-night convention, and `mjd` specifies the UTC exposure time.

## Inputs and replay

Use Python with astropy, numpy, pandas, scipy and matplotlib. Start in this directory. The checked-in ADQL scripts preserve all cuts and the raw input receipts preserve release identifiers, byte lengths and SHA-256 digests. `input_hashes.json` covers checked-in input products and selected response derivatives. Database ordering can change on a repeated unordered query. Compare full ID-keyed rows when replayed bytes differ.

```bash
python control/requery.py replay/control
python control/requery.py replay/lower lower/lower_all.adql
python control/requery.py replay/los los/losall.adql los/xpavail.adql
python control/control_rates.py
python control/parallax_screen.py
python mws_fetch.py ../desi_screen.csv replay/coadds
python los/select_leads.py
python mws_epochs.py los/ids.txt los/leads.csv los
python mws_epochs.py nearby_ids.txt queries/coadds.csv nearby
python mws_epochs.py lower/ids.txt lower/lower_all.csv lower
python mws_epochs.py tangential/ids.txt queries/coadds.csv tangential
python mws_epochs.py moderate/ids.txt ../moderate/desi/queries/coadds.csv moderate
python backup_correction.py
python distance_values.py /path/to/iron-yr1-v2.1.fits /path/to/spdist.fits
python tangential/quality.py /path/to/spdist.fits
python los/other_rvs.py
python los/dr19_counterparts.py
python lower/galaxy_check.py
for id in 1717410430709663616 4444095982221534720 2752545725694798336 6340918677332110720 4435527217525195264 4550108286326235520
do
    python los/inspect_spectrum.py /path/to/coadd_cache "$id"
done
```

The 1,020 targets in the parent's SpecDis screen return 1,055 MWS coadds. `queries/coadds.csv` preserves these rows across surveys and programs. `queries/batch_receipts.json` records the hashes of the original 16 TAP responses. `epoch_products.json` in each subdirectory names each upstream RVSpecFit FITS file and gives its SHA-256, verified against DESI's manifest. The diagnostic plots use flux, inverse variance and data masks from DESI coadded spectra. Each plot has a parent FITS digest and an extracted local pixel table. Five separate sloping-continuum Gaussian absorption fits on vacuum barycentric wavelengths supply an independent line-identification and centroid check. Model assumptions and calibration errors limit the formal centroid uncertainties.

SDSS DR18 Legacy optical SpecObjAll and LAMOST DR12 LRS returned zero coordinate matches for the four radial outliers. Unfiltered SDSS DR19 MWM BOSS/APOGEE allStar queries returned zero Gaia-ID matches for these four, the two transverse candidates, the two additional moderate-speed candidates and Gaia 2294414169958935296. LAMOST DR12 MRS also returned zero Gaia-ID matches for the four radial outliers. None of the four has a Gaia DR3 RVS velocity. Gaia 1717410430709663616 does have continuous XP spectra. A query of the published Verberne et al. extended XP-RV catalogue by exact ID returned zero rows. `los/` contains the query receipts.

DESI DR1 is released under CC BY 4.0. DESI data, survey acknowledgments, wavelength-system and VAC definitions appear in the official release documentation. This analysis used services of the NSF NOIRLab Astro Data Lab. The backup calibration follows Appendix C of Koposov et al. (2026).

- [DESI MWS DR1 stellar catalogue and exposures](https://data.desi.lbl.gov/doc/releases/dr1/vac/mws/)
- [DESI SpecDis catalogue](https://data.desi.lbl.gov/doc/releases/dr1/vac/mws-specdis/) and [SPDist catalogue](https://data.desi.lbl.gov/doc/releases/dr1/vac/mws-spdist/)
- [DESI DR1 data licence and acknowledgments](https://data.desi.lbl.gov/doc/acknowledgments/)
- [Astro Data Lab DR1 access](https://datalab.noirlab.edu/data/desi)
- [Koposov et al., DESI DR1 Stellar Catalogue](https://arxiv.org/abs/2505.14787) and [target-level backup RV correction](https://doi.org/10.5281/zenodo.15469272)
- [Verberne et al., radial velocities from Gaia BP/RP spectra](https://doi.org/10.1051/0004-6361/202348406)
