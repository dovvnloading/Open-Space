# Independent ground-based eclipse-clock check

Two ground surveys recorded the known eclipsing variable ASAS J173806−4920.9 (TIC 96100077). We tested whether the eclipses shift in the direction and at the orbital phase predicted by the Gaia combined astrometric and spectroscopic orbit and the earlier TESS timing analysis. We fixed the eccentric waveform, phase, 272.720177-day period, and signed comparison amplitude −0.801640959 before fitting the ground magnitudes.

| Ground photometry | Usable magnitudes | Camera/field-years | Free scale | 30-day score SE¹ | Improvement for fixed TESS amplitude² | Shifted phases as good / 99 |
|---|---:|---:|---:|---:|---:|---:|
| ASAS-3, 2001–2009 | 819 | 15 | −0.689 | 0.207 | 7.060 | 1 |
| ASAS-SN, 2016–2024 | 937 | 16 | −0.479 | 0.191 | 2.410 | 11 |
| Both surveys | 1,756 | 31 | −0.574 | 0.140 | 9.470 | 3 |

¹ Cluster sandwich uncertainty on the corresponding *linearized* scales (−0.676, −0.476 and −0.565). The full pooled profile and the linearized solution differ by 0.009.

² Reduction in weighted squared magnitude residuals relative to a zero-delay clock. Positive values favor the prediction. Empirical weights make this a relative fit score. The 99 neighboring phases are correlated comparisons.

## Ground measurements and clocks

The ASAS-3 archive returned 897 aperture-1 V measurements across two fields. Grades A/B, positive aperture error below 0.2 mag, and magnitude between 5 and 20 retained 841. Sky Patrol v2 returned 1,220 rows at ASAS-SN ID 515396329465, the only master-list object within 7.2 arcsec. Of 1,071 flagged G, 1,068 passed the same error and magnitude cuts. Camera/field-years needed 30 observations and coverage of nine of twelve inner phase bins. Sky Patrol data span four camera identifiers (bG, bg, bk and bo). The upper date in the table is the last year with usable phase coverage. [ASAS-3 target photometry](https://www.astrouw.edu.pl/cgi-asas/asas_cgi_get_data?173806-4920.9,asas3) · [Sky Patrol v2](https://arxiv.org/abs/2304.03791) · [ASAS-SN variable identification](https://asas-sn.osu.edu/variables/c6f95d92-ed00-5248-b30a-ecf37de90d39)

We interpreted ASAS-3 HJD and Sky Patrol JD as UTC. To obtain BJD(TDB) we added the time-scale correction and projected the barycentric Sun displacement for ASAS-3 or Earth displacement for Sky Patrol. Geocentric observer positions suffice at the present precision. For 165 of 167 legacy ASAS-SN V timestamps within 1,200 seconds of a Sky Patrol camera-bg timestamp, the median absolute HJD−JD difference falls from 207.4 to 2.9 seconds after applying the calculated heliocentric correction. These overlapping legacy observations enter only the timestamp check. [Legacy HJD photometry](https://asas-sn.osu.edu/photometry/e6d1bc3f-883d-5527-9ddf-96aec2c38018.csv)

## Fit and sensitivity

Each camera/field-year has an independent four-harmonic light curve with a free zero point and free sine/cosine coefficients, accommodating passband, level, eclipse shape and local phase shifts. The timing signal comes from changes *within* those groups. The model evaluates each observed magnitude at its inner phase corrected by the hypothesized time delay. It profiles over the Fourier coefficients at every trial and fits a linearized period drift term per survey. Within-group photometric scatter sets the scale and residual outliers receive Huber weights trained under the zero-delay null. No search over the inner period or outer phase enters the nominal fit. The replay code uses the parent dossier's Gaia displacement and time-conversion conventions. [Timing fits](ground-fits.csv) · [Phase comparisons](ground-phase-controls.csv)

- Dropping each of 18 calendar years gives pooled scales from −0.659 to −0.510. Independently fitting linear clock drift in every channel-year gives −0.830. Changing the Fourier order from four to three or five gives pooled scales −0.901 and −0.539. ASAS-SN alone gives a negative free scale in each case, but the exact fixed amplitude worsens its five-harmonic fit and its four-harmonic fit with 2017 or 2019 omitted. [Year omissions](ground-year-controls.csv) · [Timing fits](ground-fits.csv)
- A 4,999-draw sign-flip test on the residualized, linearized null score gives one-sided conditional probabilities 0.0004 for 30-day blocks, 0.0002 for 90 days and 0.006 for whole calendar years. Each assumes independently sign-symmetric blocks at its chosen scale. These diagnostics condition on the selected star, prior Gaia/TESS template and fitted nuisance model. A population-wide false-alarm rate would require rerunning the selection. [Cluster controls](cluster-controls.csv)
- The fixed Gaia-phase improvement of 9.470 exceeds the best annual (0.006), half-year (−5.594) and synodic-month (−10.864) improvements, each given 100 phases at the TESS amplitude. The best pooled Gaia-period phase is 0.98, about 5.5 days from nominal. Varying the Gaia period and periastron time by their reported marginal ±1σ errors in a 3×3 grid gives improvements from 7.71 to 9.60. This is a point-estimate stress test and omits full orbit covariance. [Calendar controls](calendar-controls.csv) · [Orbit parameter grid](orbit-uncertainty-grid.csv)
- The two halves of the inner cycle give linearized scales −0.612 ± 0.187 and −0.512 ± 0.218 when estimated together. The difference is −0.100 ± 0.292. The ASAS-SN-only halves are less consistent (−0.802 vs −0.113, difference −0.689 ± 0.338). Injecting zero and ±0.801641 into noiseless templates on the exact observing calendar recovers all three to within 0.000002. That exercise checks implementation sign and response. [Eclipse controls](eclipse-controls.csv)

The ground data independently favor a modulation with the predicted phase and sign. Survey-dependent waveform variations affect the inferred amplitude. The archival apertures include neighbors, so the test supplies timing evidence while the TESS pixel localization constrains the source position. Component-resolved spectroscopy would test the physical configuration.

## Reproduce

From the parent dossier, run the following with NumPy, SciPy, pandas and Astropy. The committed observations include the columns used in analysis from every source row, including rejected rows. Current packaged outputs were generated on 2026-11-14 UTC.

```sh
python ground-check/ground_clock.py
python ground-check/calendar_controls.py
python ground-check/eclipse_controls.py
python ground-check/uncertainty_controls.py
```

A rerun of `eclipse_controls.py` can differ from the committed `eclipse-controls.csv` in one column. That column holds remainders between one hundred-millionth and one millionth, which are round-off of the linear solve and depend on the installed NumPy. On 16 November 2026 a rerun on Windows 11 with Python 3.14.0, NumPy 2.5.4 and SciPy 1.18.1 changed them in their fifth significant figure. Every other value was identical. No number quoted in this folder changes, and the committed table is unchanged.

To regenerate the compact input from the raw products, use the URLs above for `asascgi.txt` and `asphoto_csv`. Obtain `patrollc.parquet` with the Sky Patrol client cone search at RA 264.526816°, Dec −49.349215°, radius 0.002°, catalog master_list, with download enabled and columns asas_sn_id, ra_deg, dec_deg and catalog_sources. The official client retrieves block 0 for that cone, which supplied the parquet named below. Run `python ground-check/prepare.py /path/to/downloads` with pyarrow installed. Live archive responses may gain observations. [Official Sky Patrol client](https://github.com/asas-sn/skypatrol/blob/master/pyasassn/client.py)

| Original product | Source rows | Raw SHA-256 on retrieval |
|---|---:|---|
| ASAS-3 CGI light curve | 897 | `94bd75b71e4867c92201f3343aceb9f170645c8e17f6a03abdf153a4ff74d0e9` |
| Sky Patrol v2 parquet | 1,220 | `875571318b96daedc17a0d9aaba4f9e10a83033c45a0cd392d68d9e375bcb178` |
| Legacy ASAS-SN V CSV | 167 | `2a72e19f7ebd6c67946724892295ced93a8894f38a9591cf28dfab656c51da32` |
