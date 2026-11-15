# Toruń recorded the predicted 112.76-day modulation a decade earlier

The Hitachi 17.01 km/s phase and period, fitted to its MJD 58000–61000 data, predict a positive 0.339 ± 0.095 Jy cosine amplitude in the 16.92 km/s feature monitored by a different 32 m telescope in Toruń. The archived 44 Toruń epochs run from 2009 to 2013. The first and second chronological halves independently give positive amplitudes of 0.360 and 0.331 Jy. The weaker 16.13 km/s channel points in the same direction. Both velocity channels came from the same 44 Toruń spectra and may share calibration errors.

| Toruń feature (km/s) | Hitachi counterpart (km/s) | Fixed period (days) | Fixed maximum (MJD) | Fitted cosine amplitude (Jy) | One-sided block-bootstrap p, range |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 16.92 | 17.01 | 112.75746 | 59067.65 | +0.339 ± 0.095 | 0.00050–0.00140 |
| 16.13 | 16.19 | 112.76950 | 59062.28 | +0.170 ± 0.085 | 0.0244–0.0300 |

The Toruń signed-amplitude test used the Hitachi period and maximum supplied before testing the Toruń data. A constant, a linear time trend and a cosine are the only fitted terms. Quoted uncertainties are ordinary least-squares standard errors. Under an exponential-covariance noise model with a white-noise component, the strongest line gives +0.336 ± 0.109 Jy and one-sided conditional Monte Carlo p ≈ 0.001. For the weaker line the same calculation gives +0.163 ± 0.089 Jy and p = 0.036. The second channel corroborates the direction, with much less precision. The two p-values should not be multiplied.

After the fixed-phase test, a descriptive fit allowing a sine component shifts the 16.92 peak +7.1 days from the extrapolated Hitachi maximum and the 16.13 peak −5.7 days. Fluxes are plot-derived. The data have neither the original epoch-level errors nor all information about irregular calibration noise. The period estimates in this test came entirely from Hitachi. No Toruń period scan contributed to the test.

## Recovering the earlier measurements

Szymczak et al. monitored 139 objects and described the flux variations of this particular weak source as dominated by noise. Their Table 2 lists 44 observations over 1,315 days with mean peak fluxes of 2.9 and 3.5 Jy. The independent extraction recovers 44 marker centers at each velocity, MJD 55016.8–56331.0 and means of 2.90 and 3.51 Jy, respectively. The quoted dates should be treated as approximate to one day. The PDF's eighth page, numbered 37, contains the G28.011−0.426 panel in **Figure 3, continued**. The vector marks are directly encoded in the PDF. We use the seven labeled date ticks for an affine horizontal calibration, and the 1 and 5 Jy ticks for the logarithmic vertical scale. Isolated marker outlines and the legend lie outside the color-and-panel selection. [Article and publisher supplement](https://repozytorium.umk.pl/items/bb0ff775-b503-4ebb-924f-76059d513414) · [Original figure PDF](https://repozytorium.umk.pl/bitstreams/fad49269-d41d-4bab-8601-dfd7e6581847/download) · [Toruń Table 2](https://academic.oup.com/view-large/102995749)

European VLBI Network mapping of G28.011−00.426 independently places the 16.1 km/s reference feature at J2000 **18h 43m 57.96865s, −04° 34′ 24.0938″**. The 15.8–18.0 km/s blue emission is distributed over approximately 100 × 100 milliarcseconds. The interferometric position and velocity span support the assignment to this maser source. A single-dish light curve does not isolate one specific VLBI spot. [Bartkiewicz et al. 2016, Table 2 and Appendix A](https://arxiv.org/abs/1601.03197)

A 2020 workshop presentation from Hitachi already reported a group of 20–30 new periodic sources, credited to work in preparation. A 2024 report from the Hitachi team states that it had identified roughly twenty additional periodic masers after the twenty-nine already known. Its results section gives the group total and omits the individual designations. Searches by G28.011−0.426, its rounded designation and the period did not locate a publication assigning this particular clock. The broader Hitachi report leaves open whether their team independently recognized this source earlier. The evidence here supports a source-specific ephemeris and external validation. Publication priority remains unresolved. [Yonekura 2020, monitoring summary](https://oister.kwasan.kyoto-u.ac.jp/wordpress/wp-content/uploads/2020/11/yonekura_ws11.pdf) · [Yonekura et al. 2024, section 3](https://doi.org/10.1017/S174392132300251X)

## Reproduce

Python 3 with PyMuPDF, NumPy, pandas and SciPy is sufficient. From this directory, run the following commands. The extractor downloads the institutional repository copy and verifies SHA-256 `3010a9ca4c8965dae33c28324ccd6c5bb7617d151f11d11896f43658ed0573ad`. Supply `--pdf /path/to/stx2693_supp.pdf` for offline extraction. The inference step can also read the committed CSV directly.

```bash
python extract.py
python validate.py
```

The block test resamples circular contiguous groups of two, three, four or six visits from the time-ordered null residuals, retaining actual observation times and refitting the cosine on each realization. The exponential-covariance test fits white and correlated noise to detrended fluxes under the null, uses the actual observing cadence, and draws 20,000 realizations at those estimated nuisance parameters. Its quoted standard error and p-value are conditional on that covariance fit. The times and fluxes come from the same observing run throughout, with distinct epochs as the units of replication. [Extracted photometry](torun.csv) · [Numerical results](fixed_phase.json)
