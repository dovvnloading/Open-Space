# Two archival follow-ups: DESI oxygen-line wings and an unresolved LMC eclipse test

Two auxiliary investigations from the November 2026 campaign concerned different targets. The DESI measurement extends the [SDSS–DESI galaxy superposition](SDSS_SUPERPOSITION.md) investigation. The LMC binary was selected by a separate matched-phase depth screen. The [SMC binary timing analysis](SMC_06090_ECLIPSE_TIMING.md) uses its own observations, controls and follow-up predictions.

## The profile of the brighter oxygen line predicts its weaker companion

DESI observed the northwest galaxy (target 39628534716566336, redshift 0.841615) and the foreground-centered aperture 3.86 arcseconds away with the same repositioned fiber in separate exposures on 2021 November 11. In the northwest aperture we fitted the stronger [O III] line at 5008.239 Å rest-vacuum wavelength, integrating the model over each pixel and convolving it with the published DESI resolution matrix. A narrow core plus broader wings fits better than a single Gaussian. Freezing the shape before fitting the weaker 4960.295 Å transition tests whether the result persists at a second wavelength.

| Northwest aperture | Single Gaussian χ² | Core and wings χ² | Data retained |
| --- | ---: | ---: | --- |
| Strong line, fit shape and continuum | 133.65 | 51.49 | 53 pixels, 5 and 7 fit parameters |
| Weaker line, use strong-line shape | 60.06 | 46.84 | 52 pixels, 3 fit parameters each |

The modeled wing contains 20.4% of the stronger line's flux. All 12 fits that varied the wavelength window and continuum order still favored the broader-wing prediction on the weaker line. Twelve of 49 designated sky fibers passed the quality and coverage cuts. Their highest positive wing coefficients were 1.13 and 1.21 times their formal errors at the two wavelengths. Both lines were examined during exploration and came from one exposure. A separate instrument with a measured spectral response could distinguish an astrophysical line shape from shared errors in the instrumental response. The two spectra alone leave the gas dynamics and lens geometry unresolved. [Profiles, observational provenance and replay instructions](../results/archival-followup-controls/desi-oxygen/README.md) · [DESI DR1](https://data.desi.lbl.gov/doc/releases/dr1/)

## OGLE-LMC-ECL-13803 has inadequate late eclipse coverage

Of 5,858 catalog-eligible LMC binaries, 5,814 had OGLE-IV I-band files and 1,280 passed the two training-era coverage cuts. OGLE-LMC-ECL-13803 ranked second behind the cataloged transient eclipsing binary 16023. The prediction recorded before opening the 2022–2024 holdout was that *both* eclipse cores would remain deeper than in 2010–2014. [Training ranking, original prediction and replay instructions](../results/archival-followup-controls/lmc-13803/README.md)

| OGLE-IV interval | Accepted measurements | Nights in primary / secondary core | Depth relative to quadrature, primary / secondary |
| --- | ---: | ---: | --- |
| 2010–2014 | 643 | 27 / 18 | 0.228 / 0.170 mag |
| 2017–2019 | 150 | 22 / 14 | 0.395 / 0.347 mag |
| 2022–2024 | 4,871 | 0 / 1 | Insufficient coverage |

The catalog period is 1.9964371 days. The observing cadence missed the late minima despite thousands of observations. Historical OGLE-III and EROS-2 measurements show that eclipses occurred before the training era. Their filters, blending and sampled phases limit direct depth comparisons across surveys. Photometry covering both full minima could test whether the deepening persisted. [OGLE-IV release](https://www.astrouw.edu.pl/ogle/ogle4/OCVS/lmc/ecl/) · [EROS-2 release](https://doi.org/10.1051/0004-6361/202557321)
