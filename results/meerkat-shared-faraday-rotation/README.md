# A shared Faraday rotation signal in the MeerKAT Thousand Pulsar Array

PSR J1057−5226 has a published rise in rotation measure (RM) of **+0.126 ± 0.027 rad m⁻² yr⁻¹**. An independent fit to the 59 released measurements reproduces this to the printed precision. Other pulsars observed by MeerKAT within 30° and 2.16 hours predict essentially the entire rise on the 53 visits with adequate coverage. The target rises by **+0.133 ± 0.036** on those visits. Subtracting the offset estimated exclusively from the other pulsars gives **+0.0001 ± 0.0269 rad m⁻² yr⁻¹**. Splitting the entire sample into disjoint comparison and target pulsars gives **−0.013 ± 0.032** for this target. The original study interprets its population of gradients in the context of Galactic magnetized plasma. These observations show how much of the apparent J1057−5226 gradient is shared by nearby sightlines.

The comparison spans **304 pulsars and 11,717 target visits**. Subtracting the predicted offset reduces the scatter of 232 pulsars, including 226 after each series is independently detrended. Fifteen of the 35 slopes originally above three formal standard errors in the matched sample fall below that threshold, while ten others cross upward. All counts use identical visits within each pair of fits. Changes in formal threshold membership are descriptive. Shared observing sessions correlate the targets, and the differential model also removes a genuine sky signal if it extends across the chosen comparison field.

| Pulsar | Visits with peers | Full published-data fit | Same visits before subtraction | After subtraction |
|---|---:|---:|---:|---:|
| J1057−5226 | 53 | +0.126 ± 0.027 | +0.133 ± 0.036 | +0.0001 ± 0.0269 |
| J0942−5657 | 51 | +0.216 ± 0.048 | +0.235 ± 0.056 | +0.055 ± 0.043 |
| J1047−6709 | 52 | +0.227 ± 0.063 | +0.227 ± 0.063 | +0.078 ± 0.042 |
| J0909−7212 | 53 | +0.302 ± 0.087 | +0.302 ± 0.087 | +0.112 ± 0.068 |
| J1048−5832 | 51 | +0.923 ± 0.049 | +0.923 ± 0.052 | +0.806 ± 0.036 |
| J1114−6100 | 51 | −2.177 ± 0.103 | −2.177 ± 0.103 | −2.302 ± 0.098 |

Slopes are in rad m⁻² yr⁻¹. The table illustrates both sensitive and persistent gradients. [All 304 fits](slopes.csv) · [J1057−5226 visit measurements and predictions](J1057-5226.csv)

## How the offset was measured

The input consists of the 20,365 rows in all 597 RM tables in the 2019–2023 MeerTime release. Subtract the published ionospheric prediction from each observation and add its uncertainty in quadrature with the RM measurement uncertainty. For each pulsar, take the median corrected RM over all its visits with combined uncertainty below 2 rad m⁻². Use its departures from that median as comparisons for other pulsars. An eligible comparison has uncertainty below 1.8 rad m⁻², the same UTC MJD day, a separation below 30°, and a time difference below 2.16 hours. Average duplicate contributions from each pulsar by taking their median, then take the median across at least five other pulsars. Each target has at least 18 valid input rows, at least 12 satisfying the preliminary precision cut, and at least 18 final paired visits spanning 800 days. Formally the angle and time constraints describe two different coordinates, but the telescope schedule correlates them.

A least-squares fit weights visits by the inverse of their reported combined variance. As in the published fits, the covariance is scaled by the reduced chi-square of the residuals. The quoted errors condition on the estimated comparison series. Their interpretation does not include correlated errors in the ionosphere maps, calibrator time series, or individual pulsar noise. The full-data fits use every finite row of the eligible target, before matching visits.

## Falsification and reach

- **Disjoint pulsars.** A SHA-256 split by source name assigns each pulsar exclusively to the comparison set or the target set. This leaves 166 independently assigned targets and 6,351 target visits. The predicted offset reduces scatter for 110 targets and detrended scatter for 107.
- **Independent daily measurements.** Use the same disjoint pulsar split and a separate baseline for each pulsar on each visit, the median of at least five measurements 6–180 days from that visit. On 124 dates with at least 12 distinct pulsars per half, independently measured daily median offsets correlate at **r = 0.730**. On October 8, 2022, the two halves give +0.875 and +0.580 rad m⁻² across 66 and 93 pulsars, respectively. Offsets refer to the residual after the published ionospheric correction. [All 124 dates](dates.csv)
- **Temporal agreement.** After subtracting separate linear fits, take the mean cross-product of the target residual and the peer prediction, both divided by the target residual's standard deviation. Weight each target equally. The statistic is +0.177, and +0.159 in the disjoint test. In a diagnostic reference distribution which flips the sign of each 60-day interval together across all targets, none of 9,999 random draws reach either observed value. The intervals are anchored at MJD 58000. This reference assumes approximate block-sign exchangeability and the window choices followed inspection of the archive. [Sensitivity and temporal checks](checks.json)
- **Choices of comparison set.** The J1057−5226 differential slope varies from −0.014 to +0.0001 across the 20°/1.2 h and 40°/3 h windows, the primary selection, and the ten-peer cut. The strong negative J1114−6100 gradient survives the primary subtraction. The observations come from one telescope. Incomplete ionospheric correction, an instrumental effect, and genuine spatially coherent foreground rotation can all contribute.
- **Scope of prior art.** Keith et al. measured the individual gradients and reported the J1057−5226 value in Table A2. Their reduction used ionFR and JPL total-electron-content maps. Searches for the target, the survey, differential RM and correlations between nearby pulsars found the survey report and general ionospheric calibration work. The added measurement is the contemporaneous peer prediction and its sample-wide performance. No cause or source-specific magnetic-field change follows uniquely from this analysis. [Publication](https://doi.org/10.1093/mnras/stae937) · [ionFR method](https://doi.org/10.1051/0004-6361/201220729)

## Reproduce the result

Run `python replay.py` from this directory with Python 3 and NumPy. It writes the three result CSVs and `checks.json`. The calculation uses no network. To reconstruct the input CSV, download `tpa_data.zip` from the DOI below and run `python replay.py --extract /path/to/tpa_data.zip`. This verifies the publisher's MD5 and reconstructs positions from the supplied timing models. Ecliptic positions use the J2000 obliquity. The bundled CSV preserves the original RM table columns and adds the pulsar names and directions rounded to five decimal places. [MeerTime data release, Keith et al., CC BY 4.0](https://doi.org/10.5281/zenodo.8430591)

Original ZIP SHA-256 `3c4d3993983299337d9653b728eec13c2354b056aae6deda8434b24e8d7092ef`. Transcribed CSV SHA-256 `0536aac6f163dfba2191f8728a073a2fb9b38c03ce33425376091172cbd490ab`.
