# Two LYACAL pointings and the 1,135 Å feature in New Horizons Alice

**Noteworthy reproducibility finding.** The short Alice exposure in the published LYACAL survey row points **45.32°** from the mean of the ten long 2023 exposures matching the scattering-template description. The short exposure's boresight matches the survey row's printed J2000 position to **0.93 arcseconds**. The two matched template epochs, in 2007 and 2023, have mean boresights **0.10°** apart and agree with the paper's separate Galactic template coordinates. The observing label alone can retrieve the wrong field. [Published reduction, Tables 2 and 4](https://arxiv.org/html/2501.00787v1) · [Journal version](https://doi.org/10.3847/1538-3881/ada4a4) · [Header results and product IDs](template-results.json)

## Pointing identity

| Exposure set | Histograms | Mean Galactic boresight `(l, b)` | How it was identified |
| --- | ---: | --- | --- |
| Survey row's first MET `0555002780` | 1 | (318.91°, −81.75°) | Open aperture and 900 s integration. Its FITS boresight is J2000 RA 10.09531°, Dec −35.02894°. |
| 2023 long LYACAL set | 10 | (305.99°, −36.66°) | Open aperture, same science description, 3,600 s per histogram and the published approximate direction. |
| 2007 deep sky set | 6 | (306.04°, −36.76°) | Three one-hour and three two-hour open histograms with the expected direction and total exposure. |

The published Galactic positions in the survey row refer to aperture centers. The FITS positions are boresights. The 45.32° separation compares FITS boresights directly. The individual boresights extend up to about 0.59° and 0.87° from their respective 2007 and 2023 means, and slit overlap depends on orientation. The 2007 products match the published integration time and direction, with their individual identifiers supplied here for replay. The paper gives a separate template direction and exposure count, both of which agree with the long integrations. [Legacy product receipts](legacy-sources.json) · [Alice instrument geometry](https://arxiv.org/abs/0709.4279)

## Fixed-wavelength detector checks

The July 2026 analysis reported the feature near 1,135 Å and identified it as consistent with the N I triplet. This check uses the 2023 observations shared with that analysis. The June 2026 archive supplies the raw integer histograms and exact calibrated wavelength maps. There were **164** open and **168** closed exposures across **24** fields. Of the open exposures, **144** have usable dark brackets in the same field. Each of those exposures uses the same science pixels in the pooled and bracketing reductions. [Published line report](https://arxiv.org/html/2607.20942v1) · [Collection](https://doi.org/10.26007/kb7p-pq79) · [Protocol](PROTOCOL.md)

The table gives the fitted peak count rate in units of 0.001 counts per second per pixel. The quoted count errors include reused darks with their summed weights. Template rows also include the shared uncertainty of the empirical scattering template. The estimates omit other systematic errors, including uncertainty from the different sky footprints and wavelength calibrations. [Exposure ledger](observations.json) · [Dark comparison](contrast-results.json) · [Template comparison](template-results.json)

| Reduction | 1,135 Å contrast | 830 Å comparison | 1,549 Å comparison |
| --- | ---: | ---: | ---: |
| Bracketing darks | 0.561 ± 0.025 | 0.117 ± 0.023 | 0.156 ± 0.023 |
| Pooled darks | 0.566 ± 0.025 | 0.118 ± 0.023 | 0.164 ± 0.022 |
| Full empirical template with common dark | 0.507 ± 0.057 | 0.064 ± 0.050 | 0.271 ± 0.050 |
| Full empirical template with separate epoch darks | 0.574 ± 0.079 | 0.060 ± 0.072 | 0.101 ± 0.071 |
| Early template halves with common dark | 0.473 ± 0.080 | 0.039 ± 0.069 | 0.256 ± 0.068 |
| Late template halves with common dark | 0.538 ± 0.075 | 0.087 ± 0.065 | 0.292 ± 0.065 |

The 1,135 Å rate changes by less than 1% between the two dark-only reductions. Its sign also survives the tested scattering templates, both Stem row halves before template subtraction, and the wider and curved-baseline fits. The high detector edge shows a positive contrast under the wider fit, and the nominally dark 830 Å window shows a positive contrast before scattering subtraction. Those checks expose detector or scatter structure that matters for narrow-line interpretation. The median older wavelength map is about 2 Å lower at the same detector pixels near the feature. The reductions use the same Alice photons and their count errors do not determine the feature's astrophysical significance. [Template method and limits](TEMPLATE.md) · [Detector contrasts](contrast-results.json)

## Evidence and replay

The compact [detector windows](windows.npz) preserve the selected raw count and wavelength cells and their zero-based column indices. The [observation ledger](observations.json) provides the source member and digest for every frame. Source FITS archives remain with PDS. The [checksums](SHA256SUMS) bind the committed evidence.

With Python, NumPy and Astropy, run the following from this directory. Replace the placeholders with paths to an empty output location and the original inputs. The native archive directory should contain the two tarballs in [sources.json](sources.json). Place each 2007 file at its recorded relative path in a separate directory. The first command rebuilds both compact input files. Use this result directory in place of `<window-dir>` to replay the dark comparison directly from the deposit.

```sh
python prepare.py sources.json <archive-dir> <window-dir>
python check_contrasts.py <window-dir> <dark-output.json>
python check_template.py sources.json <archive-dir> legacy-sources.json <legacy-dir> <window-dir> <template-output-dir>
```
