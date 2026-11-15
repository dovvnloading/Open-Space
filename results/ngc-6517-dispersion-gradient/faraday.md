## A second radio observable tests the electron-column gradient

Faraday rotation measures the electron density weighted by the line-of-sight magnetic field. Liu, Wang & Pan published raw and ionosphere-corrected rotation measures (RMs) for pulsars A–D. They identified the high overall RM of NGC 6517 near the Galactic plane. We compare their four measurements with the independent timing observable DM from Yin et al. The observations come from the same telescope. [Liu et al. (2025), Table 2](https://doi.org/10.1088/1674-4527/add566) · [Searchable author postprint](https://chinarxiv.org/items/chinaxiv-202506.00083) · [Yin et al. (2024), Table 1](https://arxiv.org/abs/2405.18228)

| Pair (first minus second) | Separation | ΔDM (pc cm⁻³) | ΔRM, corrected (rad m⁻²) |
|---|---:|---:|---:|
| B − D | 71.62″ | +7.865 ± 0.004 | −3 ± 2.24 |
| C − B | 2.57″ | −0.046 ± 0.0036 | +26 ± 2.24 |

Under a common background and a single extra electron column with uniform line-of-sight field, RM should be linear in DM. A weighted fit to the four corrected RMs gives χ² = 140.99 for 2 degrees of freedom using the published errors as independent Gaussian standard deviations. Using the raw RMs gives 130.47. Omitting C gives χ² = 0.0018 for 1 degree of freedom. An independent error floor of 8.51 rad m⁻² added in quadrature to each corrected measurement raises the fit probability to 0.05. The two pair contrasts share B and are descriptive. The fit uses the four original rows. The formal DM uncertainties are negligible at the fitted slope.

C was observed on 2022-12-28, three days before A and B. D was observed on 2023-09-29. C's reported integrated linear polarization is 5 ± 8 percent. The observations provide one published RM per source and lack an epoch-to-epoch covariance model. Profile-dependent RM biases, changes between epochs, additional magnetized plasma, and field variation remain possible. The near-zero B–D RM difference limits the Faraday change across the strongest sampled DM contrast. Interpreting its ratio as a magnetic field would require attributing both differences to the same material with a common foreground. A repeat polarization observation of C together with B would directly test the 26 rad m⁻² contrast.

Run `python faraday.py` alongside `fast.csv` to regenerate `faraday.json` from the four Table 2 rows in `faraday.csv`. The raw and corrected RMs use the same polarization observations.
