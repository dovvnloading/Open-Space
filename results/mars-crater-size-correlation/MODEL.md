# Comparison with a published Martian meteoroid fragmentation simulation

The Separate Fragments Model generates local similarity between neighboring crater diameters. Its archived run 198 produces a smaller effect than the HiRISE catalog under the same nearest-neighbor statistic. At the 1 m diameter cut, the observed reduction in mean rank gap is **3.00 percentile points**. The matched simulation gives **1.70 points**, a difference of **1.29 points** (95% two-sample field-bootstrap interval **0.45–2.15**). This comparison tests a spatial statistic beyond the quantities used to calibrate the original simulation. The observing team and model authors used overlapping field data. [Collins et al.](https://doi.org/10.1029/2021JE007149) · [Versioned model archive](https://doi.org/10.5281/zenodo.6566345)

| Minimum crater diameter | Measured fields | Simulated fields | Measured reduction (points) | Matched model (points) | Difference, 95% interval (points) |
|---|---:|---:|---:|---:|---:|
| 1 m | 259 | 1,600 | 3.00 | 1.70 | 1.29 [0.45, 2.15] |
| 1.5 m | 164 | 859 | 2.97 | 1.52 | 1.45 [0.37, 2.51] |
| 2 m | 99 | 478 | 2.29 | 1.06 | 1.23 [−0.11, 2.53] |
| 3 m | 37 | 209 | 1.67 | 1.31 | 0.35 [−1.86, 2.73] |

The archive supplies 115,748 model crater centers and radii. Converting radius to diameter and applying the observational selection leaves 92,474 simulated craters in the 1,600 eligible model fields for the 1 m cut. Each observed field has equal weight. For each field we average the ten nearest simulated fields in log crater count, log effective diameter (cube root of summed cubed retained diameters), and log RMS distance from the centroid. Each covariate is scaled by its simulation interquartile range. The sample medians of these quantities agree closely (22 vs. 22 craters, 8.22 vs. 8.17 m, 62.58 vs. 60.83 m). Some outlying fields have poor matches. Restricting to the 211 fields for which all ten neighbors agree within a factor of 1.5 in every covariate gives a 1.20-point difference. Matching only on crater count gives a smaller discrepancy. The JSON also reports two-covariate and unmatched results.

The two-sample percentile intervals resample whole sites independently in both libraries, repeat matching with replacement at fixed covariate scales, and use 1,999 bootstrap replicates with a recorded seed. The code also reports a fixed-simulation-library bootstrap. The alternative implementation reproduces all 559 eligible site/threshold measurements in the observational script to within 3×10⁻¹⁶ on the unit-rank scale.

The effect at the two larger diameter cuts is compatible with the matched simulation once uncertainty in the finite simulation library is included. The archived run has a 0.5 m minimum crater radius. Synthetic craters have no image-selection or diameter-measurement errors. Differences at the smallest scales could arise from image completeness, geologic target properties, sampling, or fragmentation physics. A parameter survey or image-level comparison would be necessary to distinguish these explanations.

## Reproduce

Use the observational archive and dependencies in the [measurement report](README.md). Run from this directory.

```sh
curl -fL 'https://zenodo.org/api/records/6566345/files/gsc10/Mars-clusters-v1.0.zip/content' -o model.zip
python compare_model.py cluster_data.zip model.zip rerun
```

The script verifies both archive hashes. It reads crater-level coordinates from `crater-data/craters-198.csv` and writes the eligible site measurements and comparison summaries. The DOI archive lists the model settings in `model-data/record.csv`. Other archived runs supply summary tables. Run 198 supplies the full member positions used here.
