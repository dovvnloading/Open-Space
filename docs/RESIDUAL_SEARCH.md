# Blind search of galaxy-model residual spectra

**Subsequent plate-wide test.** At the wavelength of the 9228.74 Å residual, all 82 bright comparison galaxies on the same visit showed a positive feature. The median was 10.58% of the continuum. [Independent controls and comparison plates](SDSS_SUPERPOSITION.md#a-different-result-for-the-earlier-922874-å-residual)

We searched the entire usable wavelength range of 24 galaxies on three SDSS plate visits after subtracting the survey's fitted spectral model. The 11 residual detections outside the bundled line tables produced **one same-sign rest-frame pair**, found on a single visit. **Zero groups cleared the triage rule** of at least three distinct galaxies, two plates and two observing nights. The search covered 85,881 possible feature centers across the 24 spectra. [All detections outside the tables, group members, masks and controls](../results/residual-search.json) · [Pinned targets and source hashes](../data/holdout-targets.csv)

| Whole-spectrum measurement | Original coadd | Coadd minus survey model |
| --- | ---: | ---: |
| Searchable feature centers | 85,881 | 85,881 |
| Detected narrow features | 238 | 28 |
| Within tolerance of a source-table wavelength | 68 | 17 |
| Within tolerance of a bundled observer-frame sky line | 0 | 0 |
| Outside both reference tables | 170 | 11 |

The residual run scans the spectral pixels afresh. Its 28 detections can differ in wavelength or sign from the 238 raw detections. Their count difference measures how the screen responds to model subtraction. Counting removed lines by identity would require a separate feature-matching analysis.

## Search rule

We used all 24 targets selected by archive metadata. The exploratory scan used no input candidate wavelengths. Before reading its residual detections we fixed the pipeline to the existing narrow-feature criterion, score threshold six, source and sky catalogs, and maximum group span of 150 km/s. For each target the code verifies the manifest, normalized pixel digest and raw FITS digest against the roster. It checks the model and flux grids pixel by pixel, subtracts `COADD.model` from the calibrated flux density, and retains the original nonzero-weight and AND-mask rules. Recurrence is computed separately for emission and absorption and separately in the rest and observer frames. The strongest representative from a sky position gets one vote. The output records every group with at least two targets and whether it meets the three-target, two-plate, two-night triage rule.

For each group we report the number of targets with at least one searchable center within ±150 km/s of its median wavelength, in the applicable frame. The count measures opportunity to detect a feature at that location. Intrinsic line strength, flux uncertainties and spectral resolution still affect sensitivity. This full-range search shares the cohort used for the frozen-wavelength prediction test. It is an exploratory use of those observations. [Selection and limitations](HOLDOUT_DATA.md) · [Recurrence method](RECURRENCE.md)

## Inspect the surviving pair

Two galaxies on plate 6509, MJD 56486, showed positive residual features at observed wavelengths 9227.54 and 9229.94 Å. They also met the rest-frame grouping rule at 8552.83 and 8556.12 Å, a span of 115.17 km/s. Ten of the 24 target grids had searchable coverage around the group's rest median of 8554.48 Å. The very same two detections form the observer-frame group at 9228.74 Å, where 23 of 24 grids had coverage. Their similar target redshifts make their alignment in both frames possible.

In both fibers the coadd's subtracted-sky estimate is elevated at the residual peak. At the center pixel versus the median of 16 continuum-sideband pixels it was 7.9 versus 1.9 in fiber 0128 and 6.1 versus 2.0 in fiber 0813, in the stored sky-spectrum units. One core pixel in fiber 0813 carried a nonzero OR mask. A stricter search excluding all OR flags lost searchable coverage in the neighborhood of each member. That test is inconclusive. These observations support sky subtraction as an alternative explanation. A specific sky transition remains unassigned. [SDSS sky subtraction and calibration](https://www.sdss4.org/dr17/spectro/pipeline/) · [Known spectral caveats](https://www.sdss4.org/dr17/spectro/caveats/)

| Observer-frame group (Å) | Sign | Targets with detections / targets with coverage | Observed centroid span (Å) | Plate visits with detections |
| ---: | --- | ---: | ---: | ---: |
| 8999.07 | Emission | 2 / 2 | 0.061 | 1 |
| 9228.74 | Emission | 2 / 23 | 2.400 | 1 |
| 9028.07 | Absorption | 2 / 11 | 0.066 | 1 |

For all 11 unassigned detections, applying the OR mask removed sufficient neighboring pixels to leave zero searchable centers in each original ±150 km/s window. Only three observer-frame groups and one rest-frame group formed under the default mask. An independent set of exposures and sky fibers at the same observed wavelengths would better distinguish sky residuals from persistent source features.

The spectral model was fitted to these same survey data. It can absorb a genuine feature, and imperfect modeling can create a residual. The scan preserves the input inverse variances, with no propagation of template-fit uncertainty or spectral covariance. Thus the zero qualifying result describes this detector, sample and coverage. The experiment supplies no abundance limit, per-feature false-alarm probability or evidence for a new atomic species. [SDSS spectral fitting](https://www.sdss4.org/dr17/algorithms/redshifts/)

## Reproduce

From the repository root, use a Python environment with the archive extra. The fetch step verifies or downloads the 24 pinned FITS files, manifests and normalized spectra, and the search then runs offline.

```sh
python -m pip install -e '.[archive]'
python -m examples.fetch_holdout fetch --output runs/holdout-inputs
python -m examples.residual_search runs/holdout-inputs > runs/residual-search.json
```

Feature IDs in this result refer to the residual scan's wavelength-ordered feature list. The result stores the roster and source hashes, individual target counts, all eleven out-of-table detections with quality checks and the full membership of each group.
