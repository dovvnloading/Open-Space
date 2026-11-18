# A fitted magnetic-field delay fails the SOLAR-1 and Wind transfer test

**Null result.** On 7 November 2026, a fixed phase-front direction fitted to earlier NOAA SOLAR-1 and NASA Wind data raised the five-minute magnetic-change correlation to **0.06335**. Simultaneous and radial comparisons were weakly negative. A circular time-shift control scored **0.10442**, so the fitted delay failed the frozen transfer gate. The fit used **134** eligible minute pairs from 6 November. The common delay mask left zero eligible pairs on 5 November. [Frozen protocol](PROTOCOL.md) · [Scores and coverage](summary.json)

## The two spacecraft sampled widely separated paths

The median transverse separation was **894,433 km** across eligible held-out minutes. The predicted radial separation was about 1,396 km. Using the same retrospective Wind plasma measurements for both models gives a median radial offset of **2.84 s** and a fitted offset of **−1817.47 s**, measured from the SOLAR-1 minute to the predicted Wind minute. Negative offsets place Wind earlier. The small radial offset is below the one-minute field resolution. Both spacecraft locations come from predicted ephemerides. Independently transforming the paired Wind inertial coordinates recreates its published GSE coordinates within 1.074 km per component, which checks the axis convention. It gives no estimate of orbit prediction error. [Input byte receipts](sources.json) · [NOAA SOLAR-1 archive](https://www.ncei.noaa.gov/products/space-weather/swfo)

## Held-out comparison

The score is the centered, normalized dot product of five-minute changes of the three GSE components. Each row compares the models on identical eligible minutes. The frozen slope is 1.1 in the GSE XY plane. Before fitting, the design reserved 7 November for magnetic scoring and searched a fixed range of effective phase-front directions on the earlier dates.

| UTC date | Eligible one-minute pairs | Simultaneous | Radial | Fitted |
| --- | ---: | ---: | ---: | ---: |
| 5 November | 0 | — | — | — |
| 6 November | 134 | 0.02836 | 0.01555 | 0.49872 |
| 7 November | **870** | −0.08084 | −0.08032 | **0.06335** |

The fitted slope beat both baselines in **4 of 6** qualifying two-hour blocks. The largest correlation across the **41** circular-shift controls and their three model scores was **0.10442**. The fitted score with a fixed radial flow of −450 km/s fell to −0.04554 on the same 870 minutes. The circular shifts are descriptive because these records are serially correlated and the controls reuse the same observations.

SOLAR-1 alone selected five large local rotations on 7 November. Its largest rotation was **130.16°** near 06:46:30 UTC. At the fitted prediction of 06:13:30 UTC, Wind rotated **6.23°**. The nearby time search can produce closely aligned change directions even for much smaller rotations. These diagnostic matches carry no confirmatory weight. [Rotation and shift-control records](summary.json)

## Coverage loss

The primary design required every slope in the training grid to be valid on the same minutes. On 5 November, every finite, flow-qualified minute failed the 60-minute cap for at least one trial slope. On 6 November, 1,141 such minutes failed the cap. Requiring geometry at both endpoints on that day left 282 changes. Missing SOLAR-1 changes removed four and unavailable Wind samples or day boundaries removed 144 more. This is a limitation of the model comparison design. The native magnetic products separately supplied over 1,400 quality-qualified minutes per spacecraft each day.

A post-hoc diagnostic held the fitted slope fixed and imposed shared coverage across only the three reported models. Coverage increased to 1,361, 1,299 and 1,319 pairs across the three days. The fitted scores were −0.03142, 0.12466 and 0.04269. This diagnostic changes the evaluation subset and uses the already fitted slope. [Coverage diagnostic](coverage.json)

SOLAR-1 MAG and Wind MFI measure the magnetic field independently. The Wind SWE flow includes predominantly limited-tracking retrievals under the mission quality guide. Magnetic calibration, time alignment, the predicted positions and true variations between the sampled solar-wind paths can all affect the comparison. The test covers these particular days and the fixed model family. Prior phase-front work includes a published SOLAR-1/Wind comparison during the 19 January 2026 shock and a separate ACE/DSCOVR statistical study. A targeted search found no study of this November pairwise timing window. [NASA SWE quality guide](https://cdaweb.gsfc.nasa.gov/wind_swe_quality.html) · [Weimer et al. on phase-front timing](https://doi.org/10.1029/2001JA009102) · [Raptis et al. on the January shock](https://arxiv.org/abs/2611.18314) · [ACE/DSCOVR timing study](https://doi.org/10.1029/2026SW005102)

## Reproduce

Run the offline analysis from this result folder with Python and NumPy. The compressed bundle stores quality-filtered one-minute magnetic vectors, flow medians, counts and predicted positions. An independent Windows replay from the native products reproduced the field vectors and coverage exactly, with sub-micrometre numerical differences in projected SOLAR-1 positions. All reported statistics agreed at their printed precision. The offline replay from the deposited minute bundle reproduced the summary and coverage files byte for byte in the same environment.

```sh
python analyze.py --output summary.json
python coverage_diagnostic.py --output coverage.json
```

For native replay, download the thirteen files in the [source receipts](sources.json) to a directory with the specified mission subfolders. Run [build_inputs.py](build_inputs.py) with that directory as the data root and a path for the minute bundle as the output. The decoder checks source hashes before reading any science arrays. Native replay uses NumPy, SciPy, Astropy, cdflib and netCDF4; the original products were staged on the D drive.
