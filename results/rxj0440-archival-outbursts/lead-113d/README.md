# Audit of the approximate 113-day pulsar lead

The handoff supplied a period, without a source name, dates or observable. Two interpretations can be tested from public observations. Neither search identified a previously unreported 113-day recurrence. The archival flare in the parent directory came from a separate review of GBM pulsar histories.

## X-ray brightness recurrence

A sample of 76 unique Swift/BAT daily monitor products was assembled from the detected-source list (source types containing HMXB, Pulsar, Be/X or SGXB) and the GBM accreting-pulsar source table. Three existing repository targets were excluded: GX 1+4, 4U 1538−522 and SMC X-1. After daily quality cuts and the minimum sampling requirements, 68 light curves qualified. This is a targeted sample, not the entire X-ray sky. [Swift source list](https://swift.gsfc.nasa.gov/results/transients/BAT_detected.html) · [GBM source list](https://gammaray.msfc.nasa.gov/gbm/science/pulsars.html) · [All 76 products, URLs and SHA-256 digests](bat-inputs.csv)

A sinusoidal search scanned 181 evenly spaced frequencies spanning 108–118 days. The portion before MJD 57400 selected a period and phase; later data tested brightness at that fixed phase. Swift bins required DATA_FLAG = 0 and 0 < error < 0.003 counts cm⁻² s⁻¹. Each source required 1,200 qualifying days overall and 500 on each side of the split. Brightness was centered on the early median and clipped at the early 98th percentile of its absolute deviation. Inverse-variance weights were capped at the early 80th percentile. Complex Fourier sums were normalized by the root sum of squared weighted residuals. The held-out statistic is the projection in the direction of the early Fourier sum.

Rotating each held-out calendar year's complex contribution by a random phase preserves its internally correlated data, variability and seasonal gaps. The table reports 30,000 independent phase rotations and a Bonferroni correction across 68 targets. The null models loss of phase agreement between years. This is a retrospective chronological split with an unknown original selection history. The initial exploratory screen used full-series clipping and weight caps, giving 111.30 and 117.52 days for its first two targets; the released analysis derives these preprocessing choices entirely from the earlier data.

| Source | Early period (days) | Later projection | Year-rotation p | Across-target p | Comparison |
|---|---:|---:|---:|---:|---|
| A 0535+26 | 111.30 | 8.14 | 0.00070 | 0.0476 | Published 111.0 ± 0.4-day X-ray recurrence |
| IGR J19294+1816 | 117.58 | 6.16 | 0.00160 | 0.1088 | Published 116.2 ± 0.6-day activity |
| AX J1749.1−2639 / GRO J1750−27 | 114.93 | 3.98 | 0.0513 | 1.0000 | Published orbit of 29.80 days; near-115-day brightness evidence is weak |

The known A 0535+26 recurrence was the only case with corrected p below 0.05. Its adjusted p is near the threshold and varies with finite Monte Carlo sampling. The ordered table retains the other 65 sources. [A 0535+26 observations](https://doi.org/10.1038/303681a0) · [IGR J19294+1816 observations](https://doi.org/10.1051/0004-6361/201116729) · [GRO J1750−27 orbit](https://gammaray.msfc.nasa.gov/gbm/science/pulsars/lightcurves/groj1750.html) · [All BAT scores](bat-scores.csv)

## Spin-glitch recurrence

The Jodrell Bank snapshot lists 737 glitches in 226 objects. Fourteen pulsars with at least eight published glitch epochs contributed 301 events. The analysis used the first two thirds of each pulsar's chronological events to maximize phase concentration over 401 frequencies in the 108–118-day interval. The remaining events were projected onto the selected phase direction. The greatest mean held-out projection was 0.383 on a scale from −1 to +1. The smallest of the 14 simple uniform-phase tail estimates was 0.130 before multiplicity correction. At the period nearest 113 days, J1801−2304 had a 113.054-day early optimum and a held-out projection of 0.342 from five events (raw p = 0.147). Only the point estimates of published epochs enter this diagnostic. Irregular monitoring, epoch uncertainty and missing glitches limit sensitivity, and the catalog does not supply continuous observing windows. [Jodrell catalog](https://www.jb.man.ac.uk/pulsar/glitches/gTable.html) · [Frozen eligible epochs](glitch-epochs.csv) · [All glitch results](glitch-scores.csv)

## Replay

Run from this directory with NumPy. The command reads the frozen glitch epochs, downloads missing BAT products, verifies their hashes, and writes two result CSVs. The original Jodrell HTML SHA-256 is `e811419f2037080bd2c99c1db30555dd4ee6f462f36b91796908f0fe90aabeb4`. All numerical outputs for the dated Swift snapshot are preserved. The mission overwrites live monitor products, so reproducing a result after an update requires matching cached originals.

```sh
OPENBLAS_NUM_THREADS=1 python screen.py /path/to/bat-cache /tmp/pulsar-113d-replay
diff -u glitch-scores.csv /tmp/pulsar-113d-replay/glitch-scores.csv
diff -u bat-scores.csv /tmp/pulsar-113d-replay/bat-scores.csv
```
