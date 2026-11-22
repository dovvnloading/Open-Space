# A sustained gamma-ray decline in 3C 454.3 and an incomplete periodicity search

**3C 454.3 (4FGL J2253.9+1609) sustained a decline of at least 17.99-fold across the monthly and weekly Fermi-LAT measurements used here.** It was the only source to clear the tenfold state-change threshold among 1,459 with adequate two-cadence coverage in the completed batches. The catalog-wide periodicity calculation did not finish. No source has a reportable global periodicity significance or a discovery designation from this investigation. [Measured bounds](tables/state-transitions-emitted.csv) · [Coverage](tables/counts.json)

## What Fermi measured

Two adjacent trial dates passed for 3C 454.3, from MJD 55686.44 through 55716.88 in 2011. At the stronger split, MJD 55716.875, the monthly contrast bound was 30.49 and the weekly bound was 17.99. Each of the two half-years before this split had to exceed both half-years afterward by a factor of ten. These windows establish persistence for a year on each side. They do not determine the exact onset date or show that the source ceased emitting. The cadence products share photons.

For each detected bin the census used the flux plus or minus 1.96 times its reported error, with the lower endpoint bounded by zero. An upper-limit bin ranged from zero to its catalog limit. Missing bins received open-ended bounds. A half-year needed usable measurements in at least 5 of 6 monthly or 21 of 26 weekly bins. The contrast divides the smaller bright-side lower median by the larger faint-side upper median. This is a descriptive measurement bound conditional on individual bin intervals, not a simultaneous confidence interval or a red-noise false-alarm probability. The FSSC roster identifies the counterpart as 3C 454.3. The planned ADS, ATel, Peñil et al., DR4-note and independent optical/radio checks did not run, so this entry makes no novelty claim. [All source dispositions](tables/states-emitted.csv) · [Frozen source roster](partial/lcr-index.zip)

## Timing coverage and limits

The Fermi archive's variability threshold (DR2 Variability_Index > 21.67) selected 1,525 of 5,788 sources. Of the 1,462 sources with deposited processing records, 755 had a usable timing series. Another 704 failed the detection or time-coverage requirement and three lacked usable primary data. Across the three cadences the completed batches validated 4,377 of 4,386 requested products. The 63 sources in batch 19 have no deposited timing or state results.

The blind period interval ran from six monthly bins to one-quarter of each source's baseline. Native weekly data could substitute only for an unavailable monthly file, with a 26-week minimum. The completed workers fit each eligible source's red-noise continuum, simulated 1,023 series with its measured flux distribution and observing window, searched all specified periods, and measured stability across disjoint halves. The files retain the trial maxima, independent-half checks and injection scores at 240, 540 and 1,080 days and fractional semi-amplitudes of 25%, 50% and 90%. The planned decision used each simulated catalog's maximum over all eligible sources and periods. With the final batch missing, the workflow did not calculate that threshold, a catalog-wide p-value or a recovery fraction at the threshold. The partial record supports neither a full-survey bounded null nor a binary candidate. [Emitted timing coverage](tables/period-status-emitted.csv) · [Timing model and fixed gates](periods.py)

The first workflow run checked the NASA interface. The second froze the roster and encountered an HTTP 403 on its monthly control. The third retrieved that control successfully. In the third run, the final outstanding worker reached its 100-minute acquisition and analysis job limit. GitHub skipped the dependent combination and publication jobs. The 23 completed batches were preserved from their verified Actions artifacts. The three-run budget ended there. [Final run and job records](https://github.com/dovvnloading/Open-Space/actions/runs/45054759710)

## Input checksums and reproduction

The [per-cadence table](tables/cadence-receipts.csv) gives the SHA-256 of each available FSSC response. Each completed batch ZIP contains the full URL and receipt manifest, the normalized time, flux, error and upper-limit arrays, and the simulation outputs. The SHA-256 of the original, uncompressed DR2 source-table response is `f00ab05e89dfadf270ff8ee43b28718f232abf84db87280f15d2d4cddbce3b46`. The frozen index ZIP contains that response. The [artifact manifest](partial/transport-manifest.json) records GitHub's ZIP checksums, identifiers, raw-response artifact links and expiry. The normalized curves and completed batch outputs in this directory remain usable after the Actions artifacts expire. [Fermi-LAT Light Curve Repository](https://fermi.gsfc.nasa.gov/ssc/data/access/lat/LightCurveRepository/)

The commands below verify the deposited files and replay the batch containing 3C 454.3. Run them from the repository root with Python 3.12. Substitute another deposited batch number to reproduce its measurements.

```bash
study=results/fermi-lcr-blind-survey
(cd "$study" && sha256sum -c SHA256SUMS)
python -m pip install -r "$study/requirements.txt"
work=$(mktemp -d)
unzip -q "$study/partial/lcr-index.zip" -d "$work/index"
unzip -q "$study/partial/part-12.zip" -d "$work/reference"
mkdir "$work/data"
tar -xzf "$work/reference/snapshot/curves-12.tar.gz" -C "$work/data"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python "$study/survey.py" shard --index "$work/index/catalog.json" --data "$work/data" --out "$work/replay" --shard 12 --shards 24
python -m open_space.verify
python -m unittest discover -s tests
```
