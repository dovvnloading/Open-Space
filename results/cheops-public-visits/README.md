# No new transit candidate survived in 5,344 CHEOPS visits

0 new candidates passed automated checks and publication review. 239 additional signals require evidence or review. The inventory contained 5,618 public visit identifiers. The workers returned 5,618 visit records and searched 5,344 visits covering 60,299 distinct visit/source pairs.
This is a null for validated candidates in the searched light curves. The injection table measures recovery at the photometric stage. Open signals and incomplete checks prevent a limit on the number of planets or eclipsing systems.

## Events

[The event table](outputs/events.tsv) gives the star, visit, BJD TDB epoch, fractional-depth equivalent in ppm, duration, reduction checks, image result, TESS check and disposition for every distinct signal. Events from multiple apertures at the same source and epoch are grouped.

## Coverage and exclusions

| Extraction | Distinct visit/source pairs |
| --- | ---: |
| field, field aperture | 56,223 |
| target, DRP | 4,076 |

The searched light curves have the following distinct source identifiers, 10,730 field, 658 target. Target aliases can refer to the same star across visits.

| Field extraction outcome | Visits |
| --- | ---: |
| catalogue apertures extracted within stated spatial and sampling limits | 5,561 |
| repair incomplete | 34 |
| missing calibrated subarray, catalogue or DRP reference | 23 |

Across visit catalogues, near target 6,948, outside 123,724,473, too few 118,891, extracted 66,845, unresolved core 152,114.

| Product excluded from search | Count |
| --- | ---: |
| field, insufficient_sampling | 10,622 |
| target, insufficient_sampling | 1,519 |

The first and last accepted product spans run from BJD TDB 2458901.773906 to 2461359.793711. These endpoints include gaps between and within visits.

0 inventory visits have no worker record. 274 recorded visits have no searched light curve. 5,344 visits have partial product, field or ephemeris coverage. DACE listing completeness was True. ESA listing completeness was False.

| Search limitation | Count |
| --- | ---: |
| Unresolved catalog query or known ephemeris (source-visits) | 72,440 |
| No searched PIPE light curve | 5,618 |
| No searched DRP light curve | 1,542 |
| Field photometry unavailable or incomplete | 57 |

| Signal disposition | Count |
| --- | ---: |
| catalog ephemeris overlap | 95 |
| photometric rejection | 8,030 |
| unresolved | 239 |

Unresolved signals can have multiple evidence limitations.

| Unresolved signal check | Count |
| --- | ---: |
| Comparison between reductions: single reduction family | 239 |
| Image localization: insufficient_roll_matched_pairs | 90 |
| Image localization: unstable_pixel_deficit | 70 |
| Image localization: insufficient_reference_frames | 66 |
| Image localization: consistent | 13 |

## Injection recovery

Each selected visit/source product received trapezoidal dips before detrending, at epochs drawn uniformly across its quality-filtered visit span, including gaps and known-event masks. Two trials were drawn for each depth and duration in each product. Every trial counts in the denominator. Recovery requires a consistent epoch and duration and all automatic photometric cuts. The table reports recovered/tried and the observed recovery fraction.

Each source contributes one product per visit and reduction family, preferring DEFAULT, then OPTIMAL, then filename order. DRP and PIPE can cover the same photons. Curves from distinct field stars in the same visit can share noise. Injections into the extracted field-star curves measure detectability after aperture extraction. A field star needed at least 48 valid cadences and a measurable core in the sampled frames. Field apertures used a radius of 10 pixels and a local sky annulus of 14 to 18 pixels. The smear and contamination regressors come from the visit’s DRP target telemetry. The recovery fractions measure the automated photometric stage. Catalog completeness, image localization and publication review are separate filters.

| Source and extraction | Depth (ppm) | 1 h | 3 h | 6 h |
| --- | ---: | ---: | ---: | ---: |
| field, field aperture | 500 | 0/112,446 (0.0%) | 0/112,446 (0.0%) | 1/112,446 (0.0%) |
| field, field aperture | 1,500 | 0/112,446 (0.0%) | 1/112,446 (0.0%) | 6/112,446 (0.0%) |
| field, field aperture | 5,000 | 0/112,446 (0.0%) | 13/112,446 (0.0%) | 88/112,446 (0.1%) |
| field, field aperture | 20,000 | 273/112,446 (0.2%) | 786/112,446 (0.7%) | 986/112,446 (0.9%) |
| target, DRP | 500 | 54/8,152 (0.7%) | 126/8,152 (1.5%) | 137/8,152 (1.7%) |
| target, DRP | 1,500 | 217/8,152 (2.7%) | 231/8,152 (2.8%) | 217/8,152 (2.7%) |
| target, DRP | 5,000 | 453/8,152 (5.6%) | 284/8,152 (3.5%) | 278/8,152 (3.4%) |
| target, DRP | 20,000 | 632/8,152 (7.8%) | 304/8,152 (3.7%) | 307/8,152 (3.8%) |

## Search and reproduction

The visit search masked published transit and eclipse windows, then fit a time trend, four Fourier harmonics of spacecraft roll, background, smear and contamination with iterative outlier rejection. A sliding box compared its central depth against both adjacent flanks. Durations ranged from 0.5 to 8 hours, with a nominal threshold of 8 times the larger of the white-noise floor and the empirical scatter of the trial depths at that duration. The search required at least four cadences in each interval, including both flanks.

The signal checks tested persistence after deleting the deepest cadence, baseline agreement, recurrence at the same roll angle and auxiliary changes. Independent extraction and pixel-difference centroid checks, nearby variables, Gaia neighbors, NASA planet and TOI ephemerides, available TESS data and the literature constrain identification. The recorded depth and duration describe the best box. The sampling cuts limit sensitivity at the edges of a visit and across data gaps. Field apertures require calibrated images, catalog positions, enough valid cadences and separation from the target. Unresolved blends, catalog omissions, faint uncataloged sources and signals outside the searched durations remain possible.

The frozen visit list has SHA-256 `6947c81aa1eacdc73802318251c427c49053296faf862875d240422975fbdd1e`. [Input digests](outputs/input_checksums.jsonl.gz) identify the retrieved archive responses, source masks, light curves and derived apertures. The visit records carry request URLs, product checksums, masks, diagnostics and per-event evidence. [Output checksums](outputs/SHA256SUMS) cover every deposited output. [Missing visit IDs](outputs/missing_visits.txt) and [all distinct events](outputs/events.tsv) are machine-readable.

## Reproduce the result

The [first run](https://github.com/dovvnloading/Open-Space/actions/runs/45054151730) deposited the frozen inventory. The [target survey](https://github.com/dovvnloading/Open-Space/actions/runs/45054699737) contains the original target measurements and catalog snapshot. The final field pass translated detector centroids into subarray coordinates and retained the original target measurements and injection trials. The public archive and external catalogs can change between downloads.

To regenerate the report from the deposited measurements in a clean checkout, run

```sh
study=results/cheops-public-visits
(cd "$study/outputs" && sha256sum -c SHA256SUMS)
python "$study/science_report.py"
```

To repeat the field pass against the same target measurements, use the following commands from the repository root. The base commit fixes the original visit ledgers. The final catalog seed fixes the published masks used for field photometry.

```sh
study=results/cheops-public-visits
scratch=$(mktemp -d)
git archive 6f0b77d2c85cb5e422d12e81e1c698f6f9f5ed92 "$study" | tar -x -C "$scratch"
base="$scratch/$study"
python -m pip install -r "$study/requirements.txt"
export CHEOPS_CATALOG_SEED="$PWD/$study/outputs/field-ephemerides/catalog-seed.json.gz"
python "$study/repair.py" "$base" "$scratch/patches" "$scratch/cache" 0 1
python "$study/repair.py" merge "$base" "$scratch/patches"
CHEOPS_OUTPUTS="$base/outputs" python "$study/science_report.py"
```

The workflow at commit `69fb1aa8c973b9ef1d749a2dbef5cdede3682a9d` provides the commands for the original DRP survey. The input receipts identify the actual bytes for comparison.
