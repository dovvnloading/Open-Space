# The CHEOPS search covered 4,076 visits

0 new candidates passed automated checks and publication review. 99 additional signals require evidence or review. The inventory contained 5,618 public visit identifiers. The workers returned 5,618 visit records and searched 4,076 visits covering 4,076 distinct visit/source pairs.
The available observations and unresolved checks do not establish an archive-wide astrophysical null.

## Events

[The event table](outputs/events.tsv) gives the star, visit, BJD TDB epoch, fractional-depth equivalent in ppm, duration, reduction checks, image result, TESS check and disposition for every distinct signal. Events from multiple apertures at the same source and epoch are grouped.

## Coverage and exclusions

| Extraction | Distinct visit/source pairs |
| --- | ---: |
| target, DRP | 4,076 |

The first and last accepted product spans run from BJD TDB 2458901.773906 to 2461359.793711. These endpoints include gaps between and within visits.

0 inventory visits have no worker record. 1,542 recorded visits have no searched light curve. 4,076 visits have partial product, field or ephemeris coverage. DACE listing completeness was True. ESA listing completeness was False.

| Search limitation | Count |
| --- | ---: |
| Field photometry unavailable or incomplete | 5,618 |
| No searched PIPE light curve | 5,618 |
| Unresolved catalog query or known ephemeris (source-visits) | 5,595 |
| No searched DRP light curve | 1,542 |

| Signal disposition | Count |
| --- | ---: |
| catalog ephemeris overlap | 54 |
| photometric rejection | 826 |
| unresolved | 99 |

Unresolved signals can have multiple evidence limitations.

| Unresolved signal check | Count |
| --- | ---: |
| Comparison between reductions: single reduction family | 99 |
| Image localization: insufficient_roll_matched_pairs | 66 |
| Image localization: insufficient_reference_frames | 33 |

## Injection recovery

Each selected visit/source product received trapezoidal dips before detrending, at epochs drawn uniformly across its quality-filtered visit span, including gaps and known-event masks. Two trials were drawn for each depth and duration in each product. Every trial counts in the denominator. Recovery requires a consistent epoch and duration and all automatic photometric cuts. The table reports recovered/tried and the observed recovery fraction.

Each source contributes one product per visit and reduction family, preferring DEFAULT, then OPTIMAL, then filename order. DRP and PIPE can cover the same photons. Curves from distinct field stars in the same visit can share noise. The recovery fractions measure the automated photometric stage. Catalog completeness, image localization and publication review are separate filters.

| Source and extraction | Depth (ppm) | 1 h | 3 h | 6 h |
| --- | ---: | ---: | ---: | ---: |
| target, DRP | 500 | 54/8,152 (0.7%) | 126/8,152 (1.5%) | 137/8,152 (1.7%) |
| target, DRP | 1,500 | 217/8,152 (2.7%) | 231/8,152 (2.8%) | 217/8,152 (2.7%) |
| target, DRP | 5,000 | 453/8,152 (5.6%) | 284/8,152 (3.5%) | 278/8,152 (3.4%) |
| target, DRP | 20,000 | 632/8,152 (7.8%) | 304/8,152 (3.7%) | 307/8,152 (3.8%) |

## Search and reproduction

The visit search masked published transit and eclipse windows, then fit a time trend, four Fourier harmonics of spacecraft roll, background, smear and contamination with iterative outlier rejection. A sliding box compared its central depth against both adjacent flanks. Durations ranged from 0.5 to 8 hours, with a nominal threshold of 8 times the larger of the white-noise floor and the empirical scatter of the trial depths at that duration. The search required at least four cadences in each interval, including both flanks.

The signal checks tested persistence after deleting the deepest cadence, baseline agreement, recurrence at the same roll angle and auxiliary changes. Independent extraction and pixel-difference centroid checks, nearby variables, Gaia neighbors, NASA planet and TOI ephemerides, available TESS data and the literature constrain identification. The recorded depth and duration describe the best box. The sampling cuts limit sensitivity at the edges of a visit and across data gaps. Field apertures require calibrated images, catalog positions, enough valid cadences and separation from the target. Unresolved blends, catalog omissions, faint uncataloged sources and signals outside the searched durations remain possible.

The frozen visit list has SHA-256 `6947c81aa1eacdc73802318251c427c49053296faf862875d240422975fbdd1e`. [Input digests](outputs/input_checksums.jsonl.gz) identify the retrieved archive responses, source masks, light curves and derived apertures. The visit records carry request URLs, product checksums, masks, diagnostics and per-event evidence. [Output checksums](outputs/SHA256SUMS) cover every deposited output. [Missing visit IDs](outputs/missing_visits.txt) and [all distinct events](outputs/events.tsv) are machine-readable.

From the repository root, a serial reproduction uses

```sh
study=results/cheops-public-visits
python -m pip install -r "$study/requirements.txt"
mkdir "$study/outputs-reproduction"
cp "$study/outputs/"{inventory.jsonl,archive_status.json,archive_probes.jsonl} "$study/outputs-reproduction/"
seed=$(find "$study/outputs" -name catalog-seed.json.gz -print -quit)
if [ -z "$seed" ]; then
  python "$study/vetting.py" seed "$study/outputs/inventory.jsonl" /tmp/cheops-catalog
  seed=/tmp/cheops-catalog/catalog-seed.json.gz
fi
export CHEOPS_CATALOG_SEED="$seed"
python "$study/worker.py" "$study/outputs-reproduction/inventory.jsonl" "$study/outputs-reproduction" /tmp/cheops-cache 0 1
CHEOPS_OUTPUTS="$study/outputs-reproduction" python "$study/science_report.py"
(cd "$study/outputs" && sha256sum -c SHA256SUMS)
```

Use a clean checkout. The serial reproduction reads the deposited inventory. The archive, catalog ephemerides and publication indexes can change. The original per-visit input digests allow comparison with a later run.
