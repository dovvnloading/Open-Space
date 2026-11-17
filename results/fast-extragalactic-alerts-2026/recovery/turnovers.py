"""Measure how often frozen alert epochs can resolve a rise and fade within a day."""
import csv
import hashlib
import itertools
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

from simulate import AB_ZERO_NJY, MAG_FACTOR, ROOT, load

HERE = Path(__file__).resolve().parent
PHASES = 2048
SEEDS = {'rubin': 73117, 'ztf': 73118}


def turnover_visits(records):
    """Return each distinct visit used by a same-band triple and the triple indices."""
    by_band = defaultdict(list)
    for (band, visit), (time, error) in records.items():
        by_band[band].append((time, visit, error))
    visits = []
    triples = []
    for band in sorted(by_band):
        ordered = sorted(by_band[band])
        offset = len(visits)
        visits.extend((t, e) for t, _, e in ordered)
        for start, (t, _, _) in enumerate(ordered):
            later = []
            for index in range(start + 1, len(ordered)):
                gap = ordered[index][0] - t
                if gap > 1:
                    break
                if gap >= 0.02:
                    later.append(index)
            for middle, end in itertools.combinations(later, 2):
                if ordered[end][0] - ordered[middle][0] >= 0.02:
                    triples.append((offset + start, offset + middle, offset + end))
    if not triples:
        return None
    used, remap = np.unique(np.asarray(triples).ravel(), return_inverse=True)
    times, sigmas = np.asarray(visits)[used].T
    return times, sigmas, remap.reshape(-1, 3).T


def significant_brightening(faint, bright, faint_err, bright_err, gap):
    # Arrays carry simulated phases in rows and candidate triples in columns.
    # Clipping prevents log/division errors on fluxes that already fail S/N.
    positive = (faint >= 5 * faint_err) & (bright >= 5 * bright_err)
    faint = np.maximum(faint, faint_err)
    bright = np.maximum(bright, bright_err)
    change = 2.5 * np.log10(bright / faint)
    error = MAG_FACTOR * np.hypot(faint_err / faint, bright_err / bright)
    return positive & (change >= 0.30) & (change > gap) & (change >= 3 * error)


def recovery(survey):
    first, data = load(survey)
    models = list(itertools.product((18, 20, 22), (0, 2, 4)))
    recovered = defaultdict(int)
    bracketed = defaultdict(int)
    variance = defaultdict(float)
    n_capable = 0
    n_triples = 0
    rng = np.random.default_rng(SEEDS[survey])
    for oid in sorted(first):
        geometry = turnover_visits(data[oid])
        if geometry is None:
            continue
        n_capable += 1
        time, sigma, (left, middle, right) = geometry
        time = time - first[oid]
        n_triples += len(left)
        peak = rng.uniform(-0.5, 0.5, (PHASES, 1))
        distance = np.abs(time - peak)
        noise = rng.normal(0, sigma, (PHASES, len(sigma)))
        contains_peak = (time[left] <= peak) & (peak <= time[right])
        for mag, slope in models:
            flux = 10 ** (0.4 * (AB_ZERO_NJY - mag - slope * distance)) + noise
            rise = significant_brightening(flux[:, left], flux[:, middle], sigma[left],
                                          sigma[middle], time[middle] - time[left])
            fade = significant_brightening(flux[:, right], flux[:, middle], sigma[right],
                                          sigma[middle], time[right] - time[middle])
            shape = rise & fade
            count = np.count_nonzero(shape.any(axis=1))
            recovered[mag, slope] += count
            if slope:
                bracketed[mag, slope] += np.count_nonzero((shape & contains_peak).any(axis=1))
            # Independent phases per object. Variance of this stratum's mean.
            p = count / PHASES
            variance[mag, slope] += p * (1 - p) / (PHASES - 1)
    rows = []
    for model in models:
        mag, slope = model
        count = recovered[model]
        denominator = PHASES * len(first)
        zero_upper = ''
        if count == 0:
            # Uniform phases within each capable object. The arithmetic/geometric
            # mean bound handles heterogeneous object-specific probabilities.
            zero_upper = 100 * n_capable / len(first) * (-math.expm1(math.log(0.05) / (PHASES * n_capable)))
        rows.append(dict(survey=survey, peak_AB=mag, slope_mag_day=slope, templates=len(first),
                         timing_capable=n_capable, distinct_triples=n_triples,
                         phases_per_template=PHASES, denominator=denominator,
                         measured_turnover=count, true_peak_bracketed=bracketed[model] if slope else '',
                         percent_all=100 * count / denominator,
                         percent_capable=100 * count / (PHASES * n_capable),
                         mc_se_percentage_points=100 * math.sqrt(variance[model]) / len(first) if count else '',
                         zero_count_upper95_percent_all=zero_upper))
    return rows


def observed(survey):
    if survey == 'rubin':
        path = ROOT / 'rubin/pairs.csv'
        band, a, b, t_a, t_b, dm = 'band', 'measurement1', 'measurement2', 't1', 't2', 'dm_mag'
    else:
        path = ROOT / 'ztf/derived/ztf_pairs.csv'
        band, a, b, t_a, t_b, dm = 'fid', 'candid1', 'candid2', 'mjd1', 'mjd2', 'dm'
    rising = defaultdict(list)
    fading = defaultdict(list)
    with path.open() as file:
        pairs = list(csv.DictReader(file))
    for row in pairs:
        if float(row[dm]) < 0:
            rising[row['oid'], row[band], row[b]].append(row)
        if float(row[dm]) > 0:
            fading[row['oid'], row[band], row[a]].append(row)
    triplets = []
    for key in sorted(rising.keys() & fading.keys()):
        for rise, fade in itertools.product(rising[key], fading[key]):
            if float(fade[t_b]) - float(rise[t_a]) <= 1:
                triplets.append(dict(oid=key[0], band=key[1], measurement_ids=[rise[a], key[2], fade[b]]))
    # Two distinct visits at the top of a light curve can also separate its limbs.
    ordered = set()
    for rise in pairs:
        for fade in pairs:
            if (rise['oid'] == fade['oid'] and rise[band] == fade[band]
                    and float(rise[dm]) < 0 < float(fade[dm])
                    and float(rise[t_b]) <= float(fade[t_a])
                    and float(fade[t_b]) - float(rise[t_a]) <= 1):
                ordered.add(rise['oid'])
    return dict(objects_with_ordered_rise_and_fade=sorted(ordered),
                passing_pairs=len(pairs), rising_pairs=sum(map(len, rising.values())),
                fading_pairs=sum(map(len, fading.values())), triplets=triplets,
                objects_with_turnover=len({t['oid'] for t in triplets}))


def main():
    rows = recovery('rubin') + recovery('ztf')
    with (HERE / 'turnover_recovery.csv').open('w') as file:
        writer = csv.DictWriter(file, fieldnames=rows[0], lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    files = list(ROOT.glob('rubin/data/day_*.csv.gz')) + list(ROOT.glob('ztf/data/ztf/*_detections.csv.gz'))
    files.extend([ROOT / 'rubin/pairs.csv', ROOT / 'ztf/derived/ztf_pairs.csv'])
    audit = dict(seeds=SEEDS, phases_per_template=PHASES,
                 observed={s: observed(s) for s in SEEDS},
                 sha256_input_files={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                                     for p in sorted(files)})
    (HERE / 'turnover_audit.json').write_text(json.dumps(audit, indent=2) + '\n')


if __name__ == '__main__':
    main()
