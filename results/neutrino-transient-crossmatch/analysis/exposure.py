"""Descriptive RA-exposure and season diagnostics separate from the registered tests."""
import json
import gzip
from pathlib import Path
from datetime import datetime
from collections import Counter

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def moment_score(ra):
    # First three Fourier harmonics on the circle, maximum across harmonics.
    # Scaling is common to observed and simulated values for each cohort.
    r = np.radians(ra)
    return np.maximum.reduce([np.abs(np.exp(1j * k * r).mean(axis=-1)) for k in (1, 2, 3)])


def main():
    with gzip.open(ROOT / 'ice/archival_maps.jsonl.gz', 'rt') as f:
        historic = [json.loads(row) for row in f]
    with (ROOT / 'ice/gcn_hpd90.jsonl').open() as f:
        modern = [json.loads(row) for row in f]
    rng = np.random.default_rng(5436901)
    checks = []
    for label, entries in [('archival', historic), ('GCN_native', modern)]:
        ra = np.array([r['ra_deg'] for r in entries])
        d = [datetime.fromisoformat(r['utc'].replace('Z', '+00:00')) for r in entries]
        observed = float(moment_score(ra))
        exceed = 0
        for size in [500] * 19 + [499]:
            exceed += int(np.sum(moment_score(rng.uniform(0, 360, (size, len(ra)))) >= observed))
        checks.append(dict(sample=label, alerts=len(ra), circle_max_three_harmonics=observed, uniform_ra_diagnostic_p=(exceed+1)/10000,
                      ra_quadrants=np.histogram(ra, bins=[0, 90, 180, 270, 360])[0].tolist(), calendar_months=dict(sorted(Counter(t.month for t in d).items()))))
    output = Path(__file__).with_name('exposure.json')
    output.write_text(json.dumps(checks, indent=2)+'\n', newline='\n')


if __name__ == '__main__':
    main()
