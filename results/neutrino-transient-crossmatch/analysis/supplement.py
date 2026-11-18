"""Describe independently published IR-flare hosts with W2 maxima in 2024, as a post-hoc coverage check."""
import csv
import gzip
import json
import numpy as np

from run import ROOT, days, rows, prepared_region, CONTROL_POSITIONS, xyz


def main():
    cat = rows('gammair/tdecat/observed_maxima.csv')
    gcn = rows('ice/gcn_alerts.csv')
    native = [json.loads(q) for q in (ROOT/'ice/gcn_hpd90.jsonl').read_text().splitlines()]
    with gzip.open(ROOT/'ice/archival_maps.jsonl.gz', 'rt') as f:
        native.extend(json.loads(q) for q in f)
    output = []
    for s in cat:
        t = days(s['peak_utc'])
        if t < days('2024-01-01T00:00:00Z'):
            continue
        ra, dec = float(s['ra_deg']), float(s['dec_deg'])
        if np.any(np.linalg.norm(xyz(*CONTROL_POSITIONS.T)-xyz(ra, dec), axis=1) <= 2*np.sin(np.deg2rad(5/3600)/2)):
            continue
        nearby = [m for m in native if abs(t-days(m['utc'])) <= 180]
        hits = [m['alert_id'] for m in nearby if prepared_region(m)(np.array([ra]), np.array([dec]))[0]]
        unresolved = []
        for a in gcn:
            if a['status'] != 'active' or a['region'] != 'circular_bounds_only' or abs(t-days(a['utc'])) > 180:
                continue
            dx = (ra-float(a['ra_deg'])+180) % 360-180
            dy = dec-float(a['dec_deg'])
            if -float(a['ra_minus90']) <= dx <= float(a['ra_plus90']) and -float(a['dec_minus90']) <= dy <= float(a['dec_plus90']):
                unresolved.append(a['alert_id'])
        output.append(dict(source=s['source_id'], peak_utc=s['peak_utc'], right_censored=s['right_censored'], native_time_eligible=len(
            nearby), native_pairs='|'.join(hits), legacy_envelope_pairs='|'.join(unresolved)))
    with (ROOT/'analysis/modern/supplement_tdecat.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(output[0]), lineterminator='\n')
        w.writeheader()
        w.writerows(output)
    print('Supplementary W2-host checks:', len(output))


if __name__ == '__main__':
    main()
