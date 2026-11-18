"""Measure the source coverage of the supplemental legacy-box pair, AT2022fpx."""
import csv
import hashlib
import io
import json
import math
import sys
from collections import defaultdict
from pathlib import Path
from statistics import median
from urllib.request import urlopen

from run import ROOT, days


def main():
    receipt = next(q for q in json.loads((ROOT/'gammair/tdecat/receipts.json').read_text())
                   if '/target_AT2022fpx_photometry.csv' in q['url'])
    data = Path(sys.argv[1]).read_bytes() if len(
        sys.argv) > 1 else urlopen(receipt['url'], timeout=45).read()
    assert hashlib.sha256(data).hexdigest() == receipt['sha256']
    a = next(r for r in csv.DictReader((ROOT/'ice/gcn_alerts.csv').open())
             if r['alert_id'] == 'IC240307A')
    m = days(a['utc'])+40587
    groups = defaultdict(list)
    for r in csv.DictReader(io.StringIO(data.decode()), delimiter=';'):
        try:
            t, mag, e = [float(r[k]) for k in ('MJD', 'Magnitude', 'Error')]
        except ValueError:
            continue
        if abs(t-m) <= 180 and math.isfinite(mag) and e > 0 and math.isfinite(e):
            groups[(r['Facility'], r['Filter'])].append((t, mag, e))
    out = []
    for (facility, band), data in sorted(groups.items()):
        t, mag, e = min(data, key=lambda q: abs(q[0]-m))
        out.append(dict(facility=facility, band=band, n_rows_180d=len(data), first_mjd=min(q[0] for q in data), last_mjd=max(
            q[0] for q in data), nearest_mjd=t, nearest_delta_days=round(t-m, 5), nearest_mag=mag, nearest_mag_error=e))
    with (ROOT/'analysis/modern/fpx_coverage.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=out[0], lineterminator='\n')
        w.writeheader()
        w.writerows(out)
    print('AT2022fpx light curve:', len(out), 'facility/band combinations')


if __name__ == '__main__':
    main()
