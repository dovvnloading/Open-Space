"""Recompute the fixed six-hour blocks from NASA CDAWeb HAPI CSV and metadata."""
import csv
import io
import json
from datetime import datetime, timedelta
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np

ROOT = Path(__file__).resolve().parent
BASE = 'https://cdaweb.gsfc.nasa.gov/hapi/'
CENTERS = {'A': '2014-02-02T08:00:00+00:00', 'B': '2009-11-15T02:00:00+00:00'}


def get(endpoint, **params):
    with urlopen(BASE + endpoint + '?' + urlencode(params), timeout=120) as f:
        return f.read().decode()


def series(dataset, parameters, start, end):
    stamp = lambda t: t.strftime('%Y-%m-%dT%H:%M:%SZ')
    result = get('data', id=dataset, parameters='Time,' + parameters,
                 **{'time.min': stamp(start), 'time.max': stamp(end), 'format':'csv'})
    rows = list(csv.reader(io.StringIO(result)))
    time = np.asarray([datetime.fromisoformat(r[0]).timestamp() for r in rows])
    values = np.array([[float(x) if x and float(x) > -1e30 else np.nan for x in r[1:]]
                      for r in rows])
    return time, values


def blocks(probe, center):
    center = datetime.fromisoformat(center)
    start, end = center - timedelta(days=7,hours=4), center + timedelta(days=7,hours=4)
    p = f'ST{probe}_L2_PLA_1DMAX_10MIN'
    tp, wind = series(p,'proton_number_density_10min,proton_temperature_10min',start,end)
    tb, field = series(f'ST{probe}_L2_MAGPLASMA_1M','BTOTAL',start,end)
    if probe == 'A':
        h = 'STA_L3_PLA_HEPLUS_F_VSW_01HR'
        data = json.loads(get('info',id=h))['parameters']
        bins = np.array(next(v for v in data if v['name']=='he_plus_counts')['bins'][0]['centers'])
        th, counts = series(h, 'he_plus_counts', start, end)
        counts = counts[:,(bins >= 1) & (bins <= 2)]
    for shift in range(-7,8):
        t = center + timedelta(days=shift)
        def inside(times):
            return (times >= t.timestamp()-10800) & (times < t.timestamp()+10800)
        n, temperature = wind[inside(tp)].T
        b = field[inside(tb),0]
        row = dict(probe=probe,days_from_event=shift,center_utc=t.isoformat(),
                   proton_n_cm3=round(np.nanmedian(n),3),proton_n_bins=int(np.isfinite(n).sum()),
                   temperature_K=round(np.nanmedian(temperature),1),
                   temperature_bins=int(np.isfinite(temperature).sum()),
                   B_nT=round(np.nanmedian(b),4),B_bins=int(np.isfinite(b).sum()),
                   He_counts='',He_hours='',He_counts_per_hour='')
        if probe == 'A':
            selected = counts[inside(th)]
            selected = selected[np.all(np.isfinite(selected),axis=1)]
            hours = len(selected)
            row.update(He_counts=int(selected.sum()),He_hours=hours,
                       He_counts_per_hour=round(selected.sum()/hours,3) if hours else '')
        yield row


if __name__ == '__main__':
    rows = [r for probe,t in CENTERS.items() for r in blocks(probe,t)]
    with (ROOT / 'control_blocks.csv').open('w',newline='') as f:
        w = csv.DictWriter(f, fieldnames=rows[0])
        w.writeheader()
        w.writerows(rows)
