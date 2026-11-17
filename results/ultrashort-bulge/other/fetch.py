"""Freeze public MOA/PRIME BLG alert indexes and download per-alert difference fluxes.

Run from repository root. Requires numpy and lxml. Cache is outside the results tree.
Repeated runs reuse immutable successful cache entries; choose a new cache for a new freeze.
"""
import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urljoin
from urllib.request import Request, urlopen

import numpy as np
from lxml import html

BASE = 'https://moaprime.massey.ac.nz/'
OUT = Path(__file__).resolve().parent
INDEX = {f'{s}{y}': f'alerts/index/{s}/{y}' for s in ['moa', 'prime'] for y in [2025, 2026]}
INDEX.update(moa2024='moaarchive/list/2024', prime2024idx='prime2024/idxmap')
SHA = lambda x: hashlib.sha256(x).hexdigest()


def get(url):
    for attempt in range(2):
        try:
            request = Request(url, headers={'User-Agent': 'Open-Space-research/1.0 (public alert photometry)'})
            with urlopen(request, timeout=45) as response:
                return response.read()
        except OSError:
            if attempt == 1:
                raise
            time.sleep(1)


def freeze(cache):
    rows, receipts = [], []
    for key, endpoint in INDEX.items():
        url = urljoin(BASE, endpoint)
        dest = cache / key
        if not dest.exists():
            dest.write_bytes(get(url))
        data = dest.read_bytes()
        receipts.append({'name': key, 'url': url, 'input_bytes': len(data), 'input_sha256': SHA(data),
                         'frozen_utc': datetime.fromtimestamp(dest.stat().st_mtime, timezone.utc).isoformat()})
        if key == 'moa2024':
            for event, phase, field, chip, ra, dec, note in json.loads(data)['entries']:
                if '-BLG-' in event:
                    rows.append(dict(event=event, field=field, ra=ra, dec=dec,
                                     page=urljoin(BASE, 'moaarchive/event/'+event[4:]), index=key))
        elif key == 'prime2024idx':
            for event, v in json.loads(data).items():
                if '-BLG-' in event:
                    rows.append(dict(event=event, field=v['intid'].split('-')[0], ra=v['ra'],dec=v['dec'],
                                     index_te_d=v['tE'], index_amax=v['amax'], index_t0=v['t0cal'],
                                     page=urljoin(BASE, 'prime2024/display/'+event), index=key))
        else:
            for row in html.fromstring(data).xpath('//tbody/tr'):
                event, date, ra, dec, t0, te, amax, mag = [c.text_content().strip() for c in row]
                if '-BLG-' in event:
                    rows.append(dict(event=event, alert_date=date, ra=ra, dec=dec, index_t0=t0,
                                     index_te_d=te, index_amax=amax, index_baseline_mag=mag,
                                     page=urljoin(BASE, row.xpath('.//a/@href')[0]), index=key))
    write_csv(OUT/'index_receipts.csv', receipts)
    write_csv(OUT/'index.csv', rows)
    return rows


def write_csv(path, rows):
    if not rows:
        return
    keys = list(dict.fromkeys(k for r in rows for k in r))
    tmp = path.with_suffix('.tmp')
    with tmp.open('w', newline='') as f:
        w = csv.DictWriter(f, keys, lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    tmp.replace(path)


def num(s):
    try:
        return float(s)
    except (TypeError,ValueError):
        return float('inf')


def process(row, cache):
    event = row['event']
    rec = cache/'records'/f'{event}.json'
    array = cache/'curves'/f'{event}.npz'
    if rec.exists() and array.exists():
        return json.loads(rec.read_text())
    r = {'event': event, 'event_url': row['page']}
    try:
        page = cache/'pages'/f'{event}.html'
        if not page.exists():
            page.write_bytes(get(row['page']))
        body = page.read_bytes()
        parsed = html.fromstring(body)
        r.update(page_sha256=SHA(body), page_bytes=len(body))
        meta = {c[0].text_content().strip(': '): c[1].text_content().strip() for c in parsed.xpath('//*[@id="metadata"]//tr') if len(c)>=2}
        for name, key in [('Internal ID','internal_id'),('Assessment','assessment')]:
            r[key] = meta.get(name, '')
        r['field'] = r['internal_id'].split('-')[0]
        for tr in parsed.xpath('//*[@id="micro"]//tr'):
            c = [t.text_content().strip() for t in tr]
            if len(c) >= 3 and c[0] in ['tE','tmax','u0']:
                r[{'tE':'page_te_d','tmax':'page_t0_jd','u0':'page_u0'}[c[0]]] = c[2]
        href = parsed.xpath('//a[contains(translate(string(.),"P","p"),"photometry")]/@href')
        if not href:
            raise ValueError('page missing photometry link')
        r['phot_url'] = urljoin(row['page'],href[0])
        data = get(r['phot_url'])
        r.update(input_sha256=SHA(data), input_bytes=len(data))
        raw = gzip.decompress(data) if data[:2] == b'\x1f\x8b' else data
        if raw.lstrip().startswith(b'<'):
            raise ValueError('HTML returned for photometry')
        arr = np.loadtxt(io.BytesIO(raw), usecols=(0,1,2), ndmin=2)
        mask = np.isfinite(arr).all(axis=1) & (arr[:,2] > 0)
        mask &= (arr[:,0] > 2450000) & (arr[:,0] < 2500000)
        a = arr[mask]
        a = a[np.argsort(a[:,0],kind='stable')]
        if len(a) < 2 or not 2400000 < np.median(a[:,0]) < 2500000:
            raise ValueError('unexpected time scale/empty file')
        t = a[:,0]
        t0 = num(r.get('page_t0_jd'))
        gap = np.diff(t)
        near = np.abs(t-t0) < 0.5
        r.update(raw_rows=len(arr), usable=len(a), first_jd=t[0], last_jd=t[-1],
                 points_within_half_day=int(near.sum()),
                 positive_gap_lt_halfday_median_min=round(float(np.median(gap[(gap>0)&(gap<0.5)]))*1440,4) if ((gap>0)&(gap<0.5)).any() else '',
                 retrieved_utc=datetime.now(timezone.utc).isoformat(), status='ok')
        tmp = array.with_suffix('.part')
        with tmp.open('wb') as f: np.savez_compressed(f, data=a)
        os.replace(tmp,array)
    except Exception as e:
        r.update(status='error',error=f'{type(e).__name__}: {e}')
    rec.write_text(json.dumps(r,sort_keys=True)+'\n')
    return r


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cache',type=Path,default=Path('.survey-cache'))
    p.add_argument('--workers',type=int,default=8)
    p.add_argument('--inventory-only',action='store_true')
    args=p.parse_args()
    for d in ['curves','records','pages']:
        (args.cache/d).mkdir(parents=True,exist_ok=True)
    rows=freeze(args.cache)
    if not args.inventory_only:
        rows.sort(key=lambda r:(num(r.get('index_te_d'))>=0.65, num(r.get('index_te_d'))))
        with ThreadPoolExecutor(args.workers) as pool:
            futures=[pool.submit(process,r,args.cache) for r in rows]
            out=[]
            for f in as_completed(futures):
                out.append(f.result())
                if len(out) % 25 == 0:
                    write_csv(OUT/'photometry.csv', sorted(out,key=lambda r:r['event']))
                    print(len(out),sum(r['status']=='ok' for r in out),flush=True)
        write_csv(OUT/'photometry.csv', sorted(out,key=lambda r:r['event']))
