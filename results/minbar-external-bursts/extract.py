"""Regenerate the compact measurements from the public archive products."""
import argparse
import csv
import hashlib
import io
import json
import re
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.time import Time

HERE = Path(__file__).resolve().parent
CATALOG = ('https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?'
           '-source=J%2FApJS%2F249%2F32%2Ftable11&-out.max=unlimited')
ROOT = 'https://heasarc.gsfc.nasa.gov/FTP/xte/data/archive'
GCN = 'https://gcn.gsfc.nasa.gov/notices_k'
EVENTS = [
    (2474, 4, '40033-06-02-01', 'SE7_985a0a0-985ae6c.evt.gz', '19990123.85718'),
    (2573, 4, '40048-01-09-01', 'GX_abf151e-abf1de0.evt.gz', '19990918.66011'),
    (2958, 6, '60036-01-05-00', 'GX_f2ecd60-f2ef5d0.evt.gz', '20020127.25633'),
    (3165, 8, '80091-01-08-00', 'GX_131b1c20-131b29ac.evt.gz', '20040228.00524'),
    (3234, 8, '80033-01-21-02', 'SE7_14638a4f-1463977c.evt.gz', '20041103.12510'),
    (3238, 9, '90039-01-03-05', 'SE7_14813f78-14814e5c.evt.gz', '20041125.58027'),
    (8279, 15, '96405-01-05-01', 'SE7_21b7ed80-21b80e1e.evt.gz', '20111205.47450'),
]
BINS = np.arange(-704, 704.25, 0.25)


def write_table(path, header, rows):
    with path.open('w', newline='') as f:
        writer = csv.writer(f, lineterminator='\n')
        writer.writerow(header)
        writer.writerows(rows)


def main(cache, output):
    cache.mkdir(parents=True, exist_ok=True)
    output.mkdir(parents=True, exist_ok=True)
    inputs = []

    def get(url):
        local = cache / hashlib.sha256(url.encode()).hexdigest()
        if not local.exists():
            with urllib.request.urlopen(url, timeout=120) as response:
                body = response.read()
            local.write_bytes(body)
        body = local.read_bytes()
        inputs.append({'url': url, 'sha256': hashlib.sha256(body).hexdigest()})
        return local

    cat = get(CATALOG).read_text()
    lines = [s for s in cat.splitlines() if s.strip() and not s.startswith('#')]
    rows = list(csv.DictReader(lines[:1] + lines[3:], delimiter='\t'))
    cols = ['Seq','Name','Inst','obsid','Time','sflag','kT','bpflux']
    write_table(output/'minbar.csv', cols, [[r[c].strip() for c in cols] for r in rows])
    triggers = []
    for year in range(1996, 2013):
        url = f'http://www.ioffe.ru/LEA/kw/triggers/{year}/index.html'
        page = get(url).read_text(errors='replace')
        for row in re.findall(r'<tr>.*?</tr>', page, re.S):
            cells = [re.sub(r'<[^>]+>', '', x).strip()
                     for x in re.findall(r'<td[^>]*>(.*?)</td>', row, re.S)]
            if cells and re.match(r'\d{4}-\d\d-\d\d ', cells[0]):
                triggers.append([cells[0][:19], cells[3], cells[1]])
    write_table(output/'konus_triggers.csv', ['utc','type','name'], sorted(triggers))

    profiles, rxte, flanks, metadata = [], [], [], []
    for seq, ao, obs, name, trigger in EVENTS:
        kw_file = 'konus_' + trigger
        level = 'gen' if seq == 3234 else 'auto'
        kw = get(f'{GCN}/{level}/{kw_file}')
        lines = kw.read_text().splitlines()
        trigger_s = float(lines[0].split()[3])
        band = list(map(float,lines[1].split()))
        k = np.loadtxt(io.StringIO('\n'.join(lines[4:])))
        profiles.extend((seq, format(x,'.6f'), int(y)) for x, y in k)
        date = datetime.strptime(trigger[:8], '%Y%m%d').replace(tzinfo=timezone.utc)
        # UTC date at midnight plus the trigger seconds. No leap occurred these days.
        t0 = Time(date, scale='utc').tt.mjd + trigger_s / 86400
        raw = get(f'{ROOT}/AO{ao}/P{obs[:5]}/{obs}/pca/{name}')
        with fits.open(raw) as hdus:
            h = hdus[1].header
            d = hdus[1].data
            tr = (t0-h['MJDREFI']-h['MJDREFF'])*86400 - h['TIMEZERO']
            t = d['TIME'] - tr
            p = d['PHA']
            pc = d['PCUID']
            # Different mode maps. Hard/soft channels are distinct detector tests.
            se = h['DATAMODE'].startswith('E_125us_64M')
            good = (p < (64 if se else 255)) & (pc < 5)
            hard = p >= (35 if se else 80)
            soft = p < (12 if se else 17)
            hist = [np.histogram(t[good & mask], BINS)[0] for mask in
                    [np.ones_like(good), hard, soft]]
            # Keep only bins fully within an uninterrupted exposure, avoiding
            # telemetry gaps. FITS GTI 2 describes data availability in this mode.
            valid = np.zeros(len(BINS)-1, dtype=bool)
            for start, stop in hdus[2].data:
                valid |= (BINS[:-1]+tr >= start) & (BINS[1:]+tr <= stop)
            rxte.extend([seq, f'{left:.2f}', int(a), int(b), int(c)]
                        for left,a,b,c,ok in zip(BINS[:-1], *hist,valid) if ok)
            meta = [seq,kw_file,trigger_s,*band,h['RA_PNT'],h['DEC_PNT'],h['DATAMODE']]
            metadata.append(meta)
            if seq == 8279:
                for r in rows:
                    num = int(r['Seq'])
                    if num not in (8278,8279,8280):
                        continue
                    dt = t - ((Time(float(r['Time']),format='mjd',scale='utc').tt.mjd-t0)*86400)
                    for detector in (2,4):
                        for label,bandpass in [('soft',soft),('hard',hard)]:
                            sel = good & (pc == detector) & bandpass
                            on = int(np.sum(sel & (dt>=0) & (dt<20)))
                            off = int(np.sum(sel & (dt>=-120) & (dt< -20)))
                            flanks.append([num,detector,label,on,off])
    write_table(output/'konus_profiles.csv',['seq','start_s','counts'],profiles)
    write_table(output/'rxte_profiles.csv',['seq','start_s','all','hard','soft'],rxte)
    write_table(output/'events.csv',['seq','konus_id','trigger_utc_s','konus_min_keV',
                'konus_max_keV','rxte_ra','rxte_dec','rxte_mode'],metadata)
    write_table(output/'triplet_counts.csv',['seq','pcu','band','on_20s','off_100s'],flanks)
    (output/'inputs.json').write_text(json.dumps(inputs,indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cache', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    main(args.cache, args.output)
