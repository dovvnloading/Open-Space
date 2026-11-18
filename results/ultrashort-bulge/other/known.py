"""Make a published/control roster and compare against the frozen survey indexes.

Run at repository root. pdftotext is needed to extract published numerical tables.
Only identifiers, positions, selected physical parameters and receipts are exported.
"""
import csv
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from crossmatch import coords
from fetch import OUT, SHA, get, write_csv

CACHE = Path('.survey-cache/papers')
PAPERS = {
    'mroz17': '1707.07634', 'keck': '2303.04610', 'mroz19': '1811.00441',
    'kosh': '2303.08279', 'sumi': '2303.08280', 'gould22': '2204.03269',
    'gould23': '2306.04870', 'saturn': '2601.00057', 'park26': '2602.22709',
    'ryu26': '2606.27725', 'upper': '2611.00754',
}
KMT = 'https://kmtnet.kasi.re.kr/exoKMTNet/assets/data/planets.js?v=20261105b'
URLS = {k + '.pdf': 'https://arxiv.org/pdf/' + v for k, v in PAPERS.items()}
URLS['planets.js'] = KMT


def freeze():
    CACHE.mkdir(parents=True, exist_ok=True)
    receipts = []
    for name, url in URLS.items():
        path = CACHE / name
        if not path.exists():
            path.write_bytes(get(url))
        raw = path.read_bytes()
        receipts.append(dict(file=name, url=url, sha256=SHA(raw), input_bytes=len(raw),
                             retrieved_utc=datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).isoformat()))
        if name in ('mroz17.pdf', 'keck.pdf', 'kosh.pdf'):
            subprocess.run(['pdftotext', '-layout', path, path.with_suffix('.txt')], check=True)
    write_csv(OUT/'known_sources.csv', receipts)


def ref(key):
    return URLS[key+'.pdf']


def make_row(name, ra, dec, survey_class, source, aliases='', epoch=None, te=None, theta=None, year=None):
    if isinstance(ra, str):
        ra, dec = coords(ra, dec)
    if year is None:
        if epoch is None:
            year = int(name.split('-')[1])
        else:
            year = datetime.fromtimestamp((float(epoch)-2440587.5)*86400, timezone.utc).year
    return dict(known_id=name, aliases=aliases, season=year, ra_deg=round(ra, 8), dec_deg=round(dec, 8),
                t0_jd=epoch, te_d=te, theta_e_uas=theta, literature_class=survey_class, reference=source)


def literature():
    rows = []
    text = (CACHE/'mroz17.txt').read_text()
    alias = dict(re.findall(r'^(BLG[\d.]+)\s+(OGLE-\d+-BLG-\d+|\.\.\.)',
                            (CACHE/'keck.txt').read_text(), re.M))
    for name, ra, dec, t0, te in re.findall(
            r'^\s*(BLG\d+\.\d+\.\d+)\s+([\d:]+\.\d+)\s+(-[\d:]+\.\d+)\s+(\d{7}\.\d+)\s+([\d.]+)', text, re.M):
        rows.append(make_row(name, ra, dec, '2017 ultrashort PSPL candidate; 2024 Keck target',
                             ref('mroz17')+'|'+ref('keck'), alias[name].replace('...', ''), float(t0), float(te)))
    assert len(rows) == 6
    kosh = (CACHE/'kosh.txt').read_text()
    fits = {name: (2450000+float(t0), float(te)) for name, t0, te in re.findall(
        r'^\s*(MOA-9y-\d+)\s+[\d.]+\s+(\d{4}\.\d+)\s+([\d.]+)\s+[\d.]+\s+±', kosh, re.M)}
    kosh_short = kosh.split('Table 3. Short timescale event candidates with')[1].split('Note—The list of all')[0]
    for name, ra, dec, alias_id in re.findall(
            r'^\s*(MOA-9y-\d+)\s+gb[\w-]+\s+([\d:]+\.\d+)\s+(-[\d:]+\.\d+)\s+(\d{4}-BLG-\d+|–)',
            kosh_short, re.M):
        sub = 'sub-day candidate; no significant FS' if name not in ('MOA-9y-770','MOA-9y-5919') else 'sub-day FSPL candidate'
        rows.append(make_row(name, ra, dec, 'MOA nine-year '+sub, ref('kosh')+'|'+ref('sumi'),
                             ('MOA-'+alias_id if alias_id != '–' else ''), *fits[name]))
    assert len(rows) == 18
    # Official archive values are coordinates and ADS links to the underlying event papers.
    # Selected additional values are the authors' measurements, not fits to these alert reductions.
    archive = json.loads((CACHE/'planets.js').read_text().split('window.KMTNET_PLANETS = ')[1].split(';\n')[0])
    extra = {
        'OGLE-2017-BLG-0560': ('KMT-2017-BLG-2830', 2457859.52120, 0.901, 34.90),
        'OGLE-2016-BLG-1540': ('KMT-2016-BLG-2262', 2457606.724, 0.330, 8.80),
        'KMT-2024-BLG-0792': ('OGLE-2024-BLG-0516|Gaia24cdn', 2460434.323, 0.842, 18.6),
        'KMT-2024-BLG-0816': ('OGLE-2024-BLG-0519', 2460433.5733, 0.4461, None),
        'KMT-2024-BLG-3237': ('', 2460541.902, 0.54, 6.30),
    }
    gould_members = {'OGLE-2016-BLG-1540','KMT-2017-BLG-2820','OGLE-2019-BLG-0551','KMT-2019-BLG-2073'}
    for item in archive:
        if not item['is_ffp']:
            continue
        name = item['name']
        cat = 'KMT archive published FFP entry'
        source = item['ads_link']
        if name in gould_members:
            cat += '; Gould 2022 four-FFP subset'
        if name in gould_members or name=='OGLE-2017-BLG-0560':
            source += '|'+ref('gould22')
        if name=='KMT-2024-BLG-0792':
            cat += '; Gaia parallax mass'
            source += '|'+ref('saturn')
        if name=='KMT-2024-BLG-3237':
            source += '|'+ref('park26')
        rows.append(make_row(name,item['ra'],item['dec'],cat,source,*extra.get(name,())))
    rows.extend([
        make_row('OGLE-2012-BLG-1323','18:00:18.51','-28:35:01.7','published FSPL FFP candidate',
                 ref('mroz19'),'',2456161.107,0.155,2.37),
        make_row('KMT-2025-BLG-2093','17:50:56.94','-28:39:57.10','2026 published FSPL FFP candidate; position from KMT 2025 catalog',
                 ref('ryu26'),'',2460901.848,1.36,13.1),
        make_row('KMT-2026-BLG-0357','17:54:23.48','-28:14:00.24','isolated substellar lens above theta-E desert; comparison',
                 ref('upper'),'OGLE-2026-BLG-0156',2461131.56057,1.667,36.5),
        make_row('MOA-2021-BLG-145','17:56:19.30','-29:07:05.02','isolated substellar lens near upper edge of theta-E desert; comparison',
                 ref('upper'),'KMT-2021-BLG-0950',2459358.43302,1.1658,23.1),
        make_row('KMT-2022-BLG-2397','18:02:19.73','-26:53:28.10','published upper-shore brown-dwarf comparison',
                 ref('gould23'),'',2459812.31849,1.35,24.8),
    ])
    assert len({r['known_id'] for r in rows}) == len(rows)
    return sorted(rows,key=lambda r:(r['season'],r['known_id']))


def frozen():
    result = []
    def records(file, delimiter=','):
        with file.open() as handle:
            return list(csv.DictReader(handle, delimiter=delimiter))
    peaks = {r['event']:r['page_t0_jd'] for r in records(OUT/'photometry.csv')}
    for row in records(OUT/'index.csv'):
        name = row['event']
        result.append((name,*coords(row['ra'],row['dec']),peaks.get(name)))
    for p in sorted((OUT.parent/'kmt').glob('catalog-*.tsv')):
        for row in records(p,'\t'):
            if row.get('Event'):
                epoch = row['t_0']
                result.append((row['Event'],*coords(row['RA'],row['Dec']),float(epoch)+2400000 if epoch not in ('','-') else None))
    for row in records(OUT.parent/'kmt/event-pages-2026.tsv','\t'):
        result.append((row['event'],float('nan'),float('nan'),None))
    for row in records(OUT.parent/'ogle/catalog.csv'):
        result.append(('OGLE-'+row['event'],*coords(row['ra_hms'],row['dec_dms']),row['t0_hjd']))
    for row in records(OUT/'gaia.csv'):
        epoch = datetime.fromisoformat(row['trigger_tcb']).replace(tzinfo=timezone.utc).timestamp()/86400+2440587.5
        result.append((row['event'],float(row['ra_deg']),float(row['dec_deg']),epoch))
    return result


def match(roster, alerts):
    # A later season at the same location does not repeat an isolated lens passage.
    matches = []
    names = np.array([a[0] for a in alerts])
    ra, dec = np.deg2rad([[a[1],a[2]] for a in alerts]).T
    for src in roster:
        target_ra, target_dec = np.deg2rad([src['ra_deg'],src['dec_deg']])
        h = np.sin((dec-target_dec)/2)**2+np.cos(dec)*np.cos(target_dec)*np.sin((ra-target_ra)/2)**2
        sep = np.rad2deg(2*np.arcsin(np.sqrt(h)))*3600
        named = np.isin(names, [src['known_id'],*src['aliases'].split('|')])
        for i in np.flatnonzero((sep<=3)|named):
            alert = names[i]
            season = (2000+int(alert[4:6])) if alert.startswith('Gaia') else int(alert.split('-')[1])
            delta = None
            if src['t0_jd'] and alerts[i][3]:
                delta = round(float(alerts[i][3])-src['t0_jd'],5)
            if named[i]:
                tag = 'known by published ID'
            elif season != src['season']:
                tag = 'other season; sightline only'
            elif delta is not None and abs(delta)<=5:
                tag = 'same season; position and peak consistent'
            else:
                tag = 'same-season positional neighbour; identity unverified'
            matches.append(dict(known_id=src['known_id'],alert_id=alert,literature_season=src['season'],
                                alert_season=season,sep_arcsec=round(sep[i],5) if np.isfinite(sep[i]) else '',
                                delta_alert_peak_or_trigger_d=delta,status=tag))
    return matches


if __name__ == '__main__':
    freeze()
    known = literature()
    write_csv(OUT/'known_roster.csv',known)
    write_csv(OUT/'known_matches.csv',match(known,frozen()))
