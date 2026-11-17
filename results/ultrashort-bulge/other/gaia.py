"""Inventory the main Gaia Science Alerts stream in |l|<10, |b|<7 deg, 2024 onward.

Keep every classification. Alert dates are Gaia's triggering observations; the 3-column
Gaia light curves contain TCB JD and G mag with no epoch uncertainty. The flux error
proxy is a robust off-alert per-source magnitude scatter, floored at 0.02 mag. It
supports screening only. Independent ground-based confirmation is essential.
"""
import csv
import io
import json
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from fetch import get, write_csv, SHA, OUT

CACHE = Path('.survey-cache')
URL = 'https://gsaweb.ast.cam.ac.uk/alerts/alerts.csv'
MATRIX = np.array([[-.0548755604,-.8734370902,-.4838350155],
                   [.4941094279,-.4448296300,.7469822445],
                   [-.8676661490,-.1980763734,.4559837762]])


def screen_index():
    idx = CACHE/'gaiaidx'
    if not idx.exists():
        idx.write_bytes(get(URL))
    raw=idx.read_bytes()
    rows=list(csv.DictReader(io.StringIO(raw.decode('utf-8-sig')), skipinitialspace=True))
    selected=[]
    for r in rows:
        if r['Date'][:4] not in ('2024','2025','2026'):
            continue
        a,d = np.deg2rad([float(r['RaDeg']), float(r['DecDeg'])])
        x,y,z = MATRIX @ [np.cos(a)*np.cos(d), np.sin(a)*np.cos(d), np.sin(d)]
        lon,lat = np.rad2deg([np.arctan2(y,x),np.arcsin(z)])
        if abs(lon) < 10 and abs(lat) < 7:
            selected.append(dict(event=r['#Name'],trigger_tcb=r['Date'],ra_deg=r['RaDeg'],
                                 dec_deg=r['DecDeg'],l_deg=round(lon,7),b_deg=round(lat,7),
                                 classification=r['Class'],comment=r['Comment']))
    (OUT/'gaia_index_receipt.json').write_text(json.dumps(dict(url=URL,bytes=len(raw),sha256=SHA(raw),
        retrieved_utc=datetime.fromtimestamp(idx.stat().st_mtime,timezone.utc).isoformat(),
        all_rows=len(rows),latest_observed_tcb=max(r['Date'] for r in rows),
        selected={} if not selected else {y:sum(r['trigger_tcb'][:4]==y for r in selected)
                                          for y in ['2024','2025','2026']}),indent=2)+'\n')
    return selected


def curve(r):
    event = r['event']
    url = f'https://gsaweb.ast.cam.ac.uk/alerts/alert/{event}/lightcurve.csv'
    r.update(url=url)
    try:
        path = CACHE/'gaiaraw'/f'{event}.csv'
        if not path.exists():
            path.write_bytes(get(url))
        raw=path.read_bytes()
        r.update(sha256=SHA(raw),input_bytes=len(raw))
        lines=raw.decode().splitlines()
        if lines[0].strip()!=event:
            raise ValueError('wrong name in file')
        rows = [e for e in csv.reader(lines[2:]) if e]
        phot=np.array([(float(e[1]),float(e[2])) for e in rows if e[2].strip().replace('.','',1).isdigit()]).reshape((-1,2))
        trigger=(datetime.fromisoformat(r['trigger_tcb']).replace(tzinfo=timezone.utc).timestamp()/86400)+2440587.5
        r.update(rows=len(rows),usable=len(phot))
        if len(phot)<3:
            raise ValueError('fewer than 3 numeric measurements')
        t,mag=phot.T
        bg=mag[abs(t-trigger)>30]
        if len(bg)<6:
            bg=mag
        noise=max(.02,1.4826*np.median(abs(bg-np.median(bg))))
        flux=10**(-.4*(mag-np.median(bg)))
        arr=np.column_stack([t,flux,flux*(np.log(10)/2.5)*noise])
        np.savez_compressed(CACHE/'curves'/f'{event}.npz',data=arr)
        record=dict(event=event,page_t0_jd=trigger,page_te_d=None,status='ok',time_scale='JD(TCB)',difference=False)
        (CACHE/'records'/f'{event}.json').write_text(json.dumps(record)+'\n')
        absdt=abs(t-trigger)
        r.update(rows=len(rows),usable=len(t),sha256=SHA(raw),input_bytes=len(raw),
                 first_jd_tcb=t.min(),last_jd_tcb=t.max(),n_within_0p5d=int((absdt<.5).sum()),
                 scatter_proxy_mag=noise,status='ok')
    except Exception as e:
        r.update(status='error',error=f'{type(e).__name__}: {e}')
    return r


if __name__=='__main__':
    for sub in ['gaiaraw','curves','records']:
        (CACHE/sub).mkdir(parents=True,exist_ok=True)
    with ThreadPoolExecutor(8) as pool:
        items=list(pool.map(curve,screen_index()))
    write_csv(OUT/'gaia.csv', items)
