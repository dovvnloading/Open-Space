"""Find independent KMT/OGLE cross-identifications across shared survey seasons."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re

import numpy as np
from scipy.spatial import cKDTree

from counterparts import HERE,coords


def unit(ra,dec):
    a,d=np.radians(coords(ra,dec))
    return [np.cos(a)*np.cos(d),np.sin(a)*np.cos(d),np.sin(d)]


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('cache',type=Path);cache=ap.parse_args().cache
    out=[];receipts=[]
    for yr in ['2024','2025','2026']:
        path=cache/yr/'lenses.par'
        if not path.exists():
            import urllib.request
            path.parent.mkdir(parents=True,exist_ok=True)
            data=urllib.request.urlopen(f'https://www.astrouw.edu.pl/ogle/ogle4/ews/{yr}/lenses.par',timeout=90).read()
            path.write_bytes(data)
        data=path.read_bytes()
        receipts.append(dict(year=yr,url=f'https://www.astrouw.edu.pl/ogle/ogle4/ews/{yr}/lenses.par',
                             sha256=hashlib.sha256(data).hexdigest(),bytes=len(data)))
        og={v[0]:v for l in data.decode().splitlines() if len((v:=l.split()))>=5 and '-BLG-' in v[0]}
        ids=list(og);uo=np.array([unit(og[k][3],og[k][4]) for k in ids]);tree=cKDTree(uo)
        kmt_path=HERE.parent/'kmt'/f'catalog-{yr}.tsv'
        if kmt_path.exists():
            entries=list(csv.DictReader(kmt_path.open(),delimiter='\t'))
        else:
            from html import unescape
            import urllib.request
            u=f'https://kmtnet.kasi.re.kr/ulens/event/{yr}/'
            hp=cache/yr/'kmt-index.html'
            if not hp.exists():hp.write_bytes(urllib.request.urlopen(u,timeout=90).read())
            html=hp.read_bytes()
            receipts.append(dict(year=yr,url=u,sha256=hashlib.sha256(html).hexdigest(),bytes=len(html)))
            cells=[[unescape(re.sub(r'<[^>]*>', '', x)).strip() for x in re.findall(r'<td[^>]*>(.*?)</td>', tr,re.S|re.I)]
                   for tr in re.findall(r'<tr[^>]*>(.*?)</tr>',html.decode(),re.S|re.I)]
            header=next(r for r in cells if r and r[0]=='Event')
            entries=[dict(zip(header,r)) for r in cells if r and r[0].startswith('KMT-')]
        for k in entries:
            try: uk=unit(k['RA'],k['Dec'])
            except (TypeError,ValueError):continue
            m=re.findall(r'OB(\d{2})(\d{4})', k['Related event'])
            official={'20'+year+'-BLG-'+no for year,no in m}
            near={ids[i] for i in tree.query_ball_point(uk,2/206264.806)}
            for n in sorted(official|near):
                # Preserve explicit aliases even with poor astrometry or unconstrained t0.
                a=og.get(n)
                if a is None:continue
                dr=2*np.arcsin(min(1,np.linalg.norm(np.asarray(uk)-unit(a[3],a[4]))/2))*206264.806
                try:dt=float(k['t_0'])+2400000-float(a[5])
                except (IndexError,TypeError,ValueError):dt=None
                official_id=n in official
                if not official_id and (dt is None or abs(dt)>10):continue
                out.append(dict(kmt_id=k['Event'],ogle_id='OGLE-'+n,official=int(official_id),sep_arcsec=round(dr,3) if yr!='2026' else '',
                                delta_peak_d=round(dt,4) if dt is not None and yr!='2026' else '',
                                kmt_te_d=k['t_E'] if yr!='2026' else '',
                                ogle_te_d=a[7] if len(a)>7 else ''))
    with (HERE/'kmt-ogle-aliases.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,out[0].keys(),lineterminator='\n');w.writeheader();w.writerows(out)
    (HERE/'alias-sources.json').write_text(json.dumps(receipts,indent=2)+'\n')

if __name__=='__main__':main()
