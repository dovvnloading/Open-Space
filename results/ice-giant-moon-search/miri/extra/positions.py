"""Freeze spacecraft-relative mid-exposure positions for the MIRI coverage and search."""
import argparse,csv,hashlib,json
from pathlib import Path
import numpy as np
from fetch import ROOT,products
from geometry import PLANET,MOONS,query
from miri_audit import write_csv

def run(cache):
    rows=[];receipts=[]
    for planet in ('Uranus','Neptune'):
        selected=[f for f,_ in products() if f['planet']==planet]
        for body in [PLANET[planet],*[v for v in MOONS[planet] if v!=75052]]:
            use=selected if body==PLANET[planet] else [f for f in selected if planet=='Uranus' or f['planet_on_chip']=='1']
            jd=np.array([float(f['exp_mid_mjd_utc']) for f in use])+2400000.5
            st=query(body,'JWST',jd,cache)
            exact=np.unique(jd.round(9))
            wanted=[tuple(exact[i:i+60]) for i in range(0,len(exact),60)];found={}
            for file in cache.glob(str(body)+'-*.json'):
                raw=file.read_bytes();result=json.loads(raw).get('result','')
                if '$$SOE' not in result or '$$EOE' not in result:continue
                table=list(csv.reader(result.split('$$SOE')[1].split('$$EOE')[0].strip().splitlines()))
                times=tuple(round(float(r[0]),9) for r in table)
                if times not in wanted:continue
                target=next(line.strip() for line in result.splitlines() if line.startswith('Target body name:'))
                if f'({body})' not in target:continue
                found[times]=dict(body=body,epochs=len(times),first_jd=times[0],last_jd=times[-1],target=target,
                                  url_sha256=file.stem.split('-',1)[1],response_sha256=hashlib.sha256(raw).hexdigest())
            if len(found)!=len(wanted):raise ValueError(f'{body}: missing response receipts')
            receipts.extend(found[k] for k in wanted)
            for f,vec in zip(use,st):
                rows.append(dict(planet=planet,body=body,file=f['file'],ra_deg=vec[0],dec_deg=vec[1],diameter_arcsec=vec[2],r_au=vec[3],delta_au=vec[4]))
    write_csv(Path(__file__).parent/'positions.csv',list(rows[0]),rows)
    write_csv(Path(__file__).parent/'horizons_receipts.csv',list(receipts[0]),receipts)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('cache',type=Path);run(p.parse_args().cache)
