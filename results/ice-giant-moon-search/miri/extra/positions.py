"""Freeze spacecraft-relative mid-exposure positions for the MIRI coverage and search."""
import argparse,csv,hashlib
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
            for f,vec in zip(use,st):
                rows.append(dict(planet=planet,body=body,file=f['file'],ra_deg=vec[0],dec_deg=vec[1],diameter_arcsec=vec[2],r_au=vec[3],delta_au=vec[4]))
    write_csv(Path(__file__).parent/'positions.csv',list(rows[0]),rows)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('cache',type=Path);run(p.parse_args().cache)
