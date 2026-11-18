"""Audit the union of tE<0.5d OGLE grid winners, with coverage at each model's own peak."""
import argparse
import csv
import gzip
import json
from pathlib import Path
import sys

import numpy as np

from counterparts import HERE
sys.path.insert(0,str(HERE.parents[1]/'ogle-dark-lenses'))
import lens


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('cache',type=Path);p.add_argument('grid',type=Path,help='checkpoint dir, CSV or CSV.gz')
    a=p.parse_args()
    cat={r['event']:r for r in csv.DictReader((HERE/'catalog.csv').open())}
    if a.grid.is_dir():
        rows=[json.loads(f.read_text()) for f in a.grid.glob('OGLE-*.json')]
    else:
        opener=gzip.open if a.grid.suffix=='.gz' else open
        with opener(a.grid,'rt') as f: rows=[r for r in csv.DictReader(f) if r['event'].startswith('OGLE-')]
    if {r['event'][5:] for r in rows}!=cat.keys() or len(rows)!=len(cat) or any(r.get('status',r.get('fit_status'))!='ok' for r in rows):
        raise ValueError('grid does not contain exactly one successful fit for each OGLE bulge ID')
    out=[]
    for r in sorted(rows,key=lambda z:z['event']):
        if min(float(r['ps_te']),float(r['fs_te']))>=.5:continue
        event=r['event'][5:];y,_,n=event.split('-');ts=lens.read_phot(a.cache/y/f'blg-{n}.dat')[0]+2450000
        record=dict(event=event,ews_te=float(cat[event]['tau_days']),ps_te=float(r['ps_te']),fs_te=float(r['fs_te']),
                    fs_rho=float(r['fs_rho']),fs_minus_ps_gain=round(float(r['fs_gain'])-float(r['ps_gain']),3))
        for fam in ['ps','fs']:
            peak=float(r[fam+'_t0']);delta=ts-peak
            record[fam+'_hjd_minus_ews']=round(peak-float(cat[event]['t0_hjd']),5)
            # Bins exclude a central ±2min interval to avoid mistaking a single exact
            # peak epoch or numerical rounding at that epoch for an actual wing.
            record[fam+'_before']=int(np.sum((delta>=-.5)&(delta< -2/1440)))
            record[fam+'_center']=int(np.sum(abs(delta)<=2/1440))
            record[fam+'_after']=int(np.sum((delta<=.5)&(delta>2/1440)))
        out.append(record)
    with (HERE/'short-grid-audit.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,out[0],lineterminator='\n');w.writeheader();w.writerows(out)

if __name__=='__main__':main()
