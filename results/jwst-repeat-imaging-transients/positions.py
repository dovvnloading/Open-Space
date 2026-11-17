"""Merge morphology-selected coordinates across filters and epoch pairs."""
import argparse
from pathlib import Path
import pandas as pd
import numpy as np

def merge(paths):
    result=[];links=[]
    for path in paths:
        a=pd.read_csv(path/'detections.csv')
        for r in a.query('screen_pass == 1').itertuples():
            sep=[3600*np.hypot((o['ra_deg']-r.ra_deg)*np.cos(np.deg2rad(r.dec_deg)),o['dec_deg']-r.dec_deg) for o in result]
            if sep and min(sep)<.14:
                name=result[np.argmin(sep)]['id']
            else:
                name=f'N{len(result)+1:03d}'
                result.append(dict(id=name,ra_deg=r.ra_deg,dec_deg=r.dec_deg))
            links.append(dict(id=name,run=path.name,x=r.x,y=r.y,snr=r.snr,dflux_njy=r.dflux_njy,err_njy=r.err_njy))
    return pd.DataFrame(result),pd.DataFrame(links)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('out',type=Path);p.add_argument('runs',nargs='+',type=Path);a=p.parse_args()
    a.out.mkdir(exist_ok=True,parents=True)
    pos,links=merge(a.runs);pos.to_csv(a.out/'positions.csv',index=False);links.to_csv(a.out/'links.csv',index=False)
