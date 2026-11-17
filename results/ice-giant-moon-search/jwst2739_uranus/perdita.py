"""Test masked, trailed photometry along the known moon's catalog orbit rotated in longitude."""
import argparse
import csv
import sys
from pathlib import Path
import importlib.util
import numpy as np

from screen import setup,search
spec=importlib.util.spec_from_file_location('january_perdita',Path(__file__).resolve().parents[1]/'jwst2768/perdita.py')
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)


def run(cube_file,ee,output,body,labels):
    cube,cuts,maps=setup(cube_file,ee)
    if int(cube['unmasked_control_body'])!=body:
        raise ValueError('control moon differs from supplied cube')
    angles=np.arange(-180.,180.,.25)
    x,y=old.rotate_tracks(cube,angles,body=body)
    n,w=search.sample_paths(maps,x,y,cuts)
    z=search.score(n,w,cube['visit'],cuts)
    with output.open('w',newline='') as f:
        wr=csv.writer(f);wr.writerow(['offset_deg',*[f'{label}_{field}' for field in ['snr','dropout','n'] for label in labels],'rank'])
        for i,a in enumerate(angles):
            wr.writerow([a,*z['visit_snr'][i],*z['dropout'][i],*z['n'][i],z['rank'][i]])
    for j,visit in enumerate(np.unique(cube['visit'])):
        allow=(np.abs(angles)<=15)&(z['n'][:,j]>=2)
        if allow.any():
            i=np.argmax(np.where(allow,z['visit_snr'][:,j],-np.inf))
            print(visit,angles[i],z['visit_snr'][i,j],z['dropout'][i,j],z['n'][i,j],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['cube','ee','out']:p.add_argument(k,type=Path)
    p.add_argument('--body',type=int,default=725)
    p.add_argument('--labels',nargs=2,default=['feb','sep'])
    a=p.parse_args();run(a.cube,a.ee,a.out,a.body,a.labels)
