"""Retain local positive and inverted-score orbital peaks separately in one deep visit."""
import argparse
import csv
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import numpy as np
from scipy import ndimage
base=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(base/'jwst9482'))
from link_visits import make_cube,evaluate,deep_score,geometry,search


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('cube',type=Path)
    p.add_argument('ee',type=Path)
    p.add_argument('out',type=Path)
    p.add_argument('--grid',type=float,nargs=4,default=[54500,86000,250,.75])
    a=p.parse_args()
    cube=make_cube(a.cube,a.ee,True)
    psfs=search.validate_cube(cube,min_visits=1)
    cuts=search.Cuts(psf_support=.8)
    maps=search.prepare(cube,psfs,cuts)
    lo,hi,dr,dp=a.grid
    radii=np.arange(lo,hi+dr*1e-9,dr)
    phases=np.arange(0,360,dp)
    shape=(len(radii),len(phases))
    fields={k:np.zeros(shape) for k in ['score','rank','neg_rank','eligible','early','late']}
    a.out.mkdir(parents=True,exist_ok=False)
    with TemporaryDirectory() as tmp:
        name=Path(tmp)/'b.npz'
        for start in range(0,len(radii),3):
            r=radii[start:start+3]
            geometry.generate(a.cube,name,grid=[r[0],r[-1],dr,dp])
            bank=search.load_npz(name);search.validate_tracks(bank,cube)
            num,info=evaluate(cube,bank,psfs,maps,cuts)
            good,z,e,l,elig=deep_score(num,info)
            neg=deep_score(-num,info)[0]
            for k,v in [('rank',good),('neg_rank',neg),('score',z),('early',e),('late',l),('eligible',elig)]:
                # geometry.generate uses a radius-fast meshgrid within each phase.
                fields[k][start:start+len(r)] = v.reshape(len(phases),len(r)).T
            epoch=float(bank['epoch_jd_tdb'])
    with open(a.out/'peaks.csv','w',newline='') as f:
        writer=csv.writer(f)
        writer.writerow(['sign','a_km','phase_deg','score','rank','early','late'])
        for sign,key in [(1,'rank'),(-1,'neg_rank')]:
            # Different local maxima within the same object can survive this filter.
            data=np.where(fields[key]>=1,sign*fields['score'],-1e9)
            maxima=ndimage.maximum_filter(data,(5,5),mode=('constant','wrap'),cval=-1e9)
            for i,j in np.argwhere((data>=7)&(data==maxima)):
                writer.writerow([sign,radii[i],phases[j],sign*fields['score'][i,j],fields[key][i,j],
                                 sign*fields['early'][i,j],sign*fields['late'][i,j]])
    summary=dict(hash=search.file_hash(a.cube),ee_sha256=search.file_hash(a.ee),script_sha256=search.file_hash(Path(__file__)),
        epoch_jd_tdb=epoch,radii=a.grid,templates=int(np.prod(shape)),eligible=int(fields['eligible'].sum()),
        screen_passing_signs={key:int((fields[key]>=1).sum()) for key in ['rank','neg_rank']},
        local_maximum_window=[5,5],minimum_retained_absolute_score=7)
    (a.out/'parameters.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__=='__main__':
    main()
