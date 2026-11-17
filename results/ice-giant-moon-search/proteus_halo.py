"""Measure co-moving residuals around Proteus and veto the contaminated region."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from search import Cuts, file_hash, load_npz, prepare, sample, snr, validate_cube


def run(paths,outputs,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    cubes=list(map(load_npz,paths))
    profiles=[];last=0.
    for c in cubes:
        cut=Cuts(psf_support=.8)
        maps=prepare(c,validate_cube(c,min_visits=1),cut)
        xy=c['known_xy'][:,list(c['known_id']).index(808)]
        angle=np.linspace(0,2*np.pi,72,endpoint=False)[:,None]
        for rad in np.arange(6.,41.):
            x=xy[None,:,0]+rad*np.cos(angle);y=xy[None,:,1]+rad*np.sin(angle)
            n,w=sample(maps,x,y);usable=np.sum(w>0,axis=1)>=3
            if not np.any(usable):continue
            scores=snr(n.sum(1),w.sum(1))[usable]
            if scores.max()>=5:
                last=max(last,rad)
            profiles.append([str(c['band']),rad,int(usable.sum()),*np.percentile(scores,[50,75,90,100])])
    radius=float(last+3)
    with (out/'proteus_profile.csv').open('w',newline='') as f:
        wr=csv.writer(f);wr.writerow(['band','radius_px','supported_position_angles','median_score','p75_score','p90_score','max_score']);wr.writerows(profiles)
    for c,dst in zip(cubes,outputs):
        p=c['known_path'][list(c['known_id']).index(808)]
        yy,xx=np.indices(c['image'].shape[1:])
        for j,trail in enumerate(p):
            masked=np.logical_or.reduce([(xx-x)**2+(yy-y)**2<=radius**2 for x,y in trail])
            c['valid'][j,masked]=False;c['geometry_valid'][j,masked]=False
        c['proteus_halo_px']=radius
        np.savez_compressed(dst,**c)
    (out/'mask.json').write_text(json.dumps(dict(input_sha256=[file_hash(p) for p in paths],
        cut_radius_px=radius,last_5sigma_radius_px=last,output_sha256=[file_hash(p) for p in outputs],
        guard_px=3,pa_samples=72,radial_probe_px=[6,40],required_exposures=3),indent=2)+'\n')
    print('Proteus mask radius',radius)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('f140','f210','out140','out210','out'):
        p.add_argument(key)
    a=p.parse_args();run([a.f140,a.f210],[a.out140,a.out210],a.out)
