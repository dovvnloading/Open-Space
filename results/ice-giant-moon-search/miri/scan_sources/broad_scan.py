"""Exact exposure-integrated circular-orbit survey in radial batches."""
import argparse
import csv
import json
import tempfile
from pathlib import Path
import numpy as np
from geometry import generate
from search import Cuts, file_hash, load_npz, prepare, sample_paths, score, validate_cube, validate_tracks


def run(cube, output, lo, hi, da, dphase, support):
    out=Path(output)
    out.mkdir(parents=True,exist_ok=True)
    with np.load(cube) as f:
        if str(f["unmasked_control_body"])!="None":
            raise ValueError("blind scan requires known-moon masks")
    c=load_npz(cube)
    kernels=validate_cube(c)
    cuts=Cuts(psf_support=support)
    maps=prepare(c,kernels,cuts)
    epoch=None
    progress=[]
    # Periodic radii avoid losing work if interrupted. Each output contains all
    # sky phases, including screened-out trials and their exposure counts.
    with tempfile.TemporaryDirectory() as tmp:
        for a in np.arange(lo,hi+.000001,da):
            name=out/f"r{a:g}.npz"
            if not name.exists():
                tracks_file=Path(tmp)/"tracks.npz"
                generate(cube,tracks_file,grid=(a,a,da,dphase),epoch=epoch)
                tracks=load_npz(tracks_file)
                validate_tracks(tracks,c)
                epoch=float(tracks['epoch_jd_tdb'])
                n,w=sample_paths(maps,tracks['x_path'],tracks['y_path'],cuts)
                pos=score(n,w,c['visit'],cuts)
                neg=score(-n,w,c['visit'],cuts)
                np.savez_compressed(name,phase_deg=tracks['phase_deg'],negative_rank=neg['rank'],
                                    epoch_jd_tdb=epoch,cube_sha256=file_hash(cube),support=support,**pos)
            with np.load(name) as data:
                if (str(data['cube_sha256'])!=file_hash(cube) or float(data['support'])!=support
                        or not np.array_equal(data['phase_deg'],np.arange(0,360,dphase))):
                    raise ValueError("refusing to mix incompatible scan checkpoints")
                epoch=float(data['epoch_jd_tdb'])
                rank,negative=data['rank'],data['negative_rank']
                progress.append([a,len(rank),int(data['eligible'].sum()),int((rank>=1).sum()),
                                 int((negative>=1).sum()),float(max(rank)),float(max(negative))])
            if len(progress)%10==0:
                print(a,sum(r[2] for r in progress),flush=True)
    with (out/'radii.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['a_km','trials','eligible','passing','negative_passing','max_rank','max_negative_rank']);w.writerows(progress)
    with (out/'passing.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['a_km','phase_deg','pooled_score','visit1_score','visit2_score','dropout1','dropout2','fit_jy','n_visit1','n_visit2'])
        for row in progress:
            with np.load(out/f"r{row[0]:g}.npz") as d:
                for i in np.flatnonzero(d['rank']>=1):
                    w.writerow([row[0],d['phase_deg'][i],d['snr'][i],*d['visit_snr'][i],*d['dropout'][i],d['flux'][i],*d['n'][i]])
    print(json.dumps(dict(totals=np.sum(np.array(progress)[:,1:5],axis=0).astype(int).tolist(),epoch=epoch)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('cube')
    p.add_argument('output')
    p.add_argument('--grid',nargs=4,type=float,required=True,metavar=('RMIN','RMAX','DR_KM','DPHASE_DEG'))
    p.add_argument('--support',type=float,default=.90)
    a=p.parse_args()
    run(a.cube,a.output,*a.grid,a.support)
