"""Inject one off-grid source at a time and test matching orbits in the published search bank."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from screen import setup,projector,search,subtract_ring_profile


def run(cube_file,ee,out,levels,n):
    cube,cuts,maps=setup(cube_file,ee)
    kernels=search.validate_cube(cube)
    epoch,project=projector(cube)
    mid=cube['path_mjd'].shape[1]//2
    rng=np.random.default_rng(27392023)
    r0=rng.uniform(79500,80500,120);p0=rng.uniform(0,360,120)
    x0,y0=project(r0,p0,slice(None))
    a,b=search.sample_paths(maps,x0,y0,cuts)
    supported=np.flatnonzero(search.score(a,b,cube['visit'],cuts)['eligible'])
    ids=rng.choice(supported,n,replace=False)
    rows=[]
    print('supported',len(supported),'selected',n,flush=True)
    for k in ids:
        # At most 3 degrees from the original first-visit phase can lie inside
        # the 1.5-pixel association circle at either scale in the tested radial band.
        angles=np.arange(0,360.,.75)
        use=abs((angles-p0[k]+180)%360-180)<=3
        ra,ph=[v.ravel() for v in np.meshgrid(np.arange(79500.,80500.01,.5),angles[use],indexing='ij')]
        x,y=project(ra,ph,[mid])
        match=np.hypot(x[:,:,0]-x0[k,:,mid],y[:,:,0]-y0[k,:,mid]).max(1)<=1.5
        ra,ph,x,y=ra[match],ph[match],x[match,:,0],y[match,:,0]
        model=np.zeros_like(cube['image'])
        for j,kernel in enumerate(kernels):
            stamp=search.clipped_stamp(model[j].shape,*search.trail_kernel(x0[k,j],y0[k,j],kernel))
            if stamp is not None:
                sl,ps=stamp; model[j][sl]+=ps
        for flux in levels:
            planted={**cube,'image':cube['image']+flux*model}
            pm=search.prepare(subtract_ring_profile(planted),kernels,cuts)
            n1,w1=search.sample(pm,x,y)
            ids1=np.flatnonzero(search.score(n1,w1,cube['visit'],cuts)['rank']>=.8)
            hit=False
            # Individual near-truth templates suffice: the blind search emits every pass.
            # Stop at the first matching pass; hits outside the association cannot affect it.
            for ix in np.array_split(ids1,max(1,(len(ids1)+15)//16)):
                if not len(ix):continue
                xp,yp=project(ra[ix],ph[ix],slice(None))
                nn,ww=search.sample_paths(pm,xp,yp,cuts)
                if (search.score(nn,ww,cube['visit'],cuts)['rank']>=1).any():
                    hit=True;break
            rows.append(dict(flux=flux,plant_index=int(k),a_km=r0[k],phase_deg=p0[k],recovered=bool(hit),
                             status='recovered' if hit else 'missed'))
    with out.with_suffix('.csv').open('w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=rows[0]);wr.writeheader();wr.writerows(rows)
    report=dict(random_orbits=120,observable=len(supported),n_at_each_flux=n,seed=27392023,epoch=epoch,
                injection='Trailed model in calibrated science images before elliptical-ring and glare subtraction, fixed noise and DQ.',
                association_px=1.5,summary=search.summarize_plants(rows),cube_sha256=search.file_hash(cube_file),ee_sha256=search.file_hash(ee))
    out.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
    print([(r['flux'],r['recovered']) for r in report['summary']],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['cube','ee','out']:p.add_argument(k,type=Path)
    p.add_argument('--flux',type=float,nargs='+',required=True);p.add_argument('--n',type=int,default=18)
    a=p.parse_args();run(a.cube,a.ee,a.out,a.flux,a.n)
