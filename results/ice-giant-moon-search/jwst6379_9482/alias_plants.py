"""Blind two-epoch injection recovery with orbital cycle-count ambiguity included."""
import argparse
import csv
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import numpy as np
from scipy import ndimage
here=Path(__file__).resolve().parent
sys.path.insert(0,str(here))
from link_aliases import compatible
from link_visits import remove_stationary,radial_psf,geometry,search,trials,deep_score


def detect(maps, cache, grid, cube, cuts, radii, phases):
    n,i=cache.measure(maps,len(grid['x']),cuts)
    rank,z,*_=deep_score(n,i)
    score=np.where(rank>=1,z,-1e9).reshape(len(phases),len(radii)).T
    peaks=ndimage.maximum_filter(score,(5,5),mode=('constant','wrap'),cval=-1e9)
    return [dict(a_km=radii[r],phase_deg=phases[p],index=p*len(radii)+r)
            for r,p in np.argwhere((score>=7)&(peaks==score))]


def near_peaks(peaks,grid,px,py):
    out=[]
    for pk in peaks:
        j=pk['index']
        delta=np.hypot(grid['x'][j]-px,grid['y'][j]-py)
        if np.max(delta)<=4:
            out.append(pk)
    return out


def qualifying(sets, truth, grids, paths, cubes, psfs, maps, cut, epoch, dt, tmp):
    # Geometric restriction only accelerates evaluation of output-to-truth associations.
    # Both peak maps were computed by the entire local blind bank, with local maxima cuts.
    left,right=[near_peaks(sets[j],grids[j],truth[j][0],truth[j][1]) for j in range(2)]
    vals=set()
    for p1 in left:
        for p2 in right:
            for r,ph1,ph2,_ in compatible(p1,p2,dt,.75,375):
                vals.add((r,ph1))
    if not vals:return None
    elem=np.zeros((len(vals),6))
    elem[:,0],elem[:,-1]=np.asarray(sorted(vals)).T
    keep=np.ones(len(elem),bool)
    for j in range(2):
        b=trials(tmp/'aliases.csv',paths[j],tmp/'aliases.npz',elem,epoch)
        for px,key in zip(truth[j],['x','y'],strict=True):
            keep &= np.max(abs(b[key]-px),axis=1)<=1.5
        keep &= np.max(np.hypot(b['x']-truth[j][0],b['y']-truth[j][1]),axis=1)<=1.5
        if not keep.any():return None
        ind=np.flatnonzero(keep)
        n,i=search.sample_paths(maps[j],b['x_path'][ind],b['y_path'][ind],cut)
        rr,zz,*_=deep_score(n,i)
        keep[ind]=(rr>=1)&(zz>=7)
    indices=np.flatnonzero(keep)
    if len(indices):return tuple(elem[indices[0]][[0,5]])
    return None


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['early','later','ee','out']:
        p.add_argument(key,type=Path)
    p.add_argument('--flux',type=float,nargs='+',default=[5e-8,1e-7,3e-7])
    p.add_argument('--n',type=int,default=36)
    p.add_argument('--clear-dispositions',type=Path,help='replay original planted samples only when the independent exclusion passed')
    p.add_argument('--seed',type=int,default=16027)
    args=p.parse_args();paths=[args.early,args.later]
    clear=None
    if args.clear_dispositions:
        clear={(float(r['flux']),int(r['plant_index'])) for r in csv.DictReader(args.clear_dispositions.open()) if r['guaranteed_clear']=='True'}
    rng=np.random.default_rng(args.seed)
    crops=[search.load_npz(path) for path in paths]
    for c in crops:
        ps=radial_psf(args.ee,float(c['scale_arcsec']),'F150W2')
        c.update(psf=ps,psf_fraction=ps.sum())
    cubes=[remove_stationary(c) for c in crops]
    kernels=[search.validate_cube(c,min_visits=1) for c in cubes]
    cuts=search.Cuts(psf_support=.8)
    baseline=[search.prepare(c,k,cuts) for c,k in zip(cubes,kernels,strict=True)]
    radii=np.arange(79500,80501,250);phases=np.arange(0,360,.75)
    args.out.mkdir(parents=True,exist_ok=False)
    with TemporaryDirectory() as tmpdir:
        tmp=Path(tmpdir)
        grids=[]
        for path in paths:
            geometry.generate(path,tmp/'g.npz',grid=[radii[0],radii[-1],250,.75])
            grids.append(search.load_npz(tmp/'g.npz'))
        epoch,epoch2=[float(g['epoch_jd_tdb']) for g in grids]; dt=(epoch2-epoch)*86400
        caches=[search.TemplateCache(c,g,k) for c,g,k in zip(cubes,grids,kernels,strict=True)]
        sets=[detect(b,z,g,c,cuts,radii,phases) for b,z,g,c in zip(baseline,caches,grids,cubes,strict=True)]
        elem=np.zeros((120,6));elem[:,0]=rng.uniform(radii[0],radii[-1],len(elem));elem[:,-1]=rng.uniform(0,360,len(elem))
        banks=[trials(tmp/'plant.csv',path,tmp/'p.npz',elem,epoch) for path in paths]
        use=np.ones(len(elem),bool)
        for c,b,k,m in zip(cubes,banks,kernels,baseline,strict=True):
            n,w=search.sample_paths(m,b['x_path'],b['y_path'],cuts)
            use &= deep_score(np.zeros_like(n),w)[-1]
        eligible=np.flatnonzero(use)
        if len(eligible)<args.n:raise ValueError(f'only {len(eligible)} supported trial orbits')
        saved=[]
        with open(args.out/'plants.csv','x',newline='') as ff:
            fields=['flux','plant_index','a_km','phase_epoch1','recovered','status','recovered_a_km','recovered_phase_deg']
            writer=csv.DictWriter(ff,fieldnames=fields);writer.writeheader()
            for flux in args.flux:
                for i in rng.choice(eligible,args.n,replace=False):
                    if clear is not None and (flux,int(i)) not in clear:continue
                    truth=[(b['x'][i],b['y'][i]) for b in banks]
                    row=dict(flux=flux,plant_index=int(i),a_km=elem[i,0],phase_epoch1=elem[i,-1],recovered=False,status='missed')
                    if qualifying(sets,truth,grids,paths,cubes,kernels,baseline,cuts,epoch,dt,tmp) is not None:
                        row['status']='preexisting_peak'
                    else:
                        maps=[]
                        for raw,kernel,bank in zip(crops,kernels,banks,strict=True):
                            img=raw['image'].copy()
                            for j,ps in enumerate(kernel):
                                patch=search.clipped_stamp(img[j].shape,*search.trail_kernel(bank['x_path'][i,j],bank['y_path'][i,j],ps))
                                if patch is not None:
                                    cut,stamp=patch;img[j][cut]+=stamp*flux
                            c=remove_stationary({**raw,'image':img})
                            maps.append(search.prepare(c,kernel,cuts))
                        pk=[detect(m,z,g,c,cuts,radii,phases) for m,z,g,c in zip(maps,caches,grids,cubes,strict=True)]
                        match=qualifying(pk,truth,grids,paths,cubes,kernels,maps,cuts,epoch,dt,tmp)
                        if match is not None:
                            row.update(recovered=True,status='recovered',recovered_a_km=match[0],recovered_phase_deg=match[1])
                    writer.writerow(row);ff.flush();saved.append(row)
        stat=dict(seed=args.seed,population=len(elem),observable=int(use.sum()),grid_radius_km=radii.tolist(),
                  grid_phase_step_deg=.75,epoch_jd_tdb=epoch,epoch2_jd_tdb=epoch2,
                  interval_days=dt/86400,flux_levels_Jy=args.flux,matched_tolerance_px=1.5,
                  summary=search.summarize_plants(saved),input_hashes={str(f):search.file_hash(f) for f in [*paths,args.ee,Path(__file__)]})
        (args.out/'report.json').write_text(json.dumps(stat,indent=2)+'\n')


if __name__=='__main__':
    main()
