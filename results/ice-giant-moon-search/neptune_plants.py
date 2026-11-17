"""Inject off-grid moons into the one-visit Neptune images and repeat the two-filter screen."""
import argparse
import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import numpy as np
from scipy.stats import binomtest
from geometry import generate
from neptune_scan import photometry
from search import Cuts, TemplateCache, clipped_stamp, file_hash, load_npz, prepare, trail_kernel, validate_cube


def passing(measurements):
    bands=[photometry(*x) for x in measurements]
    n1,s1,d1,_=bands[0];n2,s2,d2,_=bands[1]
    return (n1>=3)&(n2>=3)&(s1>=5)&(s2>=5)&(d1>=3)&(d2>=3)&(np.hypot(s1,s2)>=8)


def run(paths,out,nplant,levels,seed,ee=None,broad=False):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    cut=Cuts(psf_support=.8); rng=np.random.default_rng(seed)
    truth=np.zeros((nplant,6));truth[:,0]=rng.uniform(109500,110500,nplant)
    if broad:
        if nplant % 10:
            raise ValueError("broad population requires an equal number of plants in each of ten radial strata")
        truth[:,0]=np.repeat(np.arange(30000,130000,10000),nplant//10)+rng.uniform(0,10000,nplant)
    truth[:,-1]=rng.uniform(0,360,nplant)
    el=out/'injected_elements.csv'
    np.savetxt(el,truth,delimiter=',',header='a_km,e,i_deg,node_deg,peri_deg,mean_deg',comments='')
    data=[];measured=[];eligible=np.ones(nplant,bool);epoch=None
    injected_psfs=[];fractions=[]
    retained=None
    with TemporaryDirectory() as tmp:
        for source in paths:
            cube=load_npz(source)
            if str(cube['unmasked_control_body'])!='None':
                raise ValueError('blind plants require known-moon masks')
            psfs=validate_cube(cube,min_visits=1);maps=prepare(cube,psfs,cut)
            if ee:
                from jwst8975.pilot_inputs import radial_psf
                model=radial_psf(ee,float(cube["scale_arcsec"]),str(cube["band"]).split("/")[1])
                injected_psfs.append([model]*len(psfs));fractions.append(float(model.sum()))
            else:
                injected_psfs.append(psfs);fractions.append(1.)
            path=Path(tmp)/'bank.npz';generate(source,path,grid=(109500,110500,500,1),epoch=epoch)
            bank=load_npz(path);epoch=float(bank['epoch_jd_tdb'])
            generate(source,path,elements=el,model='circular',epoch=epoch);plants=load_npz(path)
            if broad:
                if retained is None:
                    from recovery_grid import compatible_rows
                    retained=compatible_rows(cube,plants,epoch)
                    np.savetxt(out/"recovery_compatible_grid.csv",retained,delimiter=",",header="a_km,e,i_deg,node_deg,peri_deg,mean_deg",comments="")
                generate(source,path,elements=out/"recovery_compatible_grid.csv",model="circular",epoch=epoch)
                bank=load_npz(path)
            cache=TemplateCache(cube,bank,psfs)
            measured.append(cache.measure(maps,len(bank['x']),cut))
            truthcache=TemplateCache(cube,plants,psfs)
            _,w=truthcache.measure(maps,len(plants['x']),cut);eligible &= (w>0).sum(1)>=3
            data.append((cube,psfs,bank,plants,cache))
    # Inject together only if they remain at least 50 pixels apart in every
    # frame of both bands, even including shutter trails and a 12-pixel glare
    # smoothing kernel. Test against the full narrow annulus blind grid.
    groups=[]
    for i in np.flatnonzero(eligible):
        for group in groups:
            separated=True
            for _,_,_,p,_ in data:
                # Relative sky separation between times varies much more slowly
                # than one trail length. Subtract each trail's full diameter.
                gap=np.hypot(p['x'][group]-p['x'][i],p['y'][group]-p['y'][i])
                widths=np.hypot(np.ptp(p['x_path'][group],axis=-1),np.ptp(p['y_path'][group],axis=-1))
                single=np.hypot(np.ptp(p['x_path'][i],axis=-1),np.ptp(p['y_path'][i],axis=-1))
                if np.min(gap-widths-single)<50:
                    separated=False;break
            if separated:
                group.append(i);break
        else:
            groups.append([i])
    base=passing(measured)
    results=[]
    for level in levels:
        row=['outside_two_filter_coverage' if not yes else 'missed' for yes in eligible]
        for group in groups:
            fitted=[]
            for band_i,(cube,psfs,bank,plants,cache) in enumerate(data):
                images=cube['image'].copy()
                for i in group:
                    for j,psf in enumerate(injected_psfs[band_i]):
                        px,py=plants['x_path'][i,j],plants['y_path'][i,j]
                        box=clipped_stamp(images.shape[1:],*trail_kernel(px,py,psf))
                        if box is not None:
                            sl,stamp=box
                            images[j][sl]+=level*1e-6*stamp
                maps=prepare({**cube,'image':images},psfs,cut)
                fitted.append(cache.measure(maps,len(bank['x']),cut))
            hits=passing(fitted)
            for i in group:
                # Compare midpoint astrometry of the fitted track to injected
                # truth in every frame of both filters, including unsupported ones.
                near=np.ones(len(hits),bool)
                for _,_,b,p,_ in data:
                    near &= np.hypot(b['x']-p['x'][i],b['y']-p['y'][i]).max(1)<=2.
                if (base&near).any():
                    row[i]='preexisting'
                elif (hits&near).any():
                    row[i]='recovered'
        results.extend((level,i,truth[i,0],truth[i,-1],row[i]) for i in range(nplant))
        print(level,row.count('recovered'),int(eligible.sum()),nplant,flush=True)
        with (out/'plants.csv').open('w',newline='') as f:
            w=csv.writer(f);w.writerow(['flat_flux_uJy','id','a_km','phase_deg','status']);w.writerows(results)
    summary=[]
    for level in levels:
        r=[s[-1] for s in results if s[0]==level];k=r.count('recovered');total=int(eligible.sum())
        ci=binomtest(k,total).proportion_ci() if total else None
        summary.append(dict(flat_uJy=level,eligible=total,total=nplant,k=k,
                 conditional_fraction=k/total if total else None,all_fraction=k/nplant,
                 conditional_binomial_95=[ci.low,ci.high] if ci else None))
    (out/'summary.json').write_text(json.dumps(dict(seed=seed,epoch_jd_tdb=epoch,
        cube_sha256=[file_hash(p) for p in paths],injection_grid_km=([30000,130000] if broad else [109500,110500]),
        blind_grid_trials=(72360 if broad else 1080),evaluated_for_recovery=len(data[0][2]["x"]),
        radius_step_km=500.,phase_step_deg=1.,match_radius_px=2.,min_shared_plant_separation_px=50.,
        geometric_groups=len(groups),flux_convention='same total flux density (Jy) in F140M and F210M',
        injection_psf=('STScI radial encircled energy' if ee else 'normalized Gaussian'),
        stamp_fractions=fractions,ee_table_sha256=file_hash(ee) if ee else None,
        recovery=summary),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('f140');p.add_argument('f210');p.add_argument('out')
    p.add_argument('--n',type=int,default=36);p.add_argument('--seed',type=int,default=2739140)
    p.add_argument('--levels-uJy',type=float,nargs='+',default=[.04,.08,.16,.32,.64])
    p.add_argument('--ee',type=Path,help='STScI short-wavelength encircled-energy table for off-model plants')
    p.add_argument('--broad',action='store_true',help='equal counts in ten radial strata across the full blind grid')
    a=p.parse_args();run([a.f140,a.f210],a.out,a.n,a.levels_uJy,a.seed,a.ee,a.broad)
