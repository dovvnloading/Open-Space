"""Blind native integration stacks and a targeted cycle-count sweep through a supplied late-visit position."""
import argparse
import csv
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import warnings

import numpy as np
from scipy import ndimage
base=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(base));sys.path.insert(0,str(base/'jwst9482'))
from link_visits import geometry,search,trials
import positive_control as native
from prepare import inverse_sky


def build(path, meta_cache, pixels, out):
    cube=search.load_npz(path)
    records=json.loads(Path(str(path)+'.inputs.json').read_text())
    yy,xx=np.indices(cube['image'].shape[1:])
    center=cube['center_xy'];scale=float(cube['scale_arcsec'])
    rho=np.hypot(xx-center[0],yy-center[1]);ann=(rho/10).astype(int)
    maps=[];support=[];times=[];ras=[];decs=[];deltas=[];ids=[];known=[];frame=[]
    for j,r in enumerate(records):
        m=native.metadata(meta_cache,r['exposure_id'].removesuffix('_cal'))
        residual,x0,y0=native.wide_images(pixels,m)
        ra,dec=m['wcs'].all_pix2world(r['center_native_x'],r['center_native_y'],0)
        wr,wd=inverse_sky(ra,dec,-(xx-center[0])*scale/geometry.ARCSEC,(yy-center[1])*scale/geometry.ARCSEC)
        x,y=m['wcs'].all_world2pix(wr,wd,0,quiet=True)
        coords=np.array([y-y0,x-x0])
        ts=cube['path_mjd'][j]
        for i,t in enumerate(m['times']):
            a=residual[i];ok=np.isfinite(a)
            good=ok.astype(float); raw=np.nan_to_num(a)
            def smooth(sigma):
                return ndimage.gaussian_filter(raw,sigma)/np.maximum(ndimage.gaussian_filter(good,sigma),.03)
            z=ndimage.map_coordinates(smooth(.85)-smooth(4),coords,order=1,mode='constant')
            usable=(ndimage.map_coordinates(good,coords,order=1,mode='constant')>=.9)&(rho>75)
            # Robust empirical per-annulus spatial noise at PSF scale.
            adjust=np.zeros_like(z)
            for ring in np.unique(ann[usable]):
                local=usable & (ann==ring)
                vals=z[local];m0=np.median(vals);s=1.4826*np.median(abs(vals-m0))
                if s>0 and len(vals)>=20:adjust[local]=(vals-m0)/s
                else:usable[local]=False
            maps.append(adjust.astype('f4'));support.append(usable)
            time=t+np.linspace(-.5,.5,5)*m['sec']/86400
            times.append(time);ras.append(np.interp(time,ts,cube['path_ra'][j]))
            decs.append(np.interp(time,ts,cube['path_dec'][j]));deltas.append(np.interp(time,ts,cube['path_delta'][j]))
            kp=cube['known_path'][:,j]
            known.append(np.array([[np.interp(t,ts,z[:,axis]) for axis in range(2)] for z in kp]))
            ids.append(f"{r['exposure_id']}_{i+1:02d}");frame.append(j)
    np.savez_compressed(out, image=np.array(maps),valid=np.array(support),frame=frame,
        path_mjd=times,path_ra=ras,path_dec=decs,path_delta=deltas,exposure_id=ids,
        scale_arcsec=scale,center_xy=center,planet=cube['planet'],known_id=cube['known_id'],known_xy=known)


def scores(c,bank):
    ntrial=len(bank['x']);n=len(c['image'])
    val=np.zeros((ntrial,n));ok=np.zeros((ntrial,n),bool)
    for j,(im,mask) in enumerate(zip(c['image'],c['valid'],strict=True)):
        coords=np.array([bank['y_path'][:,j],bank['x_path'][:,j]])
        v=ndimage.map_coordinates(im,coords,order=1,mode='constant')
        w=ndimage.map_coordinates(mask.astype('f4'),coords,order=1,mode='constant')>=.9
        count=w.sum(1);ok[:,j]=count>=3
        val[:,j]=np.where(ok[:,j],(v*w).sum(1)/np.maximum(count,1),0)
    def stack(sel):
        nn=ok[:,sel].sum(1)
        return np.divide(val[:,sel].sum(1),np.sqrt(nn),out=np.zeros(ntrial),where=nn>0)
    mid=n//2
    return stack(slice(None)),stack(slice(None,mid)),stack(slice(mid,None)),ok.sum(1),val,ok


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('cube',type=Path);p.add_argument('meta',type=Path);p.add_argument('pixels',type=Path)
    p.add_argument('cache',type=Path);p.add_argument('out',type=Path)
    p.add_argument('--grid',type=float,nargs=4,default=[73500,78000,100,.5])
    p.add_argument('--late',type=float,nargs=3,help='target phase deg, later epoch JD TDB, local phase half-width')
    a=p.parse_args()
    if not a.cache.exists():build(a.cube,a.meta,a.pixels,a.cache)
    c=search.load_npz(a.cache)
    lo,hi,dr,dp=a.grid;radii=np.arange(lo,hi+dr*1e-9,dr)
    if a.late:
        phi,epoch,half=a.late;phases=np.arange(phi-half,phi+half+.00001,dp)
    else:phases=np.arange(0,360,dp);epoch=None
    fields=np.zeros((4,len(radii),len(phases)))
    a.out.mkdir(exist_ok=False,parents=True)
    with TemporaryDirectory() as td:
        tmp=Path(td)
        for start in range(0,len(radii),max(1,1000//len(phases))):
            group=radii[start:start+max(1,1000//len(phases))]
            rr,pp=np.meshgrid(group,phases,indexing='ij')
            elem=np.zeros((rr.size,6));elem[:,0]=rr.ravel();elem[:,-1]=pp.ravel()
            b=trials(tmp/'e.csv',a.cache,tmp/'g.npz',elem,epoch)
            if epoch is None:epoch=float(b['epoch_jd_tdb'])
            z,e,l,n,*_=scores(c,b)
            for k,v in enumerate([z,e,l,n]):fields[k,start:start+len(group)]=v.reshape(rr.shape)
        z,e,l,n=fields;sel=(n>=len(c['image'])*.5)&(e>=2)&(l>=2)
        vals=np.where(sel,z,-1e9)
        maxima=ndimage.maximum_filter(vals,(3,5),mode=('constant','constant' if a.late else 'wrap'),cval=-1e9)
        positions=np.argwhere((vals>=4)&(vals==maxima))
        positions=sorted(positions,key=lambda ij:-vals[tuple(ij)])
        with open(a.out/'local_peaks.csv','w',newline='') as f:
            w=csv.writer(f);w.writerow(['a_km','phase_deg','score','early','late','n','closest_known_body','median_sep_px','min_sep_px'])
            for i,j in positions:
                elem=np.array([[radii[i],0,0,0,0,phases[j]]])
                b=trials(tmp/'e.csv',a.cache,tmp/'g.npz',elem,epoch)
                xy=np.stack([b['x'][0],b['y'][0]],axis=-1)
                sep=np.linalg.norm(xy[:,None,:]-c['known_xy'],axis=2)
                md=np.median(sep,axis=0);nearest=np.argmin(md)
                w.writerow([radii[i],phases[j],z[i,j],e[i,j],l[i,j],int(n[i,j]),
                            c['known_id'][nearest],md[nearest],np.min(sep[:,nearest])])
    stat=dict(epoch_jd_tdb=epoch,radii=a.grid,late_fixed_phase=a.late,trials=z.size,valid_half=int((n>=len(c['image'])*.5).sum()),
              local_peaks=len(positions),native_inputs_sha256=search.file_hash(a.cube),maps_sha256=search.file_hash(a.cache),
              integrations=len(c['image']),score_definition='Sum of per-integration local-MAD responses divided by sqrt(N), not Gaussian significance.')
    (a.out/'parameters.json').write_text(json.dumps(stat,indent=2)+'\n')


if __name__=='__main__':main()
