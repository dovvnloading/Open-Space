"""Compare moving-source and named-moon positions on the same individual integrations."""
import argparse
import csv
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import numpy as np
from scipy import ndimage, optimize
import warnings

sys.path.insert(0,str(Path(__file__).resolve().parent))
from native_orbits import geometry,search,trials


def stack(c,xy,radius=12):
    dy,dx=np.mgrid[-radius:radius+1,-radius:radius+1]
    images=[]
    for m,valid,(x,y) in zip(c['image'],c['valid'],xy,strict=True):
        coords=[dy+y,dx+x]
        s=ndimage.map_coordinates(m,coords,order=1,mode='constant')
        ok=ndimage.map_coordinates(valid.astype('f4'),coords,order=1,mode='constant')>.85
        images.append(np.where(ok,s,np.nan))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',RuntimeWarning)
        st=[np.nanmedian(chunk,axis=0) for chunk in [images,images[:len(images)//2],images[len(images)//2:]]]
    q=dx**2+dy**2; ring=(q>=6**2)&(q<=12**2); small=q<=3**2
    def normalize(img):
        off=np.nanmedian(img[ring]);noise=1.4826*np.nanmedian(abs(img[ring]-off))
        return (img-off)/noise
    norm=list(map(normalize,st));im=norm[0]
    valid=small&np.isfinite(im)
    if not valid.any():return dict(n=0,center=None,peak=None,dx=None,dy=None,early=None,late=None)
    p=np.unravel_index(np.nanargmax(np.where(valid,im,np.nan)),im.shape)
    return dict(n=int(np.isfinite(np.asarray(images)[:,radius,radius]).sum()),center=float(im[radius,radius]),
                peak=float(im[p]),dx=int(p[1]-radius),dy=int(p[0]-radius),early=float(norm[1][p]),late=float(norm[2][p]))


def fit_xy(c,xy):
    y0,x0=np.mgrid[-3:4,-3:4]
    out=[]
    for m,valid,(x,y) in zip(c['image'],c['valid'],xy,strict=True):
        coords=[y+y0,x+x0]
        vals=ndimage.map_coordinates(m,coords,order=1,mode='constant')
        v=ndimage.map_coordinates(valid.astype('f4'),coords,order=1,mode='constant')>.85
        valid=v&(x0*x0+y0*y0<=9)
        if not valid.any():
            out.append([np.nan,np.nan,np.nan]);continue
        iy,ix=np.unravel_index(np.argmax(np.where(valid,vals,-np.inf)),vals.shape)
        peak=vals[iy,ix]
        if peak<2:
            out.append([np.nan,np.nan,peak]);continue
        # Positive PSF-core centroid of a local thresholded peak.
        close=(x0-x0[iy,ix])**2+(y0-y0[iy,ix])**2<=2.25
        weight=np.where(v&close,np.maximum(vals-max(0,peak*.35),0),0)
        if weight.sum()<=0:
            out.append([x+x0[iy,ix],y+y0[iy,ix],peak]);continue
        out.append([x+(weight*x0).sum()/weight.sum(),y+(weight*y0).sum()/weight.sum(),peak])
    return np.asarray(out)


def intrinsic(c,xy):
    cx,cy=c['center_xy'];scale=c['scale_arcsec']/geometry.ARCSEC
    m=c['path_mjd'].shape[1]//2
    los,east,north=geometry.tangent_basis(c['path_ra'][:,m],c['path_dec'][:,m])
    h=geometry.pole('Uranus',float(np.median(c['path_mjd']))+2400000.5)
    u,v,_=geometry.basis('Uranus',float(np.median(c['path_mjd']))+2400000.5)
    e=-(xy[:,0]-cx)*scale;n=(xy[:,1]-cy)*scale
    transverse=e[:,None]*east+n[:,None]*north
    ray=los+transverse
    delta=c['path_delta'][:,m]*geometry.AU
    pos=delta[:,None]*((los@h/(ray@h))[:,None]*ray-los)
    return np.hypot(pos@u,pos@v),np.rad2deg(np.arctan2(pos@v,pos@u))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('map1',type=Path);p.add_argument('map2',type=Path)
    p.add_argument('late_astrometry',type=Path);p.add_argument('out',type=Path)
    p.add_argument('--a-km',type=float,default=75511.75)
    p.add_argument('--phase',type=float,default=254.55348)
    p.add_argument('--epoch',type=float,default=2460955.765061979)
    a=p.parse_args();a.out.mkdir(exist_ok=False,parents=True)
    late_rows=list(csv.DictReader(open(a.late_astrometry)))
    report={}
    for j,path in enumerate([a.map1,a.map2]):
        c=search.load_npz(path)
        with TemporaryDirectory() as td:
            t=Path(td);elem=np.array([[a.a_km,0,0,0,0,a.phase]])
            b=trials(t/'e.csv',path,t/'b.npz',elem,a.epoch)
        seed=np.stack([b['x'][0],b['y'][0]],axis=-1)
        fit=fit_xy(c,seed)
        if j==1:
            xy=np.full((len(c['image']),2),np.nan)
            for r in late_rows:
                mjd=float(r['mjd_utc'])
                i=np.argmin(abs(c['path_mjd'][:,2]-mjd))
                if abs(c['path_mjd'][i,2]-mjd)>1e-6:raise ValueError('wrong late astrometry')
                xy[i]=[float(r['x_sky']),float(r['y_sky'])]

        else:xy=fit[:,:2]
        measured=np.isfinite(xy).all(axis=1)
        R,phi=intrinsic(c,xy)
        bodies={}
        for body in [708,712,713,714,715,725,727,75052]:
            k=np.flatnonzero(c['known_id']==body).item()
            pos=c['known_xy'][:,k]
            rad,ang=intrinsic(c,pos)
            sep=np.linalg.norm(pos-xy,axis=1)
            diff=(phi-ang+180)%360-180
            bodies[str(body)]=dict(forced_position=stack(c,pos),median_separation_px=float(np.nanmedian(sep)),
                median_radius_difference_km=float(np.nanmedian(R-rad)),median_phase_difference_deg=float(np.nanmedian(diff)),
                phase_rms_deg=float(np.nanstd(diff)),radius_difference_mad_km=float(np.nanmedian(abs((R-rad)-np.nanmedian(R-rad)))))
        report[str(j)]=dict(orbit_seed_stack=stack(c,seed),n_measured=int(measured.sum()),
            median_deprojected_radius_km=float(np.nanmedian(R)),radius_MAD_km=float(np.nanmedian(abs(R-np.nanmedian(R)))),
            named_objects=bodies)
        if j==0:
            with open(a.out/'february_astrometry.csv','w',newline='') as f:
                w=csv.writer(f);w.writerow(['exposure_integration','mjd_utc','x_seed','y_seed','x_measured','y_measured','peak_local_MAD','deprojected_radius_km'])
                for k,id0 in enumerate(c['exposure_id']):
                    w.writerow([id0,c['path_mjd'][k,2],*seed[k],*fit[k],R[k]])
    def clean(v):
        if isinstance(v,dict):return {k:clean(x) for k,x in v.items()}
        if isinstance(v,float) and not np.isfinite(v):return None
        return v
    (a.out/'comparison.json').write_text(json.dumps(clean(report),indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
