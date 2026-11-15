"""Independent Kepler pointing-polynomial test. Uses no other pixels as regressors.
Uses only NumPy, SciPy and Astropy. Supply the same cache as replay.py.
"""
import argparse
import json
from pathlib import Path
import tarfile

import numpy as np
from astropy.io import fits
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

def image_data(cache):
    with fits.open(cache/'k2.fits.gz') as h:
        table=h[1].data
        time=table['TIME'].copy()+4833
        cad=table['CADENCENO'].copy()
        good=table['QUALITY']==0
        images=np.array(table['FLUX'],dtype=float)
        origin=np.array([h[1].header['1CRV5P'],h[1].header['2CRV5P']])
    ra,dec=271.27245833,-27.82268
    grid=np.loadtxt(cache/'g92_52')
    basis=np.array([1,ra,dec,ra*ra,ra*dec,dec*dec])
    positions=np.full((len(cad),2),np.nan)
    i=np.searchsorted(cad,grid[:,0])
    keep=i<len(cad)
    keep[keep]&=cad[i[keep]]==grid[keep,0]
    positions[i[keep]]=np.column_stack([grid[keep,4:10]@basis,grid[keep,10:16]@basis])-origin
    good&=np.isfinite(positions).all(axis=1)
    return time[good],positions[good],images[good],origin

def ground_checks(cache,t,flux):
    summary={}
    with tarfile.open(cache/'kmtobs') as archive:
        for member in archive:
            if not member.name.endswith('_I.pysis'):
                continue
            data=np.loadtxt(archive.extractfile(member))
            data=data[(data[:,5]<8)&(data[:,2]>0)&(data[:,2]<2000)]
            nightly=[]
            for day in range(7531,7546):
                ground=data[np.floor(data[:,0])==day]
                if len(ground)==0:
                    continue
                nearest=np.abs(ground[:,0,None]-t).argmin(axis=1)
                paired=abs(t[nearest]-ground[:,0])<0.04
                if np.any(paired):
                    nightly.append([np.median(-ground[paired,1]),np.median(flux[nearest[paired]])])
            if len(nightly)>=4:
                earth,space=np.array(nightly).T
                summary[member.name]={'days':len(nightly),'rho':round(float(spearmanr(earth,space).statistic),3)}
    return summary

def audit(cache):
    t,pos,imgs,origin=image_data(cache)
    yy,xx=np.indices(imgs.shape[1:])
    valid=np.isfinite(imgs).all(axis=0)
    x,y=xx[valid],yy[valid]
    images=imgs[:,valid]
    train=((t<7532)|(t>7544))&(t<7568)
    heldout=t>=7568
    burst=(t>7533)&(t<7541)
    scatter=lambda a:1.4826*np.median(np.abs(a-np.median(a)))
    dxy=(pos-np.mean(pos[train],axis=0))/np.std(pos[train],axis=0)
    dx,dy=dxy.T
    def apertures(residual):
        def flux(u,v):
            use=((x-42-u)**2+(y-1.5-v)**2)<5**2
            d2=(pos[:,0,None]+u-x[use])**2+(pos[:,1,None]+v-y[use])**2
            return (np.exp(-d2/(2*.85**2))*residual[:,use]).sum(axis=1)
        return flux
    def stats(f):
        contrast=np.median(f[burst])-np.median(f[heldout])
        # Local temporal correlation makes this a ranking statistic, not a Gaussian significance.
        return float(contrast),float(contrast/scatter(f[heldout]))
    result={}
    for degree in (1,2,3):
        terms=[dx**j*dy**(i-j) for i in range(degree+1) for j in range(i+1)]
        design=np.column_stack(terms)
        coefficients=np.linalg.lstsq(design[train],images[train],rcond=None)[0]
        residual=images-design@coefficients
        aperture=apertures(residual)
        f=aperture(0,0)
        contrast,score=stats(f)
        local=(x-42)**2+(y-1.5)**2<8**2
        j=np.argmax(np.where(local,np.median(residual[burst],axis=0),-np.inf))
        controls=[]
        for x0 in range(4,46,4):
            for y0 in range(4,53,4):
                if np.hypot(x0-42,y0-1.5)>=6:
                    controls.append(stats(aperture(x0-42,y0-1.5)))
        controls=np.array(controls)
        result[str(degree)]=dict(event_median=round(float(np.median(f[burst])),2),
            heldout_median=round(float(np.median(f[heldout])),2),heldout_rms=round(float(np.std(f[heldout])),2),
            contrast=round(contrast,2),contrast_over_holdout_MAD=round(score,3),n_controls=len(controls),
            n_controls_contrast_at_least_target=int(np.sum(controls[:,0]>=contrast)),
            n_controls_normalized_contrast_at_least_target=int(np.sum(controls[:,1]>=score)),
            brightest_pixel_in_local_radius=[int(x[j]+origin[0]),int(y[j]+origin[1])],
            KMT_site_field=ground_checks(cache,t,f))
    # Unsuccessful linear/nearest-pointing reductions are retained as model sensitivity controls.
    inner=(x-42)**2+(y-1.5)**2<8**2
    ref=(t>7545)&(t<7568)
    distances,index=cKDTree(pos[ref]).query(pos,k=8)
    diff=images[:,inner]-np.median(images[ref][:,inner][index],axis=1)
    d2=(pos[:,0,None]-x[inner])**2+(pos[:,1,None]-y[inner])**2
    f=(diff*np.exp(-d2/(2*.85**2))).sum(axis=1)
    result['nearest_8_points']=dict(event_median=round(float(np.median(f[burst])),2),
        heldout_median=round(float(np.median(f[heldout])),2),heldout_rms=round(float(np.std(f[heldout])),2),
        median_distance_to_8th_neighbor_pixels=round(float(np.median(distances[burst,-1])),4))
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache',type=Path,required=True)
    parser.add_argument('--out',type=Path,default=Path('detector-audit.json'))
    args=parser.parse_args()
    args.out.write_text(json.dumps(audit(args.cache),indent=2)+'\n')
