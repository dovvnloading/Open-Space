"""Search predefined spectroscopic galaxy positions for temporal F444W excursions."""
import argparse
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from scipy.ndimage import shift
from acquire_nexus import DATES, R

THRESHOLD = 7.
UNIT = 84.6159499407523  # nJy per mosaic pixel at 1 MJy/sr
YY, XX = np.mgrid[-9:10,-9:10]
CY,CX = np.mgrid[-R:R+1,-R:R+1]
OFFSETS=[(0,0)]+[(round(r*np.cos(a)),round(r*np.sin(a))) for r,n,phase in [(19,12,0.),(24,16,.5)] for a in (np.arange(n)+phase)*2*np.pi/n]


def weights(x,y):
    radius=np.hypot(XX-x,YY-y)
    f=np.zeros_like(radius)
    for dx in (np.arange(5)+.5)/5-.5:
        for dy in (np.arange(5)+.5)/5-.5:
            f+=np.hypot(XX+dx-x,YY+dy-y)<=3
    return f/25,(radius>=6)&(radius<=8.3)


def aperture(im,err,w,b):
    valid=np.isfinite(im)&np.isfinite(err)&(err>0)
    if np.any((w>0)&(~valid)) or np.sum(valid&b)<.75*b.sum():
        return np.nan,np.nan
    sky=np.median(im[valid&b])
    return UNIT*np.sum((im[w>0]-sky)*w[w>0]),UNIT*np.sqrt(np.sum((err[w>0]*w[w>0])**2))


def phot_one_ep(t,ep,cache):
    science=np.load(cache/f'ep{ep:02d}_data.npz'); err=np.load(cache/f'ep{ep:02d}_error.npz')['stamps']
    imgs=science['stamps']
    fx=science['px']-np.rint(science['px']);fy=science['py']-np.rint(science['py'])
    # Known zero-redshift bright objects in the control group anchor the PSF and registration.
    stars=np.flatnonzero((t.z==0)&(t.f444>.45)&(t['sample']=='control'))
    offsets=[];profiles=[]
    for j in stars:
        img=imgs[j]
        if not np.isfinite(img[R-11:R+12,R-11:R+12]).all():continue
        rad=np.hypot(CX-fx[j],CY-fy[j])
        sky=np.median(img[(rad>8)&(rad<11)])
        clean=img-sky
        w=np.where(rad<3.5,np.maximum(clean,0.),0.)
        if not w.sum():continue
        xc,yc=(CX*w).sum()/w.sum(),(CY*w).sum()/w.sum()
        if np.hypot(xc-fx[j],yc-fy[j])>1.4:continue
        aligned=shift(np.nan_to_num(clean),(-yc,-xc),order=3)
        norm=aligned[np.hypot(CX,CY)<=9].sum()
        if norm<=0:continue
        offsets.append((xc-fx[j],yc-fy[j]))
        profiles.append(aligned/norm)
    if len(profiles)<4:raise ValueError(f'Only {len(profiles)} PSF stars at epoch {ep}')
    dx,dy=np.median(offsets,axis=0)
    psf=np.median(profiles,axis=0)
    psf=np.where(np.hypot(CX,CY)<=9,psf,0.);psf/=psf.sum()
    n=len(t); flux=np.full((n,len(OFFSETS)),np.nan);error=np.full(n,np.nan)
    subpix=[]
    for j in range(n):
        x,y=fx[j]+dx,fy[j]+dy
        # Split centroid into an integer center and a phase for the photometric footprint.
        ix,iy=int(round(x)),int(round(y));x-=ix;y-=iy
        w,b=weights(x,y);subpix.append((ix,iy,x,y))
        for k,(ox,oy) in enumerate(OFFSETS):
            xx,yy=R+ix+ox,R+iy+oy
            if min(xx,yy)<9 or max(xx,yy)>2*R-9:continue
            sl=np.s_[yy-9:yy+10,xx-9:xx+10]
            flux[j,k],e=aperture(imgs[j][sl],err[j][sl],w,b)
            if k==0:error[j]=e
    return flux,error,psf,dict(epoch=ep,stars=len(profiles),dx=dx,dy=dy),subpix


def get_contrasts(flux,err):
    # Compare a contiguous run with the mean of all other covered epochs for that same host.
    # At least two reference epochs remain. Every local sky aperture uses the identical epochs.
    rows=[]
    for j in range(len(flux)):
        good=np.isfinite(flux[j,:,0])&np.isfinite(err[j])&(err[j]>0)
        for start in range(8):
            for n in (1,2,3):
                stop=start+n
                if stop>8:continue
                on=np.arange(start,stop);off=np.flatnonzero(good&~np.isin(np.arange(8),on))
                if not good[on].all() or len(off)<2:continue
                delta=flux[j,on].mean(axis=0)-flux[j,off].mean(axis=0)
                blank=delta[1:];blank=blank[np.isfinite(blank)]
                if len(blank)<14:continue
                sky=np.median(blank)
                sig=max(1.4826*np.median(np.abs(blank-sky)),
                    np.sqrt(np.sum(err[j,on]**2)/n**2+np.sum(err[j,off]**2)/len(off)**2))
                if sig>0:
                    rows.append((j,start,n,len(off),int(np.sum(2**off)),len(blank),delta[0]-sky,sig))
    return pd.DataFrame(rows,columns=['j','start','n','ref_epochs','reference_mask','blank_apertures','delta_njy','local_sigma_njy'])


def calibrate(t,flux,ct):
    b=np.abs(np.nanmedian(flux[:,:,0],axis=1))
    bins=np.nanquantile(b[t['sample']=='control'],[.25,.5,.75])
    ct['brightness_bin']=np.digitize(b[ct.j],bins)
    ct['raw_score']=ct.delta_njy/ct.local_sigma_njy
    ct['control_scale']=1.
    ctrl=np.asarray(t['sample']=='control')[ct.j]
    for (_, _),g in ct.groupby(['brightness_bin','n']):
        c=g.loc[ctrl[g.index],'raw_score'].dropna()
        scale=max(1.,float(1.4826*np.median(np.abs(c-np.median(c)))))
        ct.loc[g.index,'control_scale']=scale
    ct['sigma_njy']=ct.local_sigma_njy*ct.control_scale
    ct['score']=ct.delta_njy/ct.sigma_njy


def injections(t,ct,psfs,phases,cache,scales):
    hosts=t.index[t['sample']=='high_z']
    # In isolated positive one-epoch injections, measure the changed stamp via the real aperture
    # including the sky re-estimate. Local comparison apertures sit beyond the injected PSF.
    c=ct[(ct.n==1)&ct.j.isin(hosts)].reset_index(drop=True)
    hdus={};all_rows=[]
    for ep in range(8):
        hdus[ep]=(np.load(cache/f'ep{ep+1:02d}_data.npz')['stamps'],
                  np.load(cache/f'ep{ep+1:02d}_error.npz')['stamps'])
    for i in range(6):
        mag=25+i*.5
        intended=10**((31.4-mag)/2.5)
        for _,row in c.iterrows():
            j,k=int(row.j),int(row.start)
            ix,iy,x,y=phases[k][j]
            sl=np.s_[R+iy-9:R+iy+10,R+ix-9:R+ix+10]
            im,er=hdus[k][0][j][sl],hdus[k][1][j][sl]
            w,b=weights(x,y)
            p=shift(psfs[k],(y,x),order=3)[R-9:R+10,R-9:R+10]
            old=aperture(im,er,w,b)[0]
            flux_native=intended*scales[k]
            new=aperture(im+flux_native*p/UNIT,er,w,b)[0]
            score=(row.delta_njy+(new-old)/scales[k])/row.sigma_njy
            all_rows.append((t.id[j],k+1,mag,row.score,score,score>=THRESHOLD))
    return pd.DataFrame(all_rows,columns=['id','epoch','injected_ab','original_score','injected_score','recovered'])


def main(cache,out,comparisons):
    t=pd.read_csv(cache/'targets.csv');out.mkdir(exist_ok=True,parents=True)
    # Three MSAIDs have exact-coordinate aliases in these releases. Count each position once.
    dup=t.duplicated(['ra','dec'],keep='first')
    t.loc[dup,'sample']='duplicate'
    all_flux=[];all_err=[];psfs=[];meta=[];phases=[]
    for ep in range(1,9):
        f,e,p,m,c=phot_one_ep(t,ep,cache)
        all_flux.append(f);all_err.append(e);psfs.append(p);meta.append(m);phases.append(c)
    flux=np.stack(all_flux,axis=1);err=np.stack(all_err,axis=1)
    smask=(t.z==0)&(t.f444>.45)&(t['sample']=='control')
    sf=flux[smask,:,0]; ratios=sf/np.nanmedian(sf,axis=1)[:,None]
    scales=np.nanmedian(ratios,axis=0);scales/=np.nanmedian(scales)
    flux/=scales[None,:,None];err/=scales[None,:]
    for i,m in enumerate(meta):m['scale']=scales[i]
    pd.DataFrame(meta).to_csv(out/'epoch_calibration.csv',index=False)
    coords=SkyCoord(t.ra,t.dec,unit='deg')
    comp=pd.read_csv(comparisons); known=SkyCoord(comp.ra_deg,comp.dec_deg,unit='deg')
    near,sep,_=coords.match_to_catalog_sky(known)
    t['nearest_comparison']=comp.designation.to_numpy()[near];t['separation_arcsec']=sep.arcsec
    t['nearest_catalog']=comp.catalog.to_numpy()[near]
    ct=get_contrasts(flux,err)
    ct=ct[~np.asarray(dup)[ct.j]].reset_index(drop=True)
    calibrate(t,flux,ct)
    ct['id']=t.id.to_numpy()[ct.j];ct['sample']=t['sample'].to_numpy()[ct.j]
    ct.to_csv(out/'contrasts.csv.gz',index=False,float_format='%.7g',compression={'method':'gzip','mtime':0})
    aliases=t.loc[dup,['id','ra','dec','z']].merge(t.loc[~dup,['id','ra','dec']],on=['ra','dec'],suffixes=('_alias','_retained'))
    aliases.to_csv(out/'duplicate_catalog_ids.csv',index=False)
    hosts=t[t['sample']=='high_z'].copy()
    light=[]
    for j in hosts.index:
        for k in range(8):
            light.append([t.id[j],t.ra[j],t.dec[j],t.z[j],k+1,DATES[k],flux[j,k,0],err[j,k]])
    pd.DataFrame(light,columns=['id','ra','dec','z_spec','epoch','date','aperture_njy','pipeline_sigma_njy']).to_csv(out/'host_lightcurves.csv',index=False,float_format='%.9g')
    # Save one row for every targeted host including hosts with no usable contrast.
    best=ct[ct['sample']=='high_z'].sort_values('score',key=np.abs,ascending=False).drop_duplicates('id')
    hosts=hosts.merge(best[['id','start','n','score','delta_njy','sigma_njy']],on='id',how='left')
    hosts['start_epoch']=hosts.pop('start')+1
    hosts.to_csv(out/'host_screen.csv',index=False,float_format='%.9g')
    inj=injections(t,ct,psfs,phases,cache,scales)
    inj.to_csv(out/'injections.csv.gz',index=False,float_format='%.7g',compression={'method':'gzip','mtime':0})
    inj=inj[np.abs(inj.original_score)<THRESHOLD]
    inj.groupby('injected_ab').recovered.agg(['count','sum','mean']).to_csv(out/'recovery.csv')
    inp=[json.loads(f.read_text()) for f in sorted(cache.glob('ep*.json'))]
    (out/'mosaic_inputs.json').write_text(json.dumps(inp,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cache',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--comparisons',type=Path,required=True)
    args=p.parse_args();warnings.filterwarnings('ignore',category=RuntimeWarning)
    main(args.cache,args.out,args.comparisons)
