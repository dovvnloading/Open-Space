"""Compare September with the mean of May and August Level-3 source catalog fluxes at the host list."""
import argparse
import hashlib
from pathlib import Path
from urllib.request import urlopen

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u

BANDS=['f200w','f444w']
ROOT='https://stpubdata.s3.amazonaws.com/jwst/public/jw09263/L3/t/'


def run(targetfile,cache,out):
    frozen=pd.read_csv(Path(__file__).resolve().parent/'september/catalog_inputs.csv')
    if out.exists() and any(out.iterdir()): raise ValueError('Output directory must be empty')
    from acquire_nexus import make_targets
    if targetfile.name!='targets.csv': raise ValueError('Use targets.csv produced by frozen acquire_nexus inputs')
    make_targets(targetfile.parent)  # validates source hashes and byte-identical target order
    t=pd.read_csv(targetfile).drop_duplicates(['ra','dec']).reset_index(drop=True)
    target=SkyCoord(t.ra,t.dec,unit='deg')
    obslist=[2,4,6]; all_rows=[];summary=[]; inputs=[]; tests=[]
    for band in BANDS:
        flux=np.full((len(t),3),np.nan);err=flux.copy()
        for k,obs in enumerate(obslist):
            name=f'jw09263-o{obs:03}_t001_nircam_clear-{band}_cat.ecsv'
            url=f'{ROOT}o{obs:03}/{name}';p=cache/name
            if not p.exists():
                with urlopen(url,timeout=120) as response:p.write_bytes(response.read())
            receipt=frozen[(frozen.epoch==k+7)&(frozen.band==band)]
            if len(receipt)!=1 or receipt.iloc[0]['url']!=url or p.stat().st_size!=int(receipt.iloc[0]['bytes']) or hashlib.sha256(p.read_bytes()).hexdigest()!=receipt.iloc[0]['sha256']:
                raise ValueError(f'Frozen catalogue mismatch: {p}')
            cols=['label','sky_centroid.ra','sky_centroid.dec','aper_total_flux','aper_total_flux_err']
            a=pd.read_csv(p,sep=r'\s+',comment='#',usecols=cols).dropna()
            ref=SkyCoord(a['sky_centroid.ra'],a['sky_centroid.dec'],unit='deg')
            idx,dis,_=target.match_to_catalog_sky(ref)
            _,r2,_=target.match_to_catalog_sky(ref,nthneighbor=2)
            good=(dis.arcsec<=.15)&(r2.arcsec>.15)
            raw=a.aper_total_flux.to_numpy()[idx]*1e9;e=a.aper_total_flux_err.to_numpy()[idx]*1e9
            keep=good&(raw>0)&(e>0)
            flux[keep,k]=raw[keep];err[keep,k]=e[keep]
            inputs.append(dict(epoch=k+7,band=band,url=url,bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
            for j in np.flatnonzero(keep&(t['sample']=='high_z')):
                all_rows.append(dict(id=t.id[j],epoch=k+7,band=band,source_label=a.label.to_numpy()[idx[j]],offset_arcsec=dis.arcsec[j],flux_njy=raw[j],pipeline_err_njy=e[j]))
        ctrl=(t['sample']=='control').to_numpy()
        bright=ctrl&np.all(flux/err>15,axis=1)
        if bright.sum()<20:raise ValueError('Too few catalogue calibration sources')
        sf=flux[bright];scale=np.median(sf/np.median(sf,axis=1)[:,None],axis=0)
        scale/=np.median(scale);flux/=scale;err/=scale
        delta=flux[:,2]-flux[:,:2].mean(axis=1)
        sig=np.sqrt(err[:,2]**2+np.sum(err[:,:2]**2,axis=1)/4)
        ratio=delta/sig;lum=flux[:,:2].mean(axis=1)
        bins=np.nanquantile(lum[ctrl],[.25,.5,.75]);bins=np.digitize(lum,bins)
        stat=np.full(len(t),np.nan)
        for b in range(4):
            z=ratio[(bins==b)&ctrl&np.isfinite(ratio)]
            if len(z)<8: continue
            med=np.median(z);infl=max(1.,1.4826*np.median(np.abs(z-med)))
            stat[bins==b]=(ratio[bins==b]-med)/infl
        for j in np.flatnonzero(np.isfinite(stat)):
            tests.append(dict(id=t.id[j],sample=t['sample'][j],ra=t.ra[j],dec=t.dec[j],z=t.z[j],band=band,delta_njy=delta[j],scale_agnostic_score=stat[j]))
        summary.append(dict(band=band,sources=len(t),controls=bright.sum(),scale_ep7=scale[0],scale_ep8=scale[1],scale_ep9=scale[2]))
    out.mkdir(exist_ok=True,parents=True)
    pd.DataFrame(inputs).to_csv(out/'catalog_inputs.csv',index=False)
    pd.DataFrame(summary).to_csv(out/'calibration.csv',index=False,float_format='%.7g')
    pd.DataFrame(all_rows).to_csv(out/'host_catalog_lightcurves.csv',index=False,float_format='%.8g')
    df=pd.DataFrame(tests)
    df.to_csv(out/'catalog_tests.csv',index=False,float_format='%.8g')
    pivot=df.pivot(index=['id','sample'],columns='band',values='scale_agnostic_score').dropna()
    s=pivot.sum(axis=1);ctrl=s[pivot.index.get_level_values('sample')=='control']
    denom=max(np.sqrt(2.),1.4826*np.median(np.abs(ctrl-np.median(ctrl))))
    pivot['joint_score']=s/denom
    pivot['passes']=(pivot.f200w*pivot.f444w>0)&(pivot[BANDS].abs().min(axis=1)>=3)&(pivot.joint_score.abs()>=7)
    pivot.to_csv(out/'joint_catalog_tests.csv',float_format='%.8g')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('targets',type=Path);p.add_argument('cache',type=Path);p.add_argument('out',type=Path)
    args=p.parse_args();args.cache.mkdir(exist_ok=True,parents=True);run(args.targets,args.cache,args.out)
