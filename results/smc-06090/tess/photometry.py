"""Forced TESS photometry at the prior OGLE ephemeris and local image checks."""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import least_squares
from astropy.io import fits
from astropy.wcs import WCS
from datetime import datetime,timedelta
from fetch import HERE,RA,DEC,url

PER=1.1696719
T0=2460000.0369
LATE=(-.0061,.5337)
EARLY=(.035,.496)
W=.025

def wrap(x):return (x+.5)%1-.5

def template(ph,centers):
    return np.stack([np.exp(-.5*(wrap(ph-c)/W)**2) for c in centers],axis=-1)

def md(a):return np.median(np.abs(a-np.median(a,axis=0)),axis=0)*1.4826

def lightcurves(file,return_native=False,inject=None):
    with fits.open(file) as h:
        tab=h[1].data
        t=tab['TIME'].astype(float); flags=tab['QUALITY'];flux=tab['FLUX'].astype(float)
        wcs=WCS(h[2].header,naxis=2)
        x,y=map(float,wcs.world_to_pixel_values(RA,DEC))
        s=h[0].header['SECTOR']
        source_hash=hashlib.sha256()
        for arr,dt in ((t,'<f8'),(flags,'<i8'),(flux,'<f4')):source_hash.update(np.asarray(arr,dtype=dt).tobytes())
        good=(flags==0)&np.isfinite(t)&np.all(np.isfinite(flux),axis=(1,2))
        skycoord=[x,y];dt=float(np.nanmedian(np.diff(t)))
        meta=dict(url=url(s),sector=s,science_sha256=source_hash.hexdigest(),
                  camera=h[0].header['CAMERA'],ccd=h[0].header['CCD'],initial=len(t),
                  quality_zero_complete=int(good.sum()),cadence_days=dt,
                  target_pixel=skycoord,range_btjd=[float(t[good][0]),float(t[good][-1])])
        if inject is not None:
            signal,psf,amplitude=inject
            flux-=amplitude*signal[:,None,None]*psf[None,:,:]
        t=t[good];img=flux[good]
    yy,xx=np.indices(img.shape[1:]); dist=np.hypot(xx-x,yy-y)
    Y=img.reshape(len(t),-1)
    ref=(dist>3.1).ravel()
    refdata=Y[:,ref]
    # Estimate six detector/background modes outside the science aperture.
    r=refdata-np.median(refdata,axis=0)
    r/=np.maximum(md(r),.5)
    r=np.clip(r,-8,8)
    _,v=eigh(r.T@r,subset_by_index=(r.shape[1]-6,r.shape[1]-1))
    basis=r@v
    # Continuously varying sky and spatial systematics plus segment-specific
    # linear trends (blocks 0.65 d, much longer than one eclipse core).
    seg=np.floor((t-t[0])/.65).astype(int)
    out=np.empty_like(Y)
    for k in np.unique(seg):
        ix=seg==k;u=t[ix]-np.mean(t[ix])
        D=np.column_stack([np.ones(ix.sum()),u,basis[ix]])
        out[ix]=Y[ix]-D@np.linalg.lstsq(D,Y[ix],rcond=1e-8)[0]
    # Symmetric clipping of extreme cadences using outer pixels only.
    q=np.median(np.abs(out[:,ref])/np.maximum(md(out[:,ref]),.3),axis=1)
    keep=q<np.median(q)+6*max(md(q),.05)
    t=t[keep];out=out[keep].reshape(-1,*dist.shape)
    half=(np.abs(xx-round(x))<=1)&(np.abs(yy-round(y))<=1)
    four=dist<=np.sort(dist.ravel())[3]
    pix=dist==dist.min()
    aps=np.array([out[:,ap].sum(axis=1) for ap in (half,four,pix)]).T
    # Make coarse 30-minute samples for independent re-fitting and phase control.
    # Average templates on the native timestamps before comparing primary depths.
    delta=70./86400   # BJD(TDB) to approximate HJD(UTC); residual CCD timing <=~minute
    ph=wrap((t+2457000.-delta-T0)/PER)
    block=np.floor((t-t[0])/(1/48)).astype(int)
    # Keep all native measurements for phase-fold stacks; use per-orbit fits for uncertainty.
    cycle=np.floor((t+2457000.-delta-T0)/PER+.2).astype(int)
    if return_native:
        return t,ph,aps,cycle
    amps=[];indices=[]
    for c in np.unique(cycle):
        sel=cycle==c
        if sel.sum()*dt<.78*PER:continue
        G=template(ph[sel],LATE)
        if np.sum(G,axis=0).min()*dt < .028:continue
        D=np.column_stack([np.ones(sel.sum()),wrap(ph[sel]-.25),-G])
        co=np.linalg.lstsq(D,aps[sel],rcond=None)[0][-2:]
        amps.append(co);indices.append(int(c))
    # The phase zero and half-cycle remained labeled through both OGLE releases.
    # Every phase gets an exposure-weighted mean and count; this compact image cube
    # can re-run the phased spatial and timing tests offline.
    nb=120
    binid=np.floor((ph%1)*nb).astype(int)
    count=np.bincount(binid,minlength=nb)
    stacked=np.stack([np.mean(out[binid==i],axis=0) if count[i] else np.full(dist.shape,np.nan) for i in range(nb)])
    meta.update(kept=len(t),cycles_used=len(amps),range_utc_approx=[(datetime(2000,1,1,12)+timedelta(days=float(j+2457000-2451545))).date().isoformat() for j in (t[0],t[-1])])
    return meta,stacked,count,np.array([half,four,pix]),np.array(amps),np.array(indices)

def evaluate(fileout):
    sector_files=sorted((HERE/'raw').glob('s*.fits'))
    records={};compact={}
    for f in sector_files:
        meta,images,counts,ap,amp,cycles=lightcurves(f);s=str(meta['sector'])
        for k,v in [('images',images),('counts',counts),('masks',ap),('orbital_depths',amp),('cycles',cycles)]:compact['s'+s+'_'+k]=v
        records[s]=meta
    np.savez_compressed(HERE/'folded_pixels.npz',**compact)
    (HERE/'tess_inputs.json').write_text(json.dumps(records,indent=2)+'\n')
    return records,compact

def measure(records,compact,save=True):
    b=(np.arange(120)+.5)/120
    rng=np.random.default_rng(6090)
    rec={}
    def fit_stack(ph,y,centers,weights):
        D=np.column_stack([np.ones_like(ph),-template(ph,centers)])
        mask=np.all(np.isfinite(y),axis=1)&(weights>0)
        v=np.sqrt(weights[mask,None])
        return np.linalg.lstsq(v*D[mask],y[mask]*v,rcond=None)[0]
    for s,r in records.items():
        imgs=compact['s'+s+'_images'];cnt=compact['s'+s+'_counts'];masks=compact['s'+s+'_masks'];amps=compact['s'+s+'_orbital_depths']
        if not len(amps):continue
        A=np.stack([imgs[:,m].sum(axis=1) for m in masks],axis=1)
        coef=fit_stack(b, A, LATE,cnt)
        early=fit_stack(b,A,EARLY,cnt)
        off=[]
        for shift in np.arange(1,120)/120:
            c=np.asarray(LATE)+shift
            if min(abs(wrap(c[:,None]-np.array(LATE))).ravel())<.1:continue
            off.append(fit_stack(b,A,c,cnt)[1:].mean(axis=0))
        spatial=fit_stack(b,imgs.reshape(120,-1),LATE,cnt)[1:].mean(axis=0).reshape(11,11)
        py,px=np.unravel_index(np.argmax(spatial),spatial.shape)
        avg=amps.mean(axis=0)
        rec[s]=dict(n=len(amps),orbital_mean=avg.tolist(),orbital_se=(np.std(amps,axis=0,ddof=1)/np.sqrt(len(amps))).tolist(),
                    phase_mean=coef[1:].tolist(),early=early[1:].tolist(),control_grid=np.array(off).tolist(),
                    strongest_spatial_pixel=[int(px),int(py)],strongest_spatial_depth=float(spatial[py,px]))
    groups={}
    for label,sectors in [('TESS_post_OGLE',[94,95,96]),('TESS_2023',[67,68]),('TESS_2018_20',[1,2,13,27,28])]:
        present=[s for s in sectors if str(s) in rec]
        if not present:continue
        arr=[compact[f's{s}_orbital_depths'] for s in present]
        combined=np.concatenate(arr)
        # Phase controls use the same fixed displacement for all constituent sectors.
        ctl=np.average([rec[str(s)]['control_grid'] for s in present],weights=[len(a) for a in arr],axis=0)
        pooled=np.average([np.mean(rec[str(s)]['phase_mean'],axis=0) for s in present],weights=[len(a) for a in arr],axis=0)
        ci=np.percentile(np.mean(np.concatenate([a[rng.integers(0,len(a),(1500,len(a)))] for a in arr],axis=1),axis=1),[2.5,97.5],axis=0)
        groups[label]=dict(sectors=present,cycles=len(combined),orbit_mean=combined.mean(axis=0).tolist(),
          ci95_orbit_resampling=ci.tolist(),phase_pooled_average_eclipses=pooled.tolist(),offphase_equal_deeper=(ctl>=pooled).sum(axis=0).tolist(),offphases=len(ctl))
    out=dict(apertures=['3x3','four nearest','nearest pixel'],centers_late=LATE,centers_early=EARLY,period=PER,
             depth_units='electrons per second; positive for a dip',sectors=rec,groups=groups)
    if save:(HERE/'tess_results.json').write_text(json.dumps(out,indent=2)+'\n')
    for s,r in rec.items():
        print(s,'cycles',r['n'],'depth 3x3',np.round(np.array(r['orbital_mean'])[:,0],2),'SE',np.round(np.array(r['orbital_se'])[:,0],2),'pixel max',r['strongest_spatial_pixel'])
    for label,g in groups.items():print(label,'combined',np.round(np.array(g['orbit_mean'])[:,:2],2),'ci95',np.round(np.array(g['ci95_orbit_resampling'])[:,:,0],2),'off deeper',g['offphase_equal_deeper'],g['offphases'])
    return out

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--from-folded',action='store_true');a=parser.parse_args()
    if a.from_folded:
        with np.load(HERE/'folded_pixels.npz') as store:
            measure(json.loads((HERE/'tess_inputs.json').read_text()),store)
    else:
        records,compact=evaluate(HERE);measure(records,compact)
