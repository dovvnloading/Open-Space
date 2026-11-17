"""Off-grid recovery on the continuous bound orbit family using Gaia-star PSFs."""
import csv
import json
import numpy as np
from astropy.time import Time
from scipy.ndimage import map_coordinates
from scipy.spatial import cKDTree

from dense import BASE,CACHE,HERE,OUT,R,design,kernel,p
from mosaic import xy
from vis_pixels import MODEL,TREF,SC

MAG=np.array([19.,20.,21.,22.,23.,23.75,24.,24.5,25.,25.25,25.5,25.75,26.,26.25,26.5,26.75,27.,27.25,27.5])
SEED=20261117


def catalog():
    geom=json.loads((BASE/'geometry.json').read_text())
    geox,geoy=geom['x0'],geom['y0']
    rows=list(csv.DictReader((HERE/'field_des.response').open()))
    ra=np.array([float(q['ra']) for q in rows]); de=np.array([float(q['dec']) for q in rows])
    x,y=xy(ra,de)
    return cKDTree(np.column_stack((x-geox,y-geoy)))


def scores(d,w,src,ks,counts,seed,mag,threshold):
    # Score the same integer sky centre searched by dense.py, at the phase offset
    # and actual kinematics of each off-grid draw. Scoring only the nearest centre
    # is conservative relative to using the best response of neighbouring centres.
    b,n,h,l=d.shape
    # Each Gaussian motion template is identical over this small sky patch.
    masks=np.ones((len(ks[0])+1,n),'f4')
    masks[1:]-=np.eye(n,dtype='f4')
    filters=[]
    for mask in masks:
        ww=w*mask[None,:,None,None]
        sw=np.maximum(ww.sum(axis=1),1e-20)
        banks=[]
        for K in ks:
            centered=K[None,:,:,:]-(ww*K[None,:,:,:]).sum(axis=1)[:,None,:,:]/sw[:,None,:,:]
            div=np.sqrt(np.maximum((ww*centered**2).sum(axis=(1,2,3)),1e-20))
            banks.append((ww*centered/div[:,None,None,None]).reshape(b,-1))
        filters.append(banks)
    coeff=np.stack(filters,axis=0).transpose(2,0,1,3).reshape(b,-1,n*h*l)
    # Direct image addition includes source shot noise. counts is in electrons
    # per VIS=25 source; the registered archive pixels retain their real noise.
    rng=np.random.default_rng(seed)
    blank=(coeff@d.reshape(b,-1,1))[...,0].reshape(b,len(masks),len(ks))
    out=[]
    for m in mag:
        scale=10**(.4*(25-m))
        if np.ndim(scale):
            scale=scale[:,None,None,None]
        raw=rng.poisson(np.maximum(0,src*scale*counts[None,:,None,None])).astype('f4')
        synthetic=d+raw/counts[None,:,None,None]
        z=(coeff@synthetic.reshape(b,-1,1))[...,0].reshape(b,len(masks),len(ks))
        # The trigger must persist when any single exposure is withheld.
        out.append((z[:,0].max(axis=1)>=threshold)&(z[:,1:].max(axis=2).min(axis=1)>=3))
    return np.stack(out,axis=1),blank[:,0].max(axis=1)


def wilson(k,n):
    p0=k/n;z=1.95996398454
    c=(p0+z*z/(2*n))/(1+z*z/n)
    w=z*np.sqrt(p0*(1-p0)/n+z*z/(4*n*n))/(1+z*z/n)
    return c-w,c+w


def main():
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("--scenario",choices=["nominal","centroid_stress","july_scans","fixed_grid"],default="nominal")
    scenario=parser.parse_args().scenario
    names,ts,bank,_=design()
    ks=np.array([[np.outer(kernel(dy),kernel(dx)) for dx,dy in mot] for mot in bank])
    g=json.loads((BASE/'geometry.json').read_text())
    D=[np.load(BASE/(name+'_data.npy'),mmap_mode='r') for name in names]
    W=[np.load(BASE/(name+'_weight.npy'),mmap_mode='r') for name in names]
    field=np.load(OUT/'eligible.npy',mmap_mode='r'); hull=np.load(OUT/'support.npy',mmap_mode='r')
    psf=np.load(CACHE/'measured_psf.npy')
    receipt=list(csv.DictReader((HERE/'mosaic_wcs.csv').open()))
    # Per-frame zeropoints and exposure times; native VIS gain 3.48 e-/ADU.
    counts=np.array([3.48*np.median([float(r['ADU_at_VIS25']) for r in receipt if r['quadrant'].startswith(name)]) for name in names])
    rng=np.random.default_rng(SEED)
    draws=rng.uniform(size=(8192,3))
    offsets=None
    inp=json.loads((HERE.parent/'orbit_input.json').read_text())
    if scenario=="centroid_stress":
        # IRAS formal ellipse, AKARI uncalibrated illustrative disk.
        pa=np.radians(79.); major,minor=rng.normal(size=(2,len(draws)))
        iras=np.column_stack((18*major*np.sin(pa)+7*minor*np.cos(pa),
                              18*major*np.cos(pa)-7*minor*np.sin(pa)))
        rr=40*np.sqrt(rng.uniform(size=len(draws))); ang=rng.uniform(0,2*np.pi,len(draws))
        offsets=[iras,np.column_stack((rr*np.cos(ang),rr*np.sin(ang)))]
    elif scenario=="july_scans":
        inp['iras_coadd_interval_utc']=['1983-07-09T22:07:00','1983-07-10T01:34:00']
    samp=p.family(inp) if scenario=="fixed_grid" else p.family(inp,parameters=draws,angular_offsets=offsets)
    ra,de=p.sky(samp,TREF.tt.jd-2451545,SC(TREF.tdb.jd-2451545))
    xs,ys=xy(ra,de);xs-=g['x0'];ys-=g['y0']
    xi,yi=np.rint(xs).astype(int),np.rint(ys).astype(int)
    safe=(xi>=R)&(xi<g['width']-R)&(yi>=R)&(yi<g['height']-R)
    ii=np.clip(xi,0,g['width']-1);jj=np.clip(yi,0,g['height']-1)
    inside=safe&(hull[jj,ii]!=0)
    admissible=inside&(field[jj,ii]!=0)
    # A fixed 1.5 arcsec avoidance radius around independently measured DES
    # sources protects against PSF/registration residuals on stationary objects.
    cat=catalog()
    sep,_=cat.query(np.column_stack((xs,ys)))
    final=admissible&(sep>15)
    drifts=[]
    for t in ts:
        utc=Time(t,format='mjd',scale='utc')
        a,b=p.sky(samp,utc.tt.jd-2451545,SC(utc.tdb.jd-2451545))
        x,y=xy(a,b)
        drifts.append((x-g['x0']-xi, y-g['y0']-yi))
    dx,dy=np.meshgrid(np.arange(-R,R+1),np.arange(-R,R+1))
    accepted=np.flatnonzero(final)
    recovered=np.zeros((len(ra),len(MAG)),dtype=bool)
    stored=np.zeros(len(ra),dtype='f4')
    # A recovered injection must also beat any pre-existing local neighbor
    # that could suppress the source in the 25-pixel peak finder.
    from scipy.ndimage import maximum_filter
    observed_score=np.load(OUT/'score.npy',mmap_mode='r')
    crowd=np.maximum(8.,maximum_filter(observed_score,size=25)+0.01)
    for offset in range(0,len(accepted),48):
        ids=accepted[offset:offset+48]
        sel=(yi[ids,None,None]+dy,xi[ids,None,None]+dx)
        data=np.array([a[sel] for a in D]).transpose(1,0,2,3)
        weight=np.array([a[sel] for a in W]).transpose(1,0,2,3)
        src=[]
        for i,(x,y) in enumerate(drifts):
            coords=[dy[None,:,:]-y[ids,None,None]+psf.shape[1]//2,
                    dx[None,:,:]-x[ids,None,None]+psf.shape[2]//2]
            src.append(map_coordinates(psf[i],coords,order=1,cval=0,prefilter=False))
        src=np.array(src).transpose(1,0,2,3)
        recovered[ids],stored[ids]=scores(data,weight,src,ks,counts,SEED+offset+1,MAG,crowd[yi[ids],xi[ids]])
    prefix='offgrid' if scenario=='nominal' else 'offgrid_'+scenario
    with (HERE/(prefix+'_recovery.csv')).open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['VIS_AB','n_bound','n_eligible','n_recovered','fraction_of_bound','wilson95_lo','wilson95_hi','fraction_of_eligible'])
        for mag,ok in zip(MAG,recovered.T):
            k=int(ok.sum());wc=wilson(k,len(ra))
            w.writerow([mag,len(ra),len(accepted),k,round(k/len(ra),6),*[round(s,6) for s in wc],round(k/len(accepted),6)])
    if scenario=="fixed_grid":
        with (HERE/'fixed_grid_trials.csv').open('w',newline='') as f:
            w=csv.writer(f)
            w.writerow(['grid_id','ra_ref_deg','dec_ref_deg','in_hull','quality_overlap','eligible_after_DES_mask','blank_score']+[f'recovered_VIS{v:g}' for v in MAG])
            for j in range(len(ra)):
                w.writerow([j,f'{ra[j]:.8f}',f'{de[j]:.8f}',int(inside[j]),int(admissible[j]),int(final[j]),round(float(stored[j]),5) if final[j] else '',*map(int,recovered[j])])
    out=dict(scenario=scenario,seed=SEED,proposed=samp[-1]['tested'],n_bound=len(ra),in_searched_hull=int(inside.sum()),
             quality_overlap=int(admissible.sum()),after_DES_veto=len(accepted),
             max_float_identity_delta=float(np.max(np.abs(stored[accepted]-np.load(OUT/'score.npy',mmap_mode='r')[yi[accepted],xi[accepted]]))),
             band='Euclid VIS AB',min_loo_score=3,threshold=8)
    print(out)
    (HERE/(prefix+'_summary.json')).write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__': main()
