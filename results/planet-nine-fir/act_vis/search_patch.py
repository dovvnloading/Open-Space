"""Blind moving-PSF patch search and off-grid fake recovery on the ACT best-fit rate."""
import csv,json
import warnings
import numpy as np
from scipy.ndimage import map_coordinates,uniform_filter,maximum_filter
from scipy.spatial import cKDTree
from common import CACHE,HERE,HALF,N,T,MOTION,coord,world,COS,RA,DE
from dense import kernel,score_tile,R
from offgrid import scores,wilson

MAG=np.array([20.,21.,22.,23.,23.5,24.,24.5,25.,25.5,26.,26.5,27.])
SEED=264206

def injected_psf(profile, dx, dy, xcenter, ycenter):
    # Mean motion during each 560.52 s exposure is constant on this short arc.
    drift=(MOTION[-1]-MOTION[0])*(560.52/86400)/(T.mjd[-1]-T.mjd[0])
    stamps=[]
    for f in (-.375,-.125,.125,.375):
        stamps.append(map_coordinates(profile,[dy-ycenter-f*drift[1]+11,dx-xcenter-f*drift[0]+11],
                                     order=1,prefilter=False))
    return np.mean(stamps,axis=0)


def sources():
    pts=[]
    for name in ('ps1','ls_dr10','des_dr2','gaia'):
        for r in csv.DictReader((HERE/(name+'.response')).open()):
            if name=='ps1':
                a,b=float(r['raMean']),float(r['decMean'])
            else:a,b=float(r['ra']),float(r['dec'])
            if name=='gaia':
                years=(T.jyear.mean()-2016.0)
                a+=(float(r['pmra'] or 0)*years/3.6e6/np.cos(np.radians(b)))
                b+=float(r['pmdec'] or 0)*years/3.6e6
            x,y=coord(a,b)
            pts.append([x+HALF,y+HALF])
    return cKDTree(pts)


def psf_field(D,W,coverage):
    offs=np.arange(-11,12);dx,dy=np.meshgrid(offs,offs);rad=np.hypot(dx,dy)
    apcorr=np.median([float(q['apercorr']) for exp in coverage for q in exp['quadrants'] if q['apercorr'] is not None])
    samples=[[] for _ in D]
    stars=list(csv.DictReader((HERE/'gaia.response').open()))
    stars=[s for s in stars if 14.5 < float(s['phot_g_mean_mag']) < 21]
    for s in csv.DictReader((HERE/'ls_dr10.response').open()):
        if s['type'].strip()!='PSF' or not 1 < float(s['flux_r']) < 100:continue
        if any(np.hypot((float(s['ra'])-float(q['ra']))*COS,float(s['dec'])-float(q['dec']))<1/3600 for q in stars):continue
        stars.append(dict(ra=s['ra'],dec=s['dec'],pmra=0,pmdec=0))
    for s in stars:
        years=(T.jyear.mean()-2016.)
        d=float(s['dec']);a=float(s['ra'])+float(s['pmra'] or 0)*years/3.6e6/np.cos(np.radians(d))
        d+=float(s['pmdec'] or 0)*years/3.6e6
        gx,gy=coord(a,d)
        for i,(img,wt) in enumerate(zip(D,W)):
            x,y=gx+HALF,gy+HALF
            if not (14<x<N-15 and 14<y<N-15):continue
            for _ in range(3):
                z=map_coordinates(img,[y+dy,x+dx],order=1,prefilter=False)
                good=map_coordinates((wt>0).astype('f4'),[y+dy,x+dx],order=1,prefilter=False)>.999
                ring=good&(rad>8)
                if ring.sum()<30 or not good[rad<1.5].all():break
                bg=np.median(z[ring]);z-=bg
                core=z*(rad<=3)*good
                if core.sum()<=0:break
                ex,ey=np.sum(core*dx)/core.sum(),np.sum(core*dy)/core.sum()
                if np.hypot(ex,ey)>2:break
                x+=ex;y+=ey
            else:
                f=np.sum(z[good&(rad<=6.5)])
                # Sufficient S/N for the stellar profile and sensible concentration.
                if f<2 or not .25 < z[rad<2.1].sum()/f <1.12:continue
                samples[i].append(np.where(good,z/f,np.nan))
    all_=[p for row in samples for p in row]
    if len(all_)<2:raise ValueError(f'Insufficient local star PSFs: {len(all_)}')
    out=[]
    for row in samples:
        arr=row if len(row)>=3 else all_
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',RuntimeWarning)
            p=np.nan_to_num(np.nanmedian(arr,axis=0))
        p=np.maximum(0,p-np.median(p[rad>9]));p[rad>9]=0
        p*=10**(.4*apcorr)/p[rad<=6.5].sum()
        out.append(p)
    return np.array(out),[len(r) for r in samples],float(apcorr)


def main():
    coverage=json.loads((HERE/'coverage.json').read_bytes())
    D=np.array([np.load(CACHE/f'data_{i}.npy') for i in range(4)])
    W=np.array([np.load(CACHE/f'weight_{i}.npy') for i in range(4)])
    y,x=np.indices((N,N));patch=(abs(x-HALF)<=600)&(abs(y-HALF)<=600)
    local=np.array([uniform_filter((w>0).astype('f4'),size=7)>.7 for w in W])
    eligible=(local.sum(axis=0)>=3)&patch
    cat=sources()
    se=cat.query(np.array([x.ravel(),y.ravel()]).T)[0].reshape(x.shape)
    valid=eligible&(se>15)
    bank=MOTION[None]
    z=score_tile(D,W,bank)[0]
    z[~eligible]=0
    peak=(z==maximum_filter(z,size=25))&eligible&(z>=8.)
    triggers=[]
    for iy,ix in zip(*np.nonzero(peak)):
        ds=D[:,iy-10:iy+11,ix-10:ix+11]; ws=W[:,iy-10:iy+11,ix-10:ix+11]
        loo=[]
        for k in range(4):
            w=ws.copy();w[k]=0
            loo.append(float(score_tile(ds,w,bank)[0,10,10]))
        a,b=world(ix-HALF,iy-HALF)
        triggers.append(dict(ra=float(a),dec=float(b),x=int(ix),y=int(iy),z=float(z[iy,ix]),nearest_external_arcsec=float(se[iy,ix]*.1),min_loo=min(loo),loo=loo,survives=bool(se[iy,ix]>15 and min(loo)>=3)))
    (HERE/'triggers.json').write_text(json.dumps(triggers,indent=2)+'\n')
    print('cells',patch.sum(),eligible.sum(),valid.sum(),'triggers',len(triggers),'survivors',[(r['ra'],r['dec'],r['z']) for r in triggers if r['survives']],flush=True)
    psf,ns,ap=psf_field(D,W,coverage)
    np.save(CACHE/'field_psf.npy',psf)
    counts=np.array([np.median([q['exp_e_per_vis25'] for q in r['quadrants']]) for r in coverage])
    offsets=np.arange(-R,R+1);dx,dy=np.meshgrid(offsets,offsets)
    ks=np.array([[np.outer(kernel(oy),kernel(ox)) for ox,oy in MOTION]])
    # Include the exact central track as its own last draw.
    draws=np.vstack([np.random.default_rng(SEED).uniform(-600,600,(2000,2)),[0,0]])+HALF
    ix,iy=np.rint(draws.T).astype(int);validdraw=valid[iy,ix]
    truth=np.zeros((len(draws),len(MAG)),bool)
    ids=np.flatnonzero(validdraw)
    threshold=np.maximum(8.,maximum_filter(z,size=25)+.01)
    for j in range(0,len(ids),48):
        k=ids[j:j+48]
        local=[iy[k,None,None]+dy,ix[k,None,None]+dx]
        ds=D[:,local[0],local[1]].transpose(1,0,2,3);ws=W[:,local[0],local[1]].transpose(1,0,2,3)
        src=[]
        for i,(ox,oy) in enumerate(MOTION):
            px=draws[k,0]-ix[k]+ox;py=draws[k,1]-iy[k]+oy
            src.append(injected_psf(psf[i],dx[None],dy[None],px[:,None,None],py[:,None,None]))
        src=np.array(src).transpose(1,0,2,3)
        truth[k],_=scores(ds,ws,src,ks,counts,SEED+j+1,MAG,threshold[iy[k],ix[k]])
    with (HERE/'injections.csv').open('w',newline='') as f:
        o=csv.writer(f);o.writerow(['VIS_AB','random_draws','eligible_draws','recovered_draws','recovered_fraction','wilson95_lo','wilson95_hi','central_track_recovered'])
        for mag,ok in zip(MAG,truth.T):
            n=len(ok)-1;k=int(ok[:-1].sum());o.writerow([mag,n,validdraw[:-1].sum(),k,k/n,*[min(1,max(0,v)) for v in wilson(k,n)],int(ok[-1])])
    stats=dict(ra_reference_deg=float(RA.mean()),dec_reference_deg=float(DE.mean()),time_order_utc=list(T.isot),
               motion_arcsec_east_north=(MOTION*.1).tolist(),full_drift_arcsec=float(np.linalg.norm(MOTION[-1]-MOTION[0])*.1),
               patch_cells=int(patch.sum()),eligible_cells=int(eligible.sum()),after_catalogue_cells=int(valid.sum()),raw_peaks=len(triggers),
               surviving_peaks=sum(q['survives'] for q in triggers),nominal_valid_frames=local_count(W),nominal_score=float(z[HALF,HALF]) if eligible[HALF,HALF] else None,
               nominal_nearest_static_arcsec=float(se[HALF,HALF]*.1),nominal_eligible=bool(valid[HALF,HALF]),
               field_psf_stars_per_exposure=ns,aperture_correction=ap,counts_electrons_per_VIS25=counts.tolist(),seed=SEED)
    (HERE/'summary.json').write_text(json.dumps(stats,indent=2)+'\n')
    np.save(CACHE/'z.npy',z)
    print(stats,flush=True)

def local_count(w):
    return int(sum(uniform_filter((a>0).astype('f4'),size=7)[HALF,HALF]>.7 for a in w))

if __name__=='__main__':main()
