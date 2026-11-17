"""Two-epoch, two-sided matched-filter search on same-grid calibrated NIRCam mosaics."""
import argparse
import csv
import gzip
import json
from contextlib import ExitStack
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy import ndimage as nd, signal
import sep

TILE, HALO, R = 512, 40, 10
MAGS = np.arange(24, 30.01, .5)
THRESH = 7.0
Y, X = np.mgrid[-R:R+1, -R:R+1]
RAD = np.hypot(X, Y)
INNER = RAD <= 4
RING = RAD > 8


def table(path, rows, columns=None):
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=columns or list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def tiles(paths):
    """Read the gzip streams sequentially, emitting disjoint cores with overlap."""
    with ExitStack() as stack:
        streams = [stack.enter_context((gzip.open if p.suffix == '.gz' else open)(p, 'rb')) for p in paths]
        headers = [fits.Header.fromfile(f) for f in streams]
        h = headers[0]
        nx, ny = h['NAXIS1'], h['NAXIS2']
        corners = [[0, 0], [nx-1, ny-1], [0, ny-1], [nx-1, 0]]
        for other in headers:
            assert (other['NAXIS1'], other['NAXIS2'], other['BITPIX']) == (nx, ny, -32)
            assert other.get('BSCALE',1) == 1 and other.get('BZERO',0) == 0
            assert np.max(abs(WCS(other).all_world2pix(WCS(h).all_pix2world(corners,0),0)-corners)) < .001
        buffers = [np.empty((0,nx), dtype='>f4') for f in streams]
        stop = 0
        for y in range(0, ny, TILE):
            y0, y1 = max(0,y-HALO), min(ny,y+TILE+HALO)
            for i,f in enumerate(streams):
                raw = np.frombuffer(f.read((y1-stop)*nx*4), '>f4').reshape(-1,nx)
                buffers[i] = np.concatenate([buffers[i][-(stop-y0):] if stop>y0 else buffers[i][:0],raw], dtype='>f4')
            stop = y1
            for x in range(0,nx,TILE):
                x0, x1 = max(0,x-HALO), min(nx,x+TILE+HALO)
                # Most of the bounding rectangle is empty sky off the detectors.
                es = [b[:,x0:x1] for b in buffers[1::2]]
                if not np.any((es[0]>0)&(es[1]>0)):
                    continue
                core = (slice(y-y0,min(y+TILE,ny)-y0),slice(x-x0,min(x+TILE,nx)-x0))
                yield h, x0, y0, core, [np.ascontiguousarray(b[:,x0:x1],dtype='f4') for b in buffers]


def centroid_stamp(im,y,x):
    a = im[y-R:y+R+1,x-R:x+R+1].copy()
    if a.shape != RAD.shape or not np.all(np.isfinite(a)):
        return
    a -= np.median(a[RING])
    weight = np.maximum(a,0)*INNER
    if weight.sum() <= 0:
        return
    dx, dy = [(weight*v).sum()/weight.sum() for v in (X,Y)]
    if max(abs(dx),abs(dy))>1.8:
        return
    a = nd.shift(a,(-dy,-dx),order=3,mode='constant')
    total = a.sum()
    if total<=0:
        return
    return a/total, dx, dy, float(total)


def calibrate(paths,out):
    sources, stacks = [], []
    for h,x0,y0,core,(a,ea,b,eb) in tiles(paths):
        ok = np.isfinite(a)&np.isfinite(b)&(ea>0)&(eb>0)&np.isfinite(ea)&np.isfinite(eb)
        if ok.sum()<100: continue
        noise = np.median(ea[ok])
        sky = np.median(a[ok])
        objects = sep.extract(np.where(ok,a-sky,0).astype('f4'),8*noise,minarea=6)
        for obj in objects:
            x,y = int(round(obj['x'])),int(round(obj['y']))
            if not (core[0].start<=y<core[0].stop and core[1].start<=x<core[1].stop): continue
            if min(x,y,a.shape[1]-1-x,a.shape[0]-1-y)<R: continue
            if not ok[y-R:y+R+1,x-R:x+R+1].all(): continue
            if obj['b']<.7 or obj['a']/obj['b']>1.45 or obj['a']>3.8: continue
            t1,t2 = centroid_stamp(a,y,x),centroid_stamp(b,y,x)
            if t1 is None or t2 is None: continue
            s1,dx1,dy1,f1=t1
            s2,dx2,dy2,f2=t2
            snr=min(s1.max()*f1/ea[y,x],s2.max()*f2/eb[y,x])
            if snr<35 or not .82<f2/f1<1.22: continue
            compact=min(s1[R-1:R+2,R-1:R+2].sum(),s2[R-1:R+2,R-1:R+2].sum())
            if not 0.08<compact<.9 or np.hypot(dx1-dx2,dy1-dy2)>.85: continue
            sources.append(dict(x=x0+x,y=y0+y,dx=dx1-dx2,dy=dy1-dy2,snr=snr,ratio=f2/f1,compact=compact))
            stacks.append((s1,s2))
    if len(sources)<8: raise RuntimeError('insufficient common compact source calibrators')
    comp = np.array([r['compact'] for r in sources])
    # Select the sharp end of the common-source locus. Reject isolated sharp artifacts
    # through agreement of two independent epochs and flux consistency above.
    use = np.flatnonzero(comp>=np.percentile(comp,80))
    if len(use)<8: use=np.argsort(comp)[-8:]
    model = np.median(np.asarray(stacks)[use],axis=0)
    model[:,RAD>9.5]=0
    model = np.maximum(model,0)
    model /= model.sum(axis=(1,2),keepdims=True)
    for i,r in enumerate(sources): r['used_psf']=int(i in set(use))
    dx,dy=[float(np.median([sources[i][k] for i in use])) for k in ['dx','dy']]
    np.savez(out/'model.npz',psf=model,shift=np.array([dy,dx]))
    table(out/'calibrators.csv',sources)
    print('calibrators, PSF contributors, dx, dy',len(sources),len(use),dx,dy,flush=True)


def conv(a,k):
    return signal.oaconvolve(a,k,mode='same').astype('f4')


def local_levels(score, use):
    # Piecewise tile estimates from standardized matched-filter blank sky noise.
    vals=score[use]
    if vals.size<80: return np.nan,np.nan
    center=np.median(vals)
    return center, 1.4826022*np.median(abs(vals-center))


def morph(d,template,y,x,sign):
    r=6
    a=d[y-r:y+r+1,x-r:x+r+1]*sign
    k=template[template.shape[0]//2-r:template.shape[0]//2+r+1,template.shape[1]//2-r:template.shape[1]//2+r+1]
    a=a-a.mean()
    k=k-k.mean()
    den=np.sqrt((a*a).sum()*(k*k).sum())
    return float((a*k).sum()/den) if den>0 else 0.


def search(paths,out):
    model=np.load(out/'model.npz')
    p1,p2=model['psf'].astype('f4')
    dy,dx=model['shift']
    # b is shifted into a's frame. Include the interpolation response for positive
    # injections; both PSF models already include epoch drizzling and pixel sampling.
    template=signal.fftconvolve(p1,p2)
    template/=template.sum()
    filt=template[::-1,::-1]/(template**2).sum()
    ka=signal.fftconvolve(p2,filt)
    kb=signal.fftconvolve(p1,filt)
    rows,coverage,inj=[],[],[]
    rng=np.random.default_rng(7472026)
    for h,x0,y0,core,(a,ea,b,eb) in tiles(paths):
        valid=np.isfinite(a)&np.isfinite(b)&np.isfinite(ea)&np.isfinite(eb)&(ea>0)&(eb>0)
        if valid[core].sum()<80: continue
        # Exclude the entire finite PSF+matched-filter support around gaps and edges.
        valid_search=nd.minimum_filter(valid.astype('u1'),size=65,mode='constant')>0
        own=np.zeros_like(valid)
        own[core]=True
        good=valid_search&own
        if good.sum()<80: continue
        a=np.where(valid,a,0); b=np.where(valid,b,0)
        a-=np.median(a[valid]); b-=np.median(b[valid])
        b=nd.shift(b,(dy,dx),order=1,prefilter=False,mode='constant')
        d=conv(b,p1)-conv(a,p2)
        flux=conv(d,filt)
        v=conv(np.where(valid,ea**2,0),ka**2)+conv(np.where(valid,eb**2,0),kb**2)
        valid_search &= v>0
        good &= valid_search
        score=flux/np.sqrt(np.maximum(v,1e-24))
        med,mad=local_levels(score,good)
        if not np.isfinite(mad) or mad<=0: continue
        z=(score-med)/mad
        noise=np.sqrt(np.maximum(v,1e-24))*mad
        unit=h['PIXAR_SR']*1.e15  # MJy/sr * pixel solid angle to nJy
        wcs=WCS(h)
        ys,xs=np.nonzero(good)
        coverage.append(dict(x0=x0+core[1].start,y0=y0+core[0].start,n_pix=len(xs),
            noise_scale=mad,median_fluxerr_njy=np.median(noise[good])*unit))
        for sign_ in (1,-1):
            candidate=good&(sign_*z>=THRESH)&(sign_*z==nd.maximum_filter(sign_*z,size=7))
            for y,x in zip(*np.nonzero(candidate)):
                ra,dec=wcs.all_pix2world(x+x0,y+y0,0)
                corr=morph(d,template,y,x,sign_)
                opp=float(max(0,-(sign_*z[y-5:y+6,x-5:x+6]).min()))
                keep=corr>=.75 and opp<.55*abs(z[y,x])
                rows.append(dict(x=x+x0,y=y+y0,ra_deg=float(ra),dec_deg=float(dec),snr=z[y,x],
                    dflux_njy=(flux[y,x]-med*np.sqrt(v[y,x]))*unit,err_njy=noise[y,x]*unit,
                    shape_corr=corr,opposite_snr=opp,screen_pass=int(keep)))
        # Random area-weighted trials, with independent positions for each sign and
        # magnitude. Hold the measured registration, PSFs and noise map fixed.
        nt=rng.binomial(len(xs),1/15000)
        for mag in MAGS:
            for s in (1,-1):
                nrec=0
                for i in rng.integers(len(xs),size=nt):
                    y,x=ys[i],xs[i]
                    fy,fx=rng.uniform(-.5,.5,2)
                    # Linear response of image-level injection to cross convolution.
                    # b injections traverse the same bilinear alignment as b pixels.
                    native=nd.shift(p2 if s==1 else p1,(fy,fx),order=3)
                    if s==1:
                        native=nd.shift(native,(dy,dx),order=1,prefilter=False)
                    response=signal.convolve(native,p1 if s==1 else p2,mode='full',method='fft')
                    delta=s*10**((31.4-mag)/2.5)/unit*response
                    ds=d[y-20:y+21,x-20:x+21]+delta
                    ff=signal.fftconvolve(delta,filt,mode='same')
                    zz=z[y-20:y+21,x-20:x+21]+ff/noise[y-20:y+21,x-20:x+21]
                    sy,sx=np.unravel_index(np.argmax(s*zz[18:23,18:23]),(5,5))
                    sy+=18; sx+=18
                    peak=s*zz[sy,sx]
                    baseline=z[y+sy-20,x+sx-20]*s
                    if (peak>=THRESH and baseline<THRESH and (sy-20-fy)**2+(sx-20-fx)**2<=2.25**2
                        and peak>=np.max(s*zz[sy-3:sy+4,sx-3:sx+4]) and valid_search[y+sy-20,x+sx-20]
                        and morph(ds,template,sy,sx,s)>=.75 and -(s*zz[sy-5:sy+6,sx-5:sx+6]).min()<.55*peak): nrec+=1
                inj.append(dict(x0=x0,y0=y0,mag_ab=mag,sign=s,trials=nt,recovered=nrec))
    rows.sort(key=lambda r:-abs(r['snr']))
    table(out/'detections.csv',rows,['x','y','ra_deg','dec_deg','snr','dflux_njy','err_njy','shape_corr','opposite_snr','screen_pass'])
    table(out/'coverage.csv',coverage)
    table(out/'injections.csv',inj)
    print('eligible pixels, peaks, morphology passes',sum(r['n_pix'] for r in coverage),len(rows),sum(r['screen_pass'] for r in rows),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('out',type=Path)
    parser.add_argument('paths',nargs=4,type=Path,metavar='FITS',help='reference SCI, ERR, new SCI, ERR')
    parser.add_argument('--stage',choices=['calibrate','search','both'],default='both')
    args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    if args.stage!='search': calibrate(args.paths,args.out)
    if args.stage!='calibrate': search(args.paths,args.out)

if __name__=='__main__': main()
