"""Blind moving-PSF search in the continuous 2024-12-27 sky footprint."""
import csv
import json
from concurrent.futures import ThreadPoolExecutor
from itertools import product

from astropy.time import Time
import numpy as np
from scipy.ndimage import correlate1d, maximum_filter
from scipy.spatial import ConvexHull
from matplotlib.path import Path as Polygon

from mosaic import xy, world
from vis_pixels import CACHE,HERE,MODEL,TREF,SC,RA0,DEC0,PIXSCALE,p

BASE=CACHE/'mosaic'
OUT=CACHE/'dense'
OUT.mkdir(exist_ok=True)
R=8
SIGMA=1.0
BORDER=16


def design():
    dates={r['product']:r['utc_mid'] for r in csv.DictReader((HERE/'coverage.csv').open())}
    byexp={i.name.split('_')[0]:str(np.load(i)['product']) for i in CACHE.glob('*_stamps.npz')}
    names=sorted(byexp)
    ts=[Time(dates[byexp[n]],scale='utc') for n in names]
    tracks=[]
    for t in ts:
        a,b=p.sky(MODEL,t.tt.jd-2451545,SC(t.tdb.jd-2451545))
        tracks.append(np.column_stack(((a-RA0)*np.cos(np.radians(DEC0))*3600/PIXSCALE,
                                       (b-DEC0)*3600/PIXSCALE)))
    tr=np.array(tracks)
    speed=np.linalg.norm(tr[-1],axis=1)
    # Three full tracks spanning the speed range. Coordinate drift over
    # 2.12 hours is nearly linear. Test the actual vector mismatch below.
    chosen=[int(np.abs(speed-v).argmin()) for v in np.linspace(speed.min(),speed.max(),3)]
    bank=tr[:,chosen]
    mismatch=np.linalg.norm(tr[:,:,None,:]-bank[:,None,:,:],axis=-1).max(axis=0).min(axis=-1)
    return names,np.array([t.mjd for t in ts]),bank.transpose(1,0,2),float(mismatch.max())


def kernel(d):
    u=np.arange(-R,R+1,dtype='f4')
    k=np.exp(-0.5*((u-d)/SIGMA)**2)
    return k/k.sum()


def corr(a,kx,ky):
    return correlate1d(correlate1d(a,kx,axis=1,mode='constant'),ky,axis=0,mode='constant')


def score_tile(d,w,bank):
    """Exact profiled weighted least squares at all integer grid centers."""
    sw=np.maximum(w.sum(axis=0),1e-20)
    mean=(w*d).sum(axis=0)/sw
    residual=w*(d-mean)
    scores=[]
    for motion in bank:
        k=[(kernel(dx),kernel(dy)) for dx,dy in motion]
        num=np.zeros_like(sw);den=np.zeros_like(sw)
        for i,(kx,ky) in enumerate(k):
            num+=corr(residual[i],kx,ky)
            den+=corr(w[i]*(1.-w[i]/sw),kx*kx,ky*ky)
            for j in range(i):
                q=w[i]*w[j]/sw
                if q.any():
                    den-=2*corr(q,kx*k[j][0],ky*k[j][1])
        scores.append(num/np.sqrt(np.maximum(den,1e-20)))
    return np.array(scores)


def main():
    geom=json.loads((BASE/'geometry.json').read_text())
    x0,y0=geom['x0'],geom['y0']
    names,ts,bank,mismatch=design()
    height,width=geom['height'],geom['width']
    D=[np.load(BASE/(s+'_data.npy'),mmap_mode='r') for s in names]
    W=[np.load(BASE/(s+'_weight.npy'),mmap_mode='r') for s in names]
    xs,ys=xy(RA0,DEC0)
    coords=np.column_stack((xs,ys))
    hull=Polygon(coords[ConvexHull(coords).vertices])
    shape=(height,width)
    Z=np.lib.format.open_memmap(OUT/'score.npy',mode='w+',dtype='f4',shape=shape)
    valid=np.lib.format.open_memmap(OUT/'eligible.npy',mode='w+',dtype='u1',shape=shape)
    origin=np.lib.format.open_memmap(OUT/'support.npy',mode='w+',dtype='u1',shape=shape)
    win=np.lib.format.open_memmap(OUT/'bank.npy',mode='w+',dtype='u1',shape=shape)
    tiles=[(a,b) for a in range(0,height,384) for b in range(0,width,384)]
    def tile(job):
        j,i=job
        je,ie=min(j+384,height), min(i+384,width)
        jj,ii=np.meshgrid(np.arange(j,je),np.arange(i,ie),indexing='ij')
        mask=hull.contains_points(np.column_stack((ii.ravel()+x0,jj.ravel()+y0))).reshape(ii.shape)
        if not mask.any(): return
        ja,ia=max(0,j-BORDER),max(0,i-BORDER)
        jb,ib=min(height,je+BORDER),min(width,ie+BORDER)
        sr=np.s_[ja:jb,ia:ib];dest=np.s_[j:je,i:ie];core=np.s_[j-ja:je-ja,i-ia:ie-ia]
        w=np.array([f[sr] for f in W]);d=np.array([f[sr] for f in D])
        # Require adequate support near the predicted positions. A strict zero
        # flag on each input interpolation neighbor was enforced in mosiac.py.
        from scipy.ndimage import uniform_filter
        local=np.array([uniform_filter((wi>0).astype('f4'),size=7,mode='constant')>.7 for wi in w])
        possible=(local.sum(axis=0)>=3)&((np.max(np.where(local,ts[:,None,None],-np.inf),axis=0)-
                  np.min(np.where(local,ts[:,None,None],np.inf),axis=0))*1440>=40)
        good=possible[core]&mask
        origin[dest]=mask
        if not good.any(): return
        z=score_tile(d,w,bank)[:,core[0],core[1]]
        mx=z.max(axis=0)
        Z[dest]=np.where(good,mx,0)
        win[dest]=z.argmax(axis=0)
        valid[dest]=good
    with ThreadPoolExecutor(max_workers=3) as ex:
        list(ex.map(tile,tiles))
    for dest in (Z,valid,origin,win): dest.flush()
    test=(valid!=0)
    # Discrete maxima 2.5 arcsec apart. Never interpret their number alone as
    # a statistical significance. Static gradients and cosmic rays have tails.
    peaks=(Z==maximum_filter(Z,size=25))&test&(Z>=8.)
    yy,xx=np.nonzero(peaks)
    a,b=world(xx+x0,yy+y0)
    with (HERE/'dense_triggers.csv').open('w',newline='') as f:
        cw=csv.writer(f);cw.writerow(['x','y','ra_deg','dec_deg','score','template'])
        for k in np.argsort(Z[yy,xx])[::-1]:
            cw.writerow([xx[k],yy[k],a[k],b[k],Z[yy[k],xx[k]],win[yy[k],xx[k]]])
    out=dict(hull_cells=int(origin.sum()),eligible_cells=int(test.sum()),scale_arcsec=PIXSCALE,
             triggers_above_eight=len(xx),max_score=float(Z[test].max()),min_max_score=float(Z[test].min()),
             score_percentiles=np.percentile(Z[test],[.1,1,50,99,99.9]).tolist(),
             templates=bank.tolist(),template_mismatch_max_pixel=mismatch)
    (HERE/'dense_summary.json').write_text(json.dumps(out,indent=2)+'\n')
    print(out)


if __name__=='__main__': main()
