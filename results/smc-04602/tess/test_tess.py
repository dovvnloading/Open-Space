"""Test the archived OGLE timing in MAST TESScut sectors for OGLE-SMC-ECL-04602.

Run acquire.py 94 95 and then OPENBLAS_NUM_THREADS=1 python test_tess.py 94 95.
"""
from pathlib import Path
import hashlib, json, sys, warnings
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares
from scipy.linalg import eigh

RA,DEC=15.7324167,-72.8732778
P,T0=1.7295874,2460001.4679
C,W=np.array([.00084,.61597]),np.array([.021,.023])
def wrap(x):return (x+.5)%1-.5
def tmpl(p,centers=C,widths=W):
    return np.exp(-.5*(wrap(p[:,None]-centers)/widths)**2)
def md(a):
    return 1.4826*np.median(abs(a-np.median(a,axis=0)),axis=0)

def reduce(path):
    with fits.open(path) as h:
        d=h[1].data
        t=d['TIME'].astype('d');img=d['FLUX'].astype('d')
        good=(d['QUALITY']==0)&np.isfinite(t)&np.all(np.isfinite(img),axis=(1,2))
        t,img=t[good],img[good]
        cadence=float(np.median(np.diff(t)))
        x,y=WCS(h[2].header,naxis=2).world_to_pixel_values(RA,DEC);x,y=float(x),float(y)
        dates=[h[0].header['DATE-OBS'],h[0].header['DATE-END']]
        sector=h[0].header['SECTOR']
    p=((t+2457000-70/86400-T0)/P)%1
    yy,xx=np.indices((11,11));dist=np.hypot(xx-x,yy-y).ravel()
    Y=img.reshape(-1,121);ref=dist>3.1
    r=Y[:,ref]-np.median(Y[:,ref],axis=0);r=np.clip(r/np.maximum(md(r),.5),-8,8)
    _,v=eigh(r.T@r,subset_by_index=(r.shape[1]-6,r.shape[1]-1));basis=r@v
    block=np.floor((t-t[0])/.65).astype(int)
    clean=np.empty_like(Y)
    for j in np.unique(block):
        ix=block==j
        D=np.column_stack([np.ones(ix.sum()),t[ix]-t[ix].mean(),basis[ix]])
        clean[ix]=Y[ix]-D@np.linalg.lstsq(D,Y[ix],rcond=1e-8)[0]
    q=np.median(abs(clean[:,ref])/np.maximum(md(clean[:,ref]),.3),axis=1)
    keep=q<np.median(q)+6*max(md(q),.05);t,p,clean=t[keep],p[keep],clean[keep]
    aps=np.array([(abs(xx-round(x))<=1)&(abs(yy-round(y))<=1),(dist<=np.sort(dist)[3]).reshape(11,11),(dist==min(dist)).reshape(11,11)])
    cycles=np.floor((t+2457000-70/86400-T0)/P+.17).astype(int)
    nb=144;full=np.arange(nb);cent=(full+.5)/nb
    def bins(ix):
        ids=np.floor(p[ix]*nb).astype(int)%nb
        count=np.bincount(ids,minlength=nb)
        data=np.zeros((nb,121));np.add.at(data,ids,clean[ix])
        return data,count
    imgfold,count=bins(np.ones(len(t),bool))
    # Use summed pixel values per bin, normalized by integration counts, with an
    # independent sinusoidal baseline. Eclipses are defined as positive deficits.
    def design(ph,centers=C):
        return np.column_stack([np.ones(len(ph)),np.cos(2*np.pi*ph),np.sin(2*np.pi*ph),-tmpl(ph,centers)])
    D=design(cent)
    ww=np.sqrt(np.maximum(count,1))
    coeff=np.linalg.lstsq(D*ww[:,None],imgfold/ww[:,None],rcond=None)[0][-2:]
    pixel=[]
    for co in coeff:
        z=np.argmax(co);r,c=divmod(z,11);pixel.append({'peak':[int(c),int(r)],'peak_dist':dist[z],'nearest_amplitude':co[np.argmin(dist)],'peak_amplitude':co[z]})
    ifs=[]
    for n in np.unique(cycles):
        ix=cycles==n
        a,b=bins(ix)
        if cadence > .01:  # 30-minute FFIs cannot populate 144 bins in one orbit
            if b.sum()*cadence < .78*P or (b@tmpl(cent)).min()*cadence < .025: continue
        elif (b>0).sum()<.78*nb or (b@tmpl(cent)).min()<7: continue
        ifs.append((a@aps.reshape(3,-1).T,b))
    return {'sector':sector,'file_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'dates':dates,'n_kept':len(t),'xy':[x,y],'per_cycle_count':len(ifs),'pixel':pixel},ifs

def evalf(arrs,apsel=0,centers=C,freewidth=False):
    fold,count=arrs[0].sum(axis=0),arrs[1].sum(axis=0)
    phase=(np.arange(len(count))+.5)/len(count)
    y=fold[:,apsel]/np.maximum(count,1);w=np.sqrt(count)
    def model(b):
        cc=b[5:7]; ww=b[7:9] if freewidth else W
        return b[0]+b[1]*np.cos(2*np.pi*phase)+b[2]*np.sin(2*np.pi*phase)-tmpl(phase,cc,ww)@b[3:5]
    b=[0.,0,0,max(1.,(np.percentile(y,80)-np.percentile(y,5))),max(1.,(np.percentile(y,80)-np.percentile(y,5))),*centers];lo=[-np.inf]*5+list(centers-.055);hi=[np.inf]*5+list(centers+.055)
    if freewidth: b.extend(W);lo.extend([.010,.010]);hi.extend([.050,.050])
    return least_squares(lambda b:(model(b)-y)*w,b,bounds=(lo,hi)).x

def summary(ifs,n=399):
    arrays=[np.stack(v) for v in zip(*ifs)]
    p=(np.arange(144)+.5)/144
    amps=[]
    for a,count in ifs:
        w=np.sqrt(np.maximum(count,1))
        D=np.column_stack([np.ones(144),np.cos(2*np.pi*p),np.sin(2*np.pi*p),-tmpl(p)])
        amps.append(np.linalg.lstsq(D*np.sqrt(count[:,None]),a/w[:,None],rcond=None)[0][-2:])
    amps=np.stack(amps)
    v=evalf(arrays);vf=evalf(arrays,freewidth=True);rng=np.random.default_rng(4602)
    bs=[evalf([a[idx] for a in arrays]) for idx in rng.integers(len(ifs),size=(n,len(ifs)))]
    sec=np.array(bs)[:,6]-np.array(bs)[:,5]
    errs=np.array(bs)[:,5:7]
    return {'fixed_template_depth_3x3':amps[:,:,0].mean(axis=0).tolist(),'sem_orbits_3x3':(amps[:,:,0].std(axis=0,ddof=1)/np.sqrt(len(amps))).tolist(),'positive_orbits':(amps[:,:,0]>0).sum(axis=0).tolist(),
    'free_fit_depth':v[3:5].tolist(),'centers':v[5:7].tolist(),'centers95':np.quantile(errs,[.025,.975],axis=0).tolist(),'separation':v[6]-v[5],'separation95':np.quantile(sec,[.025,.975]).tolist(),'freewidth_centers':vf[5:7].tolist(),'free_widths':vf[7:9].tolist(),'four_pixel_centers':evalf(arrays,apsel=1)[5:7].tolist(),'single_pixel_centers':evalf(arrays,apsel=2)[5:7].tolist()}

if __name__=='__main__':
    out={};allifs=[]
    for s in map(int,sys.argv[1:]):
        meta,ifs=reduce(Path(f's{s:03}.fits'));allifs+=ifs
        np.savez_compressed(f'folds-{s}.npz',flux=np.stack([f for f,c in ifs]),count=np.stack([c for f,c in ifs]))
        out[str(s)]=meta|summary(ifs)
    out['combined']=summary(allifs,1200)
    print(json.dumps(out,indent=2))
