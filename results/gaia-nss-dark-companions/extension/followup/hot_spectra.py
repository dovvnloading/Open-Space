"""Resolve H-alpha wing dilution and hot-template RV sensitivity for the 470-day SB1."""
import argparse
import csv
import gzip
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter1d, percentile_filter

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'vetting'))
from spectra import C, URL, read_header, table_dtype
from inject import BASE, FILES

SID='275582835635329024'
PAIRS=[('764807168','719508155'),('776307168','707107168')]
MODELS={**{'wave':FILES['wave']}, **{'hot'+str(t):f'PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte{t}-4.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits' for t in [10000,12000]},
        **{'giant'+str(t):f'PHOENIX-ACES-AGSS-COND-2011/Z{z}/lte0{t}-2.50{z}.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits' for t,z in [(5000,'-0.0'),(4800,'-0.5')]}}


def fetch(cache,name,url):
    p=cache/name
    if not p.exists():p.write_bytes(urllib.request.urlopen(url,timeout=100).read())
    return p,dict(file=name,url=url,sha256=hashlib.sha256(p.read_bytes()).hexdigest())


def raw_red(p):
    with gzip.open(p) as f:
        ph=read_header(f)
        assert ph['HELIO']=='T' and ph['VACUUM']=='T' and ph['NAXIS']=='0'
        h=read_header(f); f.read((int(h['NAXIS1'])+2879)//2880*2880)
        h=read_header(f)
        assert h['EXTNAME']=='COADD_R'
        r=np.frombuffer(f.read(int(h['NAXIS1'])),dtype=table_dtype(h))[0]
    q=(r['ORMASK']==0)&(r['IVAR']>0)&np.isfinite(r['FLUX'])
    return r['WAVELENGTH'][q].astype(float),r['FLUX'][q].astype(float)


def scaled(x,y,yp):
    # Reference points in flanking continuum bands, on the visible-star rest wavelength.
    anchors=[]
    for lo,hi in [(6526,6540),(6588,6602)]:
        use=(x>=lo)&(x<=hi)
        anchors.append(np.percentile(y[use],yp))
    m=(anchors[1]-anchors[0])/(6595-6533)
    return y/(anchors[0]+m*(x-6533))


def run(cache):
    meta=[];arr={}
    for key,suffix in MODELS.items():
        p,m=fetch(cache,key+'.fits',BASE+suffix);meta.append(m)
        with p.open('rb') as f:
            h=read_header(f)
            arr[key]=np.fromfile(f,dtype=('>f8' if key=='wave' else '>f4'),count=int(h['NAXIS1']))
    grid=np.exp(np.arange(np.log(6490),np.log(6635),2/C))
    tpl={}
    for k in MODELS:
        if k=='wave': continue
        y=np.interp(grid,arr['wave'],arr[k])
        for rot in ([5] if k.startswith('giant') else [15,150]):
            v=np.arange(-int(rot/2),int(rot/2)+1)*2/rot
            s=2*.4*np.sqrt(1-v*v)+.5*np.pi*.6*(1-v*v)
            tpl[k,rot]=gaussian_filter1d(np.convolve(y,s/s.sum(),mode='same'),C/7500/2.355/2)
    cat=[r for r in csv.DictReader((ROOT/'extension/lamost/exposures.csv').open()) if r['source_id']==SID]
    fits=[]; pix=[]
    for pair in PAIRS:
        for obsid in pair:
            p,m=fetch(cache,obsid+'.fits.gz',URL+obsid);meta.append(m)
            obsw,obs=raw_red(p)
            v1=np.median([float(r['rv_kms']) for r in cat if r['obsid']==obsid])
            v2=-51.526-1.25*(v1+51.526)
            rest=obsw/(1+v1/C)
            use=(rest>=6526)&(rest<=6602)
            rest,obs,obsw=rest[use],obs[use],obsw[use]
            for qp in [85,95]:
                y=scaled(rest,obs,qp)
                comp={k:scaled(rest,np.interp(obsw/(1+(v1 if k[0].startswith('giant') else v2)/C),grid,v),qp)
                      for k,v in tpl.items()}
                for g in ('giant4800','giant5000'):
                    cool=comp[g,5]
                    for hot in ('hot10000','hot12000'):
                        for rot in [15,150]:
                            second=comp[hot,rot]
                            dl=np.abs(rest-6564.61)
                            # Exclude the H-alpha core and obvious metal blends in either giant model.
                            wing=(dl>2.5)&(dl<16)&(comp['giant4800',5]>.94)&(comp['giant5000',5]>.94)
                            z=second[wing]-cool[wing]; diff=y[wing]-cool[wing]
                            beta=np.dot(z,diff)/np.dot(z,z)
                            slope=lambda b:float(np.sqrt(np.mean((diff-b*z)**2)))
                            fits.append(dict(obsid=obsid,primary_rv=v1,secondary_rv_assumed=v2,primary=g,hot=hot,
                                vsini=rot,continuum_percentile=qp,wing_pixels=int(wing.sum()),
                                fitted_halpha_continuum_fraction=beta,rms_f0=slope(0),rms_fit=slope(beta),
                                rms_f05=slope(.05),rms_f10=slope(.1),rms_f20=slope(.2),rms_f30=slope(.3)))
                            if qp==95 and g=='giant5000' and hot=='hot10000' and rot==15:
                                pix.extend(dict(obsid=obsid,primary_rest_A=round(w,4),observed=round(a,5),
                                           cool_model=round(c,5),hot_model=round(hotf,5),wing_pixel=int(ok))
                                           for w,a,c,hotf,ok in zip(rest,y,cool,second,wing))
    with (HERE/'halpha_pixels.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,list(pix[0]),lineterminator='\n');w.writeheader();w.writerows(pix)
    return dict(source_id=SID,inputs=meta,assumed_m1_over_m2=1.25,assumed_systemic_kms=-51.526,
                models=fits)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cache',required=True,type=Path)
    a=p.parse_args();a.cache.mkdir(parents=True,exist_ok=True)
    (HERE/'halpha_fits.json').write_text(json.dumps(run(a.cache),indent=2)+'\n')
