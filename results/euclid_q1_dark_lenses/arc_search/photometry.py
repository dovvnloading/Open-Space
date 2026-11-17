"""Fetch native Q1 VIS/NISP science, flags and RMS; compare presumed image colours."""
import argparse
import csv
import itertools
import json
import math
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
from xml.etree import ElementTree as ET

import numpy as np
from PIL import Image,ImageDraw
from scipy import ndimage as ndi, optimize, stats
from search import HERE,digest,sigma,world
sys.path.insert(0,str(HERE.parent/'imaging'))

BANDS=('VIS','Y','J','H')
# Approximate effective FWHM in native 0.1 arcsec pixels. All four are smoothed to
# 0.50 arcsec only for colour fitting. Raw core and Einstein apertures use native PSFs.
WIDTH=dict(VIS=1.7,Y=3.,J=3.2,H=3.6)


def save_csv(path,rows):
    if not rows:
        return
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def read(url):
    for attempt in range(4):
        try:
            with urlopen(url,timeout=60) as f:
                return f.read()
        except OSError:
            if attempt==3:raise
            time.sleep(2**attempt)


def pix(h,ra,dec):
    a,d,ra0,d0=np.deg2rad([ra,dec,float(h['CRVAL1']),float(h['CRVAL2'])])
    delta=a-ra0
    den=np.sin(d)*np.sin(d0)+np.cos(d)*np.cos(d0)*np.cos(delta)
    offset=np.rad2deg([np.cos(d)*np.sin(delta)/den,(np.sin(d)*np.cos(d0)-np.cos(d)*np.sin(d0)*np.cos(delta))/den])
    cd=np.array([[float(h[f'CD{i}_{j}']) for j in (1,2)] for i in (1,2)])
    return np.linalg.solve(cd,offset)+[float(h['CRPIX1'])-1,float(h['CRPIX2'])-1]


def keys_by_band(tile):
    keys={}
    for inst in ('NISP','VIS'):
        url='https://nasa-irsa-euclid-q1.s3.us-east-1.amazonaws.com/?'+urlencode({'list-type':2,'prefix':f'q1/MER/{tile}/{inst}/'})
        for el in ET.fromstring(read(url)).findall('.//{*}Key'):
            k=el.text
            for band in BANDS:
                suffix='VIS' if band=='VIS' else 'NIR-'+band
                for prod,term in [('image','BGSUB-MOSAIC-'+suffix+'_'),('rms','MOSAIC-'+suffix+'-RMS_'),('flag','MOSAIC-'+suffix+'-FLAG_')]:
                    if term in k:keys[band,prod]=k
    assert len(keys)==12,keys
    return keys


def get(t,key,cache):
    size=40 if t['radius_arcsec']>=2 else 24
    url='https://irsa.ipac.caltech.edu/ibe/data/euclid/'+key+f"?center={t['ra']:.8f},{t['dec']:.8f}deg&size={size}arcsec&gzip=false"
    file=cache/(t['id']+('_wide_' if size==40 else '_')+key.split('/')[-1])
    if not file.exists():
        file.write_bytes(read(url))
    raw=file.read_bytes()
    # The IRSA flag product stores signed integer bit patterns, sometimes with
    # a FITS unsigned offset. Read physical values to avoid discarding clean pixels.
    header,arr=read_primary(raw)
    receipt=dict(id=t['id'],product=key,url=url,bytes=len(raw),sha256=digest(raw))
    return header,arr,receipt


def read_primary(raw):
    h={}
    for pos in range(0,len(raw),80):
        c=raw[pos:pos+80].decode('ascii')
        if c[:8].strip()=='END':
            break
        if c[8:10]=='= ':
            val=c[10:].strip()
            h[c[:8].strip()]=val.split("'")[1].strip() if val.startswith("'") else val.split('/')[0].strip()
    off=(pos+80+2879)//2880*2880
    assert h['SIMPLE']=='T' and h['NAXIS']=='2'
    dt={'-32':'>f4','-64':'>f8','32':'>i4','16':'>i2','8':'u1'}[h['BITPIX']]
    shape=int(h['NAXIS2']),int(h['NAXIS1'])
    a=np.frombuffer(raw,dtype=dt,count=math.prod(shape),offset=off).reshape(shape)
    return h,a.astype('f8')*float(h.get('BSCALE',1))+float(h.get('BZERO',0))


def phot(a,rms,mask,cx,cy,r):
    q=int(np.ceil(r))+1
    yy,xx=np.mgrid[-q:q+1,-q:q+1]
    sub=(yy-(cy-round(cy)))**2+(xx-(cx-round(cx)))**2<=r*r
    yy,xx=yy[sub],xx[sub]
    def total(x,y): return a[y+yy,x+xx].sum()
    def good(x,y): return mask[y+yy,x+xx].all()
    core_ok=good(round(cx),round(cy))
    step=max(8,int(2*r)+2)
    blank=[total(x,y) for x in range(q+1,a.shape[1]-q,step) for y in range(q+1,a.shape[0]-q,step)
           if np.hypot(x-cx,y-cy)>max(50,3*r) and good(x,y)]
    if len(blank)<12 or not core_ok:
        return np.nan,np.nan,len(blank)
    bg,sc=sigma(np.asarray(blank))
    # Native RMS diagonal gives a lower bound on uncertainty for positively
    # correlated mosaics. Blanks carry resampling and confusion in the local field.
    err=max(sc,np.linalg.norm(rms[round(cy)+yy,round(cx)+xx]))
    return total(round(cx),round(cy))-bg,err,len(blank)


def analyze(t,files):
    name=t['id'];R=t['radius_arcsec']*10
    band={}
    for b in BANDS:
        h,im,_=files[b,'image']
        hr,rms,_=files[b,'rms']
        hf,flag,_=files[b,'flag']
        x,y=pix(h,t['ra'],t['dec'])
        assert np.hypot(*(pix(hr,t['ra'],t['dec'])-[x,y]))<.02
        assert np.hypot(*(pix(hf,t['ra'],t['dec'])-[x,y]))<.02
        ok=np.isfinite(im+rms) & (rms>0) & (rms<1e8)
        # MER propagates input bit flags by OR. A nonzero coadd flag can coexist
        # with good contributing exposures; retain per-aperture diagnostics.
        smooth=ndi.gaussian_filter(np.where(ok,im,0),math.sqrt(5**2-WIDTH[b]**2)/2.355)
        band[b]=h,im,rms,ok,smooth,x,y
    h,im,rms,ok,sm,cx,cy=band['VIS']
    rgrid=np.hypot(*[v-c for v,c in zip(np.indices(im.shape),[cy,cx])])
    sm=ndi.gaussian_filter(np.where(ok,im,0),1.3)
    base,noise=sigma(sm[ok&(rgrid>=45)&(rgrid<95)])
    peaks=np.argwhere((sm==ndi.maximum_filter(sm,7)) & (sm>base+4*noise) & ok &
                      (rgrid >= .65*R) & (rgrid < 1.4*R))
    peaks=sorted(peaks,key=lambda yx:-sm[tuple(yx)])
    keep=[]
    for y,x in peaks:
        if all(np.hypot(y-yy,x-xx)>=5 for yy,xx in keep):
            keep.append((int(y),int(x)))
    points=[('center',0,0)]+[(f'k{j+1}',x-cx,y-cy) for j,(y,x) in enumerate(keep)]
    measurements=[];colors={}
    for label,dx,dy in points:
        for b,(h,im,rms,ok,smooth,x,y) in band.items():
            for typ,r in ([('native_core',3.),('native_E',R)] if label=='center' else [('matched_035',3.5),('matched_050',5.)]):
                native=typ.startswith('native')
                # Expand the mask for the support of the smoothing kernel on colours.
                valid=ok if native else ndi.minimum_filter(ok.astype('u1'),size=13).astype(bool)
                f,e,n=phot(im if native else smooth,rms,valid,x+dx,y+dy,r)
                scale=1e6*3631*10**(-.4*float(h['MAGZERO']))
                f,e=f*scale,e*scale
                row=dict(id=name,point=label,dx_arcsec=round(dx*.1,4),dy_arcsec=round(dy*.1,4),band=b,method=typ,
                         radius_arcsec=round(.1*r,4),n_sky=n,flux_uJy=round(f,5),sigma_uJy=round(e,5),
                         snr=round(f/e,2),flux_upper5_uJy=round(max(0,f)+5*e,5))
                iy,ix=np.ogrid[:im.shape[0],:im.shape[1]]
                ap=(ix-(x+dx))**2+(iy-(y+dy))**2<=r*r
                flags=files[b,'flag'][1][ap].astype('int64')
                row['flag_or']=int(np.bitwise_or.reduce(flags))
                row['flag_nonzero_fraction']=round(float(np.mean(flags!=0)),4)
                row['ra'],row['dec']=world(h,x+dx,y+dy)
                measurements.append(row)
                if not native:colors[label,b,typ]=f,e
    pairs=[]
    for a,b in itertools.combinations([v[0] for v in points[1:]],2):
        ia,ib=int(a[1:])-1,int(b[1:])-1
        ang=math.degrees(math.acos(np.clip(np.dot(np.array(keep[ia])-[cy,cx],np.array(keep[ib])-[cy,cx])/
                           (rgrid[tuple(keep[ia])]*rgrid[tuple(keep[ib])]),-1,1)))
        for typ in ('matched_035','matched_050'):
            va=np.array([colors[a,k,typ] for k in BANDS]);vb=np.array([colors[b,k,typ] for k in BANDS])
            valid=np.isfinite(va+vb).all(axis=1);va,vb=va[valid],vb[valid]
            n=len(va)
            if n<2:continue
            def chi(logratio):
                ratio=np.exp(logratio)
                var=va[:,1]**2+(vb[:,1]*ratio)**2 + .05**2*(va[:,0]**2+(vb[:,0]*ratio)**2)
                return np.sum((va[:,0]-ratio*vb[:,0])**2/var)
            fit=optimize.minimize_scalar(chi,bounds=(-5,5),method='bounded')
            pairs.append(dict(id=name,a=a,b=b,opening_angle=round(ang,2),method=typ,nband=n,
                chi2=round(fit.fun,3),fitted_scale=round(np.exp(fit.x),4),approx_p=round(stats.chi2.sf(fit.fun,n-1),7)))
    board=Image.new('RGB',(5*151,174),'white');draw=ImageDraw.Draw(board)
    for j,b in enumerate(('VIS','Y','J','H','VIS')):
        _,a,_,ok,_,x,y=band[b]
        size=44
        stamp=a[round(y)-size:round(y)+size+1,round(x)-size:round(x)+size+1]
        base,noise=sigma(stamp.ravel())
        scaled=np.clip((np.arcsinh((stamp-base)/(2*noise))+1)/5,0,1)
        board.paste(Image.fromarray(np.uint8(255*scaled[::-1])).resize((151,151)),(151*j,0))
        draw.text((151*j+3,154),b if j!=4 else 'VIS positions',fill='black')
        if j==4:
            for label,dx,dy in points:
                xx=(size+.5+dx+x-round(x))*151/(2*size+1)+151*j
                yy=(size+.5-dy-y+round(y))*151/(2*size+1)
                col='orange' if label=='center' else 'cyan'
                draw.ellipse((xx-3,yy-3,xx+3,yy+3),outline=col)
                draw.text((xx+3,yy-8),'+' if label=='center' else label,fill=col)
    return measurements,pairs,board


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--cache',type=Path,required=True)
    p.add_argument('--out',type=Path,default=HERE/'validation')
    args=p.parse_args()
    args.out.mkdir(exist_ok=True);args.cache.mkdir(exist_ok=True,parents=True)
    ts=[v for f in sorted((HERE/'profiles').glob('*.json')) for v in json.loads(f.read_text())['profiles'] if v['pass_vis']]
    extra=json.loads((HERE/'unmodelled.json').read_text())
    if (HERE/'fallback.json').exists():
        extra+=json.loads((HERE/'fallback.json').read_text())
    for v in extra:
        for j,q in enumerate(v.get('positions',[])):
            if q['pass_vis']:
                tile=v.get('tile',v['id_str'].split('_')[0])
                ts.append({'id':tile+'_target_'+str(v['id_str'])+'_'+str(j),'ra':q['ra'],'dec':q['dec'],'radius_arcsec':q['r']})
    manifest=[];meas=[];pairs=[];keys={}
    for t in ts:
        tile=t['id'].split('_')[0]
        if tile not in keys:keys[tile]=keys_by_band(tile)
        with ThreadPoolExecutor(8) as ex:
            files=dict(zip(keys[tile],ex.map(lambda key:get(t,key,args.cache),keys[tile].values())))
        m,c,im=analyze(t,files)
        manifest.extend(d[2] for d in files.values());meas.extend(m);pairs.extend(c)
        im.save(args.out/(t['id']+'.png'))
        save_csv(args.out/'photometry.csv',meas);save_csv(args.out/'image_pairs.csv',pairs);save_csv(args.out/'native_inputs.csv',manifest)
        print(t['id'],len(m),len(c),flush=True)
if __name__=='__main__':
    main()
