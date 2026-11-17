"""Two-filter El Gordo screen. Input: locally fetched i2d windows (MJy/sr)."""
import argparse
import csv
from pathlib import Path
import json
import warnings

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy import ndimage as nd
from scipy.signal import fftconvolve
from scipy.spatial import cKDTree

HERE=Path(__file__).resolve().parent
CENTER=(15.735,-49.252)
SIDE=1101
SCALE=0.06
ARCSEC2_SR=(np.pi/(180*3600))**2
BANDS=('f150w','f444w')
RNG=534566
ALERT=5.
WIN=8

def write_csv(file,rows,fields):
    with (HERE/file).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)

def mad(a):
    return 1.4826*np.median(np.abs(a-np.median(a)))

def getgrid():
    w=WCS(naxis=2);w.wcs.ctype=['RA---TAN','DEC--TAN']
    w.wcs.crpix=[(SIDE+1)/2]*2;w.wcs.crval=CENTER;w.wcs.cdelt=[-SCALE/3600,SCALE/3600]
    return w

def render(path,w,aff=None):
    with fits.open(path) as h:
        assert h['SCI'].header['BUNIT']=='MJy/sr'
        nw=WCS(h['SCI'].header,naxis=2)
        yy,xx=np.indices((SIDE,SIDE),dtype=float)
        if aff is not None:
            dx,dy=aff@np.array([xx.flatten()-SIDE/2,yy.flatten()-SIDE/2,np.ones(SIDE*SIDE)])
            xx+=dx.reshape(xx.shape);yy+=dy.reshape(yy.shape)
        sx,sy=nw.world_to_pixel_values(*w.pixel_to_world_values(xx,yy))
        g=[sy,sx]
        sci=nd.map_coordinates(h['SCI'].data,g,order=1,mode='constant',cval=np.nan)
        err=nd.map_coordinates(h['ERR'].data,g,order=1,mode='constant',cval=np.nan)
        weight=nd.map_coordinates(h['WHT'].data,g,order=0,mode='constant',cval=0)
        con=h['CON'].data
        if con.ndim==2: con=con[None]
        # A context-plane bit stands for a contributing exposure.
        n=np.sum([np.bitwise_count(np.asarray(p,dtype=np.uint32)) for p in con],axis=0)
        n=nd.map_coordinates(n.astype(float),g,order=0,mode='constant',cval=0)
        good=(weight>0)&np.isfinite(sci)&(err>0)&np.isfinite(err)&(n>=2)
        img=np.where(good,sci,0.)
        return img,good,dict(product=path.stem,date=h[0].header['DATE-BEG'],mjd=h['SCI'].header['MJD-BEG'],area=float(good.sum()*SCALE**2/3600),median_con=float(np.median(n[good])))

def sources(img,good):
    enc=nd.gaussian_filter(img,1.1)-nd.gaussian_filter(img,4)
    sig=mad(enc[good])
    m=nd.binary_erosion(good,iterations=15)&(enc>8*sig)&(enc==nd.maximum_filter(enc,11))
    y,x=np.where(m);order=np.argsort(enc[y,x])[::-1][:700];centers=[]
    for a,b in zip(x[order],y[order]):
        c,d=float(a),float(b)
        for _ in range(3):
            yy,xx=np.mgrid[b-4:b+5,a-4:a+5];weight=np.maximum(img[yy,xx]-np.median(img[b-6:b+7,a-6:a+7]),0)*np.exp(-((xx-c)**2+(yy-d)**2)/6)
            if weight.sum()==0: break
            c,d=(np.sum(weight*t)/np.sum(weight) for t in (xx,yy))
        centers.append([c,d])
    return np.array(centers).reshape(-1,2)

def affine_fit(ref,test):
    dist,ix=cKDTree(test).query(ref);hits=dist<12
    a,b=ref[hits],test[ix[hits]]
    assert len(a)>=8
    design=np.c_[a-SIDE/2,np.ones(len(a))]
    delta=b-a
    use=np.linalg.norm(delta-np.median(delta,axis=0),axis=1)<3
    for _ in range(4):
        fit=np.linalg.lstsq(design[use],delta[use],rcond=None)[0]
        res=np.linalg.norm(delta-design@fit,axis=1)
        use=res<max(.28,min(2.,np.median(res[use])+3*mad(res[use])))
    return fit.T,int(sum(use)),float(np.sqrt(np.mean(res[use]**2)))

def psf(img,good,stars,band):
    kernel=[]; widths=[]
    for x,y in stars:
        xi,yi=round(x),round(y)
        if not good[yi-WIN:yi+WIN+1,xi-WIN:xi+WIN+1].all(): continue
        yy,xx=np.mgrid[-WIN:WIN+1,-WIN:WIN+1]
        a=nd.map_coordinates(img,[y+yy,x+xx],order=1)
        radius=np.hypot(xx,yy)
        a-=np.median(a[radius>=7]);noise=mad(a[radius>=7]);total=a[radius<5].sum()
        if total<=0 or a[WIN,WIN]<20*noise:continue
        ratio=a[radius<=1.5].sum()/total
        if ratio<.2 or ratio>.98:continue
        kernel.append(a/total);widths.append(ratio)
    if len(kernel)<4: raise ValueError(f'Insufficient PSF samples for {band}: {len(kernel)}')
    # Use the most concentrated compact objects selected in both images. Trim low
    # S/N wings and high outliers through pixelwise median coaddition.
    inds=np.argsort(widths)[-min(12,max(4,len(widths)//5)):]
    p=np.median(np.stack(kernel)[inds],axis=0)
    p=np.maximum(p,0)*np.clip((8-radius)/2,0,1)
    return p/p.sum(),dict(compact=len(kernel),used=len(inds),core_ratio=float(np.median(np.array(widths)[inds])))

def convolve(a,p): return fftconvolve(a,p,mode='same')

def local_noise(img,good):
    tile=60;ny=int(np.ceil(SIDE/tile))
    coarse=np.empty((ny,ny))
    for i,j in np.ndindex(coarse.shape):
        sl=np.s_[i*tile:(i+1)*tile,j*tile:(j+1)*tile]
        v=img[sl][good[sl]]
        coarse[i,j]=mad(v) if len(v)>180 else np.nan
    globalnoise=mad(img[good]);coarse=np.where(np.isfinite(coarse),coarse,globalnoise)
    y,x=np.indices(img.shape)/tile-.5
    return nd.map_coordinates(coarse,[y,x],order=1,mode='nearest')

def additions(shape,positions,amps):
    im=np.zeros(shape)
    for (cx,cy),amp in zip(positions,amps):
        j,i=np.floor([cx,cy]).astype(int); fx,fy=cx-j,cy-i
        im[i:i+2,j:j+2]+=amp*np.outer([1-fy,fy],[1-fx,fx])
    return im

def response(diff,effective):
    # Remove diffuse structure before smoothing by the common point-source template.
    filtered=diff-nd.gaussian_filter(diff,6)
    return convolve(filtered,effective[::-1,::-1])

def signed_peaks(score,good):
    out=[]
    for sign in (1,-1):
        u=sign*score
        y,x=np.where(good&(u>=ALERT)&(nd.maximum_filter(u,size=7)==u))
        out.extend((int(i),int(j),int(sign),float(u[j,i])) for j,i in zip(y,x))
    return out

def analyze(band,paths,w):
    im0,m0,meta0=render(paths[0],w)
    raw,m1,_=render(paths[1],w)
    c0,c1=sources(im0,m0),sources(raw,m1)
    aff,nanchor,rms=affine_fit(c0,c1)
    im1,m1,meta1=render(paths[1],w,aff)
    good=m0&m1
    stars0,stars1=sources(im0,good),sources(im1,good)
    dist,ix=cKDTree(stars1).query(stars0)
    stars0,stars1=stars0[dist<.65],stars1[ix[dist<.65]]
    p0,pmeta0=psf(im0,good,stars0,band);p1,pmeta1=psf(im1,good,stars1,band)
    common=fftconvolve(p0,p1)
    epoch0,epoch1=convolve(im0,p1),convolve(im1,p0)
    # Bright matched-source residuals are inspected separately. Exclude only the
    # upper percentile of scene brightness plus a PSF wing margin.
    bright=np.maximum(epoch0-np.median(epoch0[good]),epoch1-np.median(epoch1[good]))
    threshold=np.percentile(bright[good],99)
    brightmask=nd.maximum_filter(bright>threshold,size=21)
    good=nd.binary_erosion(good,iterations=24)&~brightmask
    d=epoch1-epoch0
    amp=response(d,common)
    amp-=np.median(amp[good]);noise=local_noise(amp,good)
    score=amp/noise
    rows=signed_peaks(score,good)
    # The matched filter output for a unit integrated flux distributed as the
    # empirical PSF after cross-convolution.
    pulse=np.zeros((81,81));pulse[40,40]=1
    resp=float(response(convolve(pulse,common),common)[40,40])
    to_njy=1e15*ARCSEC2_SR*SCALE**2
    meta=dict(epoch=[meta0,meta1],anchors=nanchor,astrom_rms_arcsec=SCALE*rms,affine=aff.tolist(),
              psf=[pmeta0,pmeta1],bright_mjysr=float(threshold),common_area_arcmin2=float(good.sum()*SCALE**2/3600),
              peaks_5sigma=len(rows),median_formal_5sigma_njy=float(5*np.median(noise[good])*to_njy/resp))
    print(band,meta,flush=True)
    return dict(good=good,score=score,rows=rows,noise=noise,psfs=[p0,p1],effective=common,resp=resp,
                to_njy=to_njy,meta=meta,im0=im0,im1=im1,diff=d)

def joint(rows,bands,w):
    reports=[]
    data0,data1=bands
    for x,y,sign,snr in rows:
        r=data1['rows']
        match=[rj for rj in r if rj[2]==sign and np.hypot(rj[0]-x,rj[1]-y)<=2.5]
        if not match:continue
        j=max(match,key=lambda rj:rj[3])
        ra,dec=w.pixel_to_world_values(x,y)
        reports.append(dict(ra=float(ra),dec=float(dec),x=x,y=y,sign=sign,snr_f150w=snr,snr_f444w=j[3],
                            dx_f444w=j[0]-x,dy_f444w=j[1]-y,
                            delta_njy_f150w=sign*snr*data0['noise'][y,x]*data0['to_njy']/data0['resp'],
                            delta_njy_f444w=sign*j[3]*data1['noise'][j[1],j[0]]*data1['to_njy']/data1['resp']))
    return reports

def select_positions(good,rng,nsample):
    # Spatially independent trials away from boundaries, separated by >36 output pixels.
    yy,xx=np.where(nd.binary_erosion(good,iterations=4))
    inds=rng.permutation(len(xx));positions=[]
    for k in inds:
        p=np.array([xx[k],yy[k]],dtype=float)+rng.uniform(-.5,.5,2)
        if all(np.linalg.norm(p-p0)>36 for p0 in positions):positions.append(p)
        if len(positions)==nsample:return np.array(positions)
    return np.array(positions)

def inject(bands):
    rng=np.random.default_rng(RNG)
    overlap=bands[0]['good']&bands[1]['good']
    scales=np.arange(23,30.01,.5)
    trials=[]
    for mag in scales:
        jy_nano=10**((31.4-mag)/2.5)
        for sign in (1,-1):
            positions=select_positions(overlap,rng,90)
            stats=[]
            for d in bands:
                px=additions(overlap.shape,positions,np.full(len(positions),sign*jy_nano/d['to_njy']))
                # Add a source into the registered native epoch, pass it through its
                # own and the opposite epoch PSF, then the full same detection filter.
                model=convolve(convolve(px,d['psfs'][int(sign==1)]),d['psfs'][int(sign==-1)])
                score=d['score']+response(model,d['effective'])/d['noise']
                detected=np.array([(x,y) for x,y,s,n in signed_peaks(score,d['good']) if s==sign]).reshape(-1,2)
                stats.append(detected)
            count=0
            for px in positions:
                a,b=[points[np.linalg.norm(points-px,axis=1)<=2.5] for points in stats]
                if len(a) and len(b) and np.min(np.linalg.norm(a[:,None]-b[None,:],axis=-1))<=2.5: count+=1
            trials.append(dict(ab_each=float(mag),polarity=sign,injected=len(positions),recovered=count))
        print("injections",mag,[r['recovered'] for r in trials[-2:]],flush=True)
    return trials

def main():
    parser=argparse.ArgumentParser();parser.add_argument('cache',type=Path);args=parser.parse_args()
    w=getgrid();bands=[]
    for b in BANDS:
        paths=sorted(args.cache.glob('*'+b+'.fits'));assert len(paths)==2
        bands.append(analyze(b,paths,w))
    overlap=bands[0]['good']&bands[1]['good']
    detections=joint(bands[0]['rows'],bands,w)
    fields=['ra','dec','x','y','sign','snr_f150w','snr_f444w','dx_f444w','dy_f444w','delta_njy_f150w','delta_njy_f444w']
    write_csv('joint_detections.csv',detections,fields)
    peaks=[]
    for b,d in zip(BANDS,bands):
        for x,y,s,n in d['rows']:
            ra,dec=w.pixel_to_world_values(x,y)
            peaks.append(dict(band=b,x=x,y=y,sign=s,snr=n,ra=float(ra),dec=float(dec),other_band_mask_valid=bool(overlap[y,x])))
    write_csv('single_band_peaks.csv',peaks,['band','x','y','sign','snr','ra','dec','other_band_mask_valid'])
    meta=dict(bands={b:d['meta'] for b,d in zip(BANDS,bands)},area_joint_arcmin2=float(overlap.sum()*SCALE**2/3600),two_band_detections=len(detections),grid_size=SIDE,grid_arcsec=SCALE,center_ra_dec=CENTER)
    (HERE/'screen.json').write_text(json.dumps(meta,indent=2)+'\n')
    print('joint',len(detections),flush=True)
    trials=inject(bands)
    write_csv('injections.csv',trials,['ab_each','polarity','injected','recovered'])
if __name__=='__main__': main()
