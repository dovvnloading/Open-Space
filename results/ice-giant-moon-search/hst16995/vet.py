"""Audit the strongest Neptune orbit, compare holdout F763M images and scrambled banks."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy import ndimage as ndi
from reduce import (AU, ARCSEC, ROOT, DOMAINS, NAMES, basis, filt, grid, orbit, peak_center, prep,
                    read_csv, samples, score, tangent_basis, write_csv)


def holdout(cache,best,epoch):
    rows=read_csv('neptune/f763m_inputs.csv'); out=[]
    for r in rows:
        name=r['obs_id']; fn=cache/(name+'_flt.fits')
        if not fn.exists():
            from urllib.request import urlopen
            with urlopen(r['url'],timeout=120) as response: fn.write_bytes(response.read())
        if hashlib.sha256(fn.read_bytes()).hexdigest()!=r['sha256']: raise ValueError(name)
        with fits.open(fn) as h:
            s=h['SCI'];im=s.data.copy();dq=h['DQ'].data.copy();w=WCS(s.header,fobj=h)
            t=(h[0].header['EXPSTART']+h[0].header['EXPEND'])/2
            if h[0].header['FILTER']!='F763M': raise ValueError(name)
        ra,dec,delta,lt=[float(r[k]) for k in ['ra899','dec899','delta','lt']]
        p=w.all_world2pix([[ra,dec]],0,quiet=True)[0]
        offsets=w.all_world2pix([[ra+1/np.cos(np.deg2rad(dec))/3600,dec],[ra,dec+1/3600]],0,quiet=True)-p
        jac=offsets.T*ARCSEC;_,east,north=tangent_basis(ra,dec);u,v,_=basis('Neptune',2400000.5+t)
        project=jac@np.array([[east@u,east@v],[north@u,north@v]])/(delta*AU)
        expected=np.pi*24764**2*abs(np.linalg.det(jac))/(delta*AU)**2
        limb=im>np.nanpercentile(im,99.98)*.3
        if not .6<limb.sum()/expected<1.4:
            out.append(dict(rootname=name,visit=name[4:6],status='excluded_disk_area',x='',y='',score='',control_score=''));continue
        yy,xx=np.where(limb);c=np.array([np.median(xx),np.median(yy)])
        prot=c+w.all_world2pix([[float(r['ra808']),float(r['dec808'])]],0,quiet=True)[0]-p
        det=peak_center(im,prot)
        if det is not None and np.linalg.norm(det-prot)<3:
            c+=det-prot;prot=det
        from geometry import frequency
        a,phase=best;theta=phase+frequency('Neptune',a)*(t-lt/1440-epoch)*86400
        pos=c+project@[a*np.cos(theta),a*np.sin(theta)]
        yy,xx=np.indices(im.shape);rr=np.hypot(xx-c[0],yy-c[1])
        good=(dq==0)&(rr>np.sqrt(expected/np.pi)+14)
        filtered=filt(im,good)
        def snr(xy):
            y,x=xy[::-1];zone=(abs(rr-np.linalg.norm(xy-c))<7)&good
            sky=filtered[zone]
            if len(sky)<100 or not(8<y<504 and 8<x<504): return None
            med=np.median(sky);scale=max(.7,1.4826*np.median(np.abs(sky-med)))
            if ndi.map_coordinates(good.astype('f4'),[[y],[x]],order=1)[0]<.75:return None
            return round(float((ndi.map_coordinates(filtered,[[y],[x]],order=1)[0]-med)/scale),3)
        out.append(dict(rootname=name,visit=name[4:6],status='tested',x=round(pos[0],3),y=round(pos[1],3),
                        score=snr(pos),control_score=snr(prot)))
    write_csv(ROOT/'neptune/f763m_predictions.csv',out)


def forensics(cache,frames,reg,bank):
    best=bank[:1]
    x,y=orbit(frames,'Neptune',best[:,0],best[:,1]);raw,ok=samples(frames,'Neptune',best)
    ep={r['rootname']:r for r in read_csv('ephemerides.csv') if r['body']=='899'}
    audit=[]
    for j,(f,r) in enumerate(zip(frames,[r for r in reg if r['status']=='used'],strict=True)):
        name=r['rootname'];pos=np.array([x[0,j],y[0,j]])
        src=peak_center(f['image'],pos)
        with fits.open(cache/(name+'_flt.fits')) as h:
            w=WCS(h['SCI'].header,fobj=h)
            center=w.all_pix2world([f['c']],0)[0]
            if src is not None: sky=w.all_pix2world([src],0)[0]-center+np.array([float(ep[name]['ra']),float(ep[name]['dec'])])
        audit.append(dict(rootname=name,visit=f['visit'],x=round(pos[0],4),y=round(pos[1],4),
                          frame_score=round(float(raw[0,j]),3),supported=int(ok[0,j]),peak_x=round(src[0],4) if src is not None else '',
                          peak_y=round(src[1],4) if src is not None else '',peak_ra=round(sky[0],8) if src is not None else '',
                          peak_dec=round(sky[1],8) if src is not None else ''))
    write_csv(ROOT/'neptune/f845m_forced.csv',audit)


def main(args):
    fs,reg,_,_=prep(args.cache,'Neptune')
    best=read_csv('neptune/highest_scores.csv')[0];best=np.array([[float(best['a_km']),float(best['phase_rad'])]])
    forensics(args.cache,fs,reg,best)
    if (ROOT/'neptune/f763m_inputs.csv').exists():
        holdout(args.cache,best[0],np.mean([f['mid'] for f in fs]))
    if not args.scrambles:return
    bank=grid('Neptune',DOMAINS['Neptune']); visits=np.array([f['visit'] for f in fs])
    unique=np.unique(visits);mid=np.array([f['mid'] for f in fs]);centers=np.array([mid[visits==v].mean() for v in unique])
    controls=[];rng=np.random.default_rng(78054)
    for trial in range(args.scrambles):
        perm=rng.permutation(len(unique))
        while np.any(perm==np.arange(len(unique))):perm=rng.permutation(len(unique))
        shifted=mid.copy()
        for i,v in enumerate(unique):shifted[visits==v]+=centers[perm[i]]-centers[i]
        ctrl=[dict(f,mid=t) for f,t in zip(fs,shifted)];peak=0.;cross=drop=0
        for start in range(0,len(bank),25000):
            raw,ok=samples(ctrl,'Neptune',bank[start:start+25000])
            r,_=score(raw,ok,visits);v,_=score(raw,ok,visits,True)
            peak=max(peak,float(r.max()));cross+=int((r>=1).sum());drop+=int((v>=1).sum())
        controls.append(dict(scramble=trial,visit_order=' '.join(unique[perm]),max_screen_rank=round(peak,5),
                             screen_crossings=cross,dropout_crossings=drop))
        write_csv(ROOT/'neptune/scrambles.csv',controls)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('cache',type=Path)
    ap.add_argument('--scrambles',type=int,default=8);main(ap.parse_args())
