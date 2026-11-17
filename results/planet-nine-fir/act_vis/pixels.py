"""Download, validate and register the four DR1 exposures in the ACT4_06 square."""
import csv,io,json,sys
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates
from common import CACHE,HALF,HERE,N,RA0,DE0,COS,T,MOTION,rows,world,coord
import vis_pixels as v

REC=CACHE/'receipts'
REC.mkdir(exist_ok=True)
v.CACHE=CACHE;v.REC=REC
HDU=v.HDU
KEYS=[x.text for x in ET.fromstring((HERE/'s3_listing.response').read_bytes()).iter() if x.tag.endswith('}Key')]


def get(stuff):
    i,row,prod,exp=stuff
    qc=next(k for k in KEYS if k.endswith('.json') and exp in k)
    centers=json.loads(v.read(v.BUCKET+qc,exp+'_qc'))
    yy,xx=np.indices((N,N),dtype=float);ra,de=world(xx-HALF,yy-HALF)
    data=np.zeros((N,N),'f4');weight=np.zeros_like(data);names=[]
    for key,ra0 in centers.items():
        if not key.endswith('.QUAD_RA'):continue
        d=centers[key[:-2]+'DEC']; x,y=coord(ra0,d)
        if max(abs(x),abs(y))>HALF+1500:continue
        quad=key.removesuffix('.QUAD_RA')
        a,b=map(int,quad.split('.')[0].split('-'))
        idx=((a-1)*6+b-1)*4+'EFGH'.index(quad[-1])
        base=6*v.BLOCK+idx*v.QUAD
        url=v.BUCKET+'dr1/VIS_QUAD/6578/'+prod
        stem=exp+'_'+quad
        h,hbytes=v.qheader(url,stem+'SCI',base,quad+'.SCI')
        assert h['BITPIX']==-32
        wcs=WCS(h)
        corner=np.array([[-HALF,-HALF],[HALF,-HALF],[HALF,HALF],[-HALF,HALF],[0,0]])
        u,z=wcs.all_world2pix(*world(*corner.T),0)
        if min(u)>v.NX-2 or max(u)<1 or min(z)>v.NY-2 or max(z)<1:continue
        lo=max(0,int(np.floor(min(z)))-8);hi=min(v.NY,int(np.ceil(max(z)))+9)
        if hi<=lo:continue
        # Interpolate TPV WCS only after verifying internal test points.
        from scipy.interpolate import RectBivariateSpline
        g=np.linspace(-HALF,HALF,17);xx0,yy0=np.meshgrid(g,g)
        us=wcs.all_world2pix(*world(xx0,yy0),0);spl=[RectBivariateSpline(g,g,a) for a in us]
        g2=(g[1:]+g[:-1])/2;test=wcs.all_world2pix(*world(*np.meshgrid(g2,g2)),0)
        err=max(abs(sp(g2,g2)-q).max() for sp,q in zip(spl,test))
        if err>.005:raise ValueError(f'WCS residual {err}')
        u,z=[sp(np.arange(-HALF,HALF+1),np.arange(-HALF,HALF+1)) for sp in spl]
        inside=(u>=2)&(u<v.NX-3)&(z>=max(2,lo+1))&(z<min(v.NY-3,hi-2))
        if not inside.any():continue
        planes=[]
        for j,(ext,dt) in enumerate([('SCI','>f4'),('RMS','>f4'),('FLG','>i4')]):
            pbase=base+j*HDU; hh,nbytes=v.qheader(url,stem+ext,pbase,quad+'.'+ext)
            assert hh['BITPIX']==(32 if ext=='FLG' else -32)
            b=v.read(url,stem+ext+'_rows',(pbase+nbytes+lo*v.NX*4,pbase+nbytes+hi*v.NX*4-1))
            planes.append(np.frombuffer(b,dtype=dt).reshape(hi-lo,v.NX))
        sci,rms,flag=planes;bg=float(np.nanmedian(sci[::12,::12][(flag[::12,::12]&v.BAD_MASK)==0]))
        ct=float(h['EXPTIME'])*10**(.4*(float(h['MAGZEROP'])-25))
        xy=[z[inside]-lo,u[inside]]
        ff=map_coordinates(((flag&v.BAD_MASK)!=0).astype('f4'),xy,order=1,prefilter=False)
        dd=map_coordinates(sci,xy,order=1,prefilter=False)
        rr=map_coordinates(rms,xy,order=1,prefilter=False)
        good=(ff==0)&np.isfinite(dd)&np.isfinite(rr)&(rr>0)
        iy,ix=np.where(inside)
        iy=iy[good];ix=ix[good]
        data[iy,ix]=(dd[good]-bg)/ct
        weight[iy,ix]=(ct/rr[good])**2
        names.append(dict(quadrant=quad,row_start=lo,row_stop=hi,wcs_error_pix=err,adu_per_vis25=ct,exp_e_per_vis25=ct*float(h['GAIN']),gain_e_per_adu=float(h['GAIN']),sky_adu=bg,zp=h['MAGZEROP'],exposure=h['EXPTIME'],apercorr=h.get('APERCORR'),unmasked_mosaic_pixels=len(iy)))
    np.save(CACHE/f'data_{i}.npy',data);np.save(CACHE/f'weight_{i}.npy',weight)
    print(i,exp,len(names),np.count_nonzero(weight),flush=True)
    return dict(index=i,product=prod,utc=T[i].isot,nominal_drift_pix=MOTION[i].tolist(),nominal_positive_weight=bool(weight[HALF,HALF]),quadrants=names)


def main():
    jobs=[]
    for i,row in enumerate(rows):
        prod=row['access_url'].split('file_name=')[-1]
        exp=prod.split('DET-')[1].split('-0000000')[0];jobs.append((i,row,prod,exp))
    with ThreadPoolExecutor(max_workers=4) as pool:
        info=list(pool.map(get,jobs))
    (HERE/'coverage.json').write_text(json.dumps(info,indent=2)+'\n')
    manifest={q.stem:json.loads(q.read_bytes()) for q in sorted(REC.glob('*.json'))}
    (HERE/'ranges.json').write_text(json.dumps(manifest,indent=1)+'\n')

if __name__=='__main__':main()
