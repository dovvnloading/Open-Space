"""Original CAL DQ and surface-brightness aperture photometry for COSMOS-022."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from io import BytesIO
from pathlib import Path
import json,sys
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
RA,DEC=150.122033,2.323386208

def measure(rec):
    stem=rec.exposure; visit=stem.split('_')[0]
    url=f'https://stpubdata.s3.amazonaws.com/jwst/public/{visit[:7]}/{visit}/{stem}_cal.fits'
    receipts=[]
    def read(start, size):
        body,r=byte_range(url,start,size);receipts.append(r)
        return body
    offset=0;out=dict(mosaic=rec.mosaic,exposure=stem,mjd=rec.mjd,band=rec.band,detector=stem.split('_')[-1])
    while True:
        buf=BytesIO(read(offset,57600));hdr=fits.Header.fromfile(buf);start=offset+buf.tell()
        dims=[hdr[f'NAXIS{k+1}'] for k in range(hdr['NAXIS'])]
        nbytes=(int(np.prod(dims)) if dims else 0)*(abs(hdr['BITPIX'])//8)+hdr.get('PCOUNT',0)
        offset=start+2880*((nbytes+2879)//2880);name=hdr.get('EXTNAME','PRIMARY')
        if name=='SCI':
            assert hdr['BUNIT']=='MJy/sr'
            w=WCS(hdr,naxis=2);xf,yf=map(float,w.world_to_pixel_values(RA,DEC));width,height=dims
            x0=max(0,int(xf)-30);x1=min(width,int(xf)+31);y0=max(0,int(yf)-30);y1=min(height,int(yf)+31)
            yy,xx=np.mgrid[y0:y1,x0:x1];a,d=w.pixel_to_world_values(xx,yy)
            radius=3600*np.hypot((a-RA)*np.cos(np.deg2rad(DEC)),d-DEC)
            pscale=np.sqrt(float(hdr['PIXAR_SR']))*206264.806
            out.update(xdet=xf,ydet=yf,edge_arcsec=pscale*min(xf,yf,width-1-xf,height-1-yf))
        if name in ('SCI','DQ'):
            r=read(start+y0*width*4,(y1-y0)*width*4)
            arr=np.frombuffer(r,dtype='>f4' if name=='SCI' else '>i4').reshape(-1,width)[:,x0:x1]
            if name=='SCI': sci=arr
            else:
                dq=(arr.astype('i8')*hdr.get('BSCALE',1)+hdr.get('BZERO',0)).astype('u4');break
    core=radius<.16;good=np.isfinite(sci)&((dq&1)==0)
    for label,bit in [('DO_NOT_USE',1),('SATURATED',2),('JUMP_DET',4),('DROPOUT',8),('OUTLIER',16),('PERSISTENCE',32),('HOT',2048),('DEAD',1024)]:
        out[label]=int(((dq[core]&bit)!=0).sum())
    out['core_pixels']=int(core.sum());out['finite_core']=int(np.isfinite(sci[core]).sum())
    for rad in (.16,.28):
        sky=np.median(sci[(radius>rad+.12)&(radius<rad+.32)&good])
        wt=np.clip((rad-radius)/pscale+.5,0,1)
        out[f'flux_r{int(rad*100):02d}_njy']=float(np.sum((sci[good]-sky)*wt[good])*(pscale/206264.806)**2*1e15)
        # Exclude every flagged pixel; report the weighted good-area fraction.
        clean=good&(dq==0);out[f'clean_r{int(rad*100):02d}_njy']=float(np.sum((sci[clean]-sky)*wt[clean])*(pscale/206264.806)**2*1e15)
        out[f'clean_r{int(rad*100):02d}_fraction']=float(wt[clean].sum()/wt.sum())
    return out,dict(exposure=stem,url=url,ranges=receipts)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('campaign',type=Path);p.add_argument('output',type=Path)
    args=p.parse_args();sys.path.insert(0,str(args.campaign.resolve()));from pixel_check import byte_range
    a=pd.read_csv(args.campaign/'individual_apertures.csv');t=a[a.candidate=='cosmos_medium-022']
    args.output.mkdir(exist_ok=True,parents=True)
    with ThreadPoolExecutor(6) as pool: rows=list(pool.map(measure,t.itertuples()))
    table=pd.DataFrame([r for r,_ in rows])
    table.to_csv(args.output/'native_cal.csv',index=False)
    table['epoch']=table.mosaic.str.split('_').str[0]
    table['dithers']=(table.exposure.str.split('_').str[2].astype(int)%2).map({0:'even',1:'odd'})
    table.groupby(['band','epoch','dithers']).flux_r16_njy.agg(count='size',median_njy='median',minimum_njy='min',maximum_njy='max').reset_index().to_csv(args.output/'native_splits.csv',index=False)
    (args.output/'native_receipts.json').write_text(json.dumps([r for _,r in rows],indent=1)+'\n')
