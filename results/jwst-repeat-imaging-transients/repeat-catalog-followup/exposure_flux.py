"""Verify isolated COSMOS coincidences in the actual constituent exposures."""
from concurrent.futures import ThreadPoolExecutor
from io import BytesIO
import json
import argparse
from pathlib import Path

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import pandas as pd

from campaign import HERE,product_url
from pixel_check import aperture,byte_range,stamp,XX,YY,STEP


def image_context(product):
    url=product_url(product,'i2d.fits');offset=0;extensions={};receipts=[]
    while True:
        raw,r=byte_range(url,offset,57600);receipts.append(r)
        while not any(raw[j:j+8]==b'END     ' for j in range(0,len(raw),80)):
            b,r=byte_range(url,offset+len(raw),57600);raw+=b;receipts.append(r)
        reader=BytesIO(raw);header=fits.Header.fromfile(reader);begin=offset+reader.tell()
        dims=tuple(header[f'NAXIS{i+1}'] for i in range(header['NAXIS']))
        size=abs(header['BITPIX'])//8*(int(np.prod(dims)) if dims else 0)+header.get('PCOUNT',0)
        name=header.get('EXTNAME','PRIMARY')
        if name in ['SCI','CON']:extensions[name]=(header,begin)
        if name=='HDRTAB':
            parts=[]
            for i in range(0,size,1024**2):
                b,r=byte_range(url,begin+i,min(size-i,1024**2));parts.append(b);receipts.append(r)
            f=BytesIO(raw[:reader.tell()]+b''.join(parts)+bytes((-size)%2880))
            table=fits.BinTableHDU.readfrom(f).data.copy()
            return url,extensions,table,receipts
        if name=='ASDF':raise ValueError('missing context table')
        offset=begin+2880*((size+2879)//2880)


def covering(meta,ra,dec):
    url,extensions,tab,_=meta;hdr,begin=extensions['CON']
    w=WCS(extensions['SCI'][0],naxis=2)
    x,y=[int(np.rint(float(q))) for q in w.world_to_pixel_values(ra,dec)]
    nx,ny=hdr['NAXIS1'],hdr['NAXIS2']
    if not (0<=x<nx and 0<=y<ny):return [],[]
    names=[];receipts=[]
    for plane in range(hdr.get('NAXIS3',1)):
        b,r=byte_range(url,begin+4*((plane*ny+y)*nx+x),4);receipts.append(r)
        value=int.from_bytes(b,'big')
        for bit in range(32):
            if value & (1<<bit):
                row=tab[32*plane+bit];names.append((str(row['FILENAME']).strip().split('_o')[0],float(row['MJD-BEG']),int(row['PATT_NUM'])))
    return names,receipts


def photometry(task):
    hit,band,product,gain,stem,mjd,dither=task
    visit=stem.split('_')[0]
    url=f'https://stpubdata.s3.amazonaws.com/jwst/public/{visit[:7]}/{visit}/{stem}_i2d.fits'
    raw,r=byte_range(url,0,57600);reader=BytesIO(raw);fits.Header.fromfile(reader);header=fits.Header.fromfile(reader)
    im,recs=stamp((url,header,reader.tell(),r),hit.ra,hit.dec)
    # Report fixed-WCS photometry and a separately centered aperture. Track every
    # shift, polarity and fractional coverage to expose registration/edge effects.
    valid=np.isfinite(im);core=np.hypot(XX,YY)<.22
    ann=(np.hypot(XX,YY)>.36)&(np.hypot(XX,YY)<.64)
    sx=sy=np.nan
    if valid[ann].any() and valid[core].all():
        weights=np.where(core&valid,np.maximum(im-np.nanmedian(im[ann]),0),0)
        if weights.sum()>0:
            sx=np.sum(XX*weights)/weights.sum();sy=np.sum(YY*weights)/weights.sum()
    row=dict(candidate=hit.candidate,band=band,mosaic=product,exposure=stem,mjd=mjd,dither=dither,
             coverage=float(valid[core].mean()),east_arcsec=sx,north_arcsec=sy,
             fixed_njy=aperture(im,0,0,.16)*gain,centered_njy=aperture(im,sx,sy,.16)*gain if np.isfinite(sx) else np.nan,
             wide_njy=aperture(im,0,0,.28)*gain)
    return row,dict(url=url,ranges=[r]+recs)


def main():
    hits=pd.read_csv(HERE/'coincidences.csv').set_index('candidate',drop=False)
    pairs=pd.read_csv(HERE/'pairs.csv').merge(pd.read_csv(HERE/'screen_counts.csv'),on=['pair','band'])
    tests={'cosmos_medium-012':['f140m','f210m'],'cosmos_medium-022':['f162m','f182m'],
           'cosmos_medium-010':['f162m','f182m'],'cosmos_medium-014':['f162m','f210m']}
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bright-source',action='store_true',help='Check the bright, unconstrained COSMOS-007 source separately')
    args=parser.parse_args()
    if args.bright_source:
        tests={'cosmos_medium-007':['f140m','f162m']}
    prefix='individual_bright' if args.bright_source else 'individual'
    plans=[]
    for key,bands in tests.items():
        h=hits.loc[key]
        for band in bands:
            row=pairs[(pairs.pair==h.pair)&(pairs.band==band)].iloc[0]
            plans.extend([(h,band,row.reference,row.flux_scale),(h,band,row.search,1.)])
    products=sorted({p[2] for p in plans})
    with ThreadPoolExecutor(8) as pool:
        meta=dict(zip(products,pool.map(image_context,products)))
        contexts=list(pool.map(lambda t:covering(meta[t[2]],t[0].ra,t[0].dec),plans))
        tasks=[(*p,*exp) for p,(names,_) in zip(plans,contexts) for exp in names]
        results=list(pool.map(photometry,tasks))
    pd.DataFrame([r for r,_ in results]).to_csv(HERE/(prefix+'_apertures.csv'),index=False,float_format='%.10g')
    receipts=[dict(url=meta[p][0],ranges=meta[p][3]) for p in products]
    receipts.extend(dict(url=meta[t[2]][0],candidate=t[0].candidate,ranges=rs) for t,(_,rs) in zip(plans,contexts))
    receipts.extend(r for _,r in results)
    (HERE/(prefix+'_inputs.json')).write_text(json.dumps(receipts,separators=(',',':'))+'\n')
    print(len(results),'constituent exposure apertures')


if __name__=='__main__':main()
