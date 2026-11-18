"""Inspect the released Euclid and HST image planes at additional large-source leads."""
import csv
import hashlib
import io
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import gaussian_filter, maximum_filter, median_filter

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'other'))
from euclid_scene import CACHE, clusters, coords, measure, repeated_pairs


def hst_image(job):
    event,ra,dec,row=job
    obs=row['observation']
    url='https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/'+obs+'_flt.fits'
    receipt={'event':event,'observation':obs,'mjd':row['start_mjd'],'url':url}
    try:
        with urlopen(url,timeout=40) as reply: data=reply.read()
        receipt.update(sha256=hashlib.sha256(data).hexdigest(),bytes=len(data))
        with fits.open(io.BytesIO(data)) as hdus:
            science=hdus['SCI'];w=WCS(science.header)
            x,y=w.all_world2pix(ra,dec,0)
            nx,ny=science.header['NAXIS1'],science.header['NAXIS2']
            receipt.update(obstype=hdus[0].header.get('OBSTYPE',''),exptime=hdus[0].header.get('TEXPTIME'),
                           target_xy_pixels=[float(x),float(y)],image_shape=[ny,nx])
            if not (0<=x<nx and 0<=y<ny):raise ValueError('position outside image')
            left=int(np.clip(round(float(x))-50,0,nx-101))
            down=int(np.clip(round(float(y))-50,0,ny-101))
            im=science.data[down:down+101,left:left+101].astype(float)
            good=(hdus['DQ'].data[down:down+101,left:left+101]==0)&np.isfinite(im)
            smooth=gaussian_filter(im,.85)
            contrast=smooth-median_filter(smooth,13)
            med=np.median(contrast[good])
            scatter=1.4826*np.median(abs(contrast[good]-med))
            ys,xs=np.where(good&(contrast>med+7*scatter)&(maximum_filter(contrast,7)==contrast))
            peaks=[]
            for yy,xx in zip(ys,xs):
                if min(yy,xx,100-yy,100-xx)<4:continue
                mask=good[yy-2:yy+3,xx-2:xx+3]
                cut=im[yy-2:yy+3,xx-2:xx+3]
                if mask.sum()<18:continue
                weight=np.where(mask,np.maximum(0.,cut-np.median(cut[mask])),0.)
                if not weight.sum():continue
                cy,cx=np.indices((5,5))
                px=left+xx-2+np.average(cx,weights=weight);py=down+yy-2+np.average(cy,weights=weight)
                a,d=w.all_pix2world(px,py,0)
                dx=(float(a)-ra)*np.cos(np.deg2rad(dec))*3600;dy=(float(d)-dec)*3600
                sep=np.hypot(dx,dy)
                if sep<=2:peaks.append(dict(east_arcsec=round(dx,4),north_arcsec=round(dy,4),
                                            separation_arcsec=round(sep,4)))
            receipt.update(status='ok',valid_pixels=int(good.sum()),total_pixels=im.size,
                           target_dq_clear=bool(good[int(round(float(y)))-down,int(round(float(x)))-left]),
                           distance_to_detector_edge_pix=round(min(float(x),nx-1-float(x),float(y),ny-1-float(y)),2),
                           peaks=sorted(peaks,key=lambda p:p['separation_arcsec']))
    except Exception as err:
        receipt.update(status='error',error=str(err))
    return receipt


if __name__=='__main__':
    CACHE.mkdir(parents=True,exist_ok=True)
    targets={}
    with (HERE/'kmt/catalog-2024.tsv').open() as f:
        for r in csv.DictReader(f,delimiter='\t'):
            targets[r['Event']]=coords(r['RA'],r['Dec'])
    with (HERE/'ogle/catalog.csv').open() as f:
        for r in csv.DictReader(f):targets['OGLE-'+r['event']]=coords(r['ra_hms'],r['dec_dms'])
    blocks=json.loads((HERE/'large-source-archive-queries.json').read_text())
    euclid=[];hst=[];names=set()
    for block in blocks:
        q=block['query'];event=q['event'];ra,dec=targets[event]
        if q['service']=='IRSA':
            data={r['access_url']:r for r in block['products'] if '-DET-' in (r.get('access_url') or '')}
            series=sorted(data.values(),key=lambda r:float(r['start_mjd']))
            if series:
                names.add(event)
                indices=sorted(set(np.linspace(0,len(series)-1,min(7,len(series)),dtype=int)))
                euclid.extend((event,ra,dec,series[i]) for i in indices)
        else:
            hst.extend((event,ra,dec,r) for r in block['products'] if r['collection']=='HST')
    with ThreadPoolExecutor(4) as pool:
        eu=list(pool.map(measure,euclid));hs=list(pool.map(hst_image,hst))
    det=[r for _,rows,_ in eu for r in rows]
    sources=clusters(det,sorted(names))
    result=dict(euclid_queries=[q for q,_,_ in eu],euclid_sources=sources,
                euclid_pairs=repeated_pairs(det,sources,sorted(names)),
                euclid_detections=det,euclid_range_receipts=[r for _,_,rows in eu for r in rows],
                hst_exposures=hs)
    (HERE/'large-source-images.json').write_text(json.dumps(result,indent=2)+'\n')
