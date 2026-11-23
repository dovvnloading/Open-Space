"""Post-inspection F277W control diagnostic for the February GOODS-S epoch.

Uses the existing motion_dithers.py compact/isolated SNR>28 rule. Select eight
nearest controls within 60 arcsec from the October catalog before examining
February fluxes. A small control sample cannot calibrate source-specific PSFs.
"""
import argparse
from datetime import datetime, timezone
import hashlib
from io import BytesIO
import json
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates

BASE="https://stpubdata.s3.amazonaws.com/"
CAT="jwst/public/jw12577/L3/t/o001/jw12577-o001_t001_nircam_clear-f277w_cat.ecsv"
PRODUCTS={"february":"jwst/public/jw05324/L3/t/o027/jw05324-o027_t022_nircam_clear-f277w_i2d.fits",
          "october":"jwst/public/jw12577/L3/t/o001/jw12577-o001_t001_nircam_clear-f277w_i2d.fits"}


def main(cache,out):
    cache.mkdir(exist_ok=True,parents=True);out.mkdir(exist_ok=True,parents=True)
    frozen_path=Path(__file__).resolve().parent/'transients-goods5324-controls-receipts.json'
    expected={}
    if frozen_path.exists():
        frozen=json.loads(frozen_path.read_text())
        expected={(r['url'],r['first_byte'],r['last_byte']):r['sha256'] for r in frozen['ranges']}
    receipt={"run_utc":datetime.now(timezone.utc).isoformat(),"stage":"post-inspection diagnostic",
             "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             "selection":"October F277W is_extended false, catalog SNR>28, nn_dist>5; 0.8<separation<=60arcsec; nearest8",
             "ranges":[],"products":{}}
    def fetch(key,start=None,end=None):
        url=BASE+key;identity=f"{url}:{start}:{end}"
        path=cache/(hashlib.sha256(identity.encode()).hexdigest()+".bin")
        if path.exists():raw=path.read_bytes()
        else:
            req=Request(url,headers={"Range":f"bytes={start}-{end}"} if start is not None else {})
            with urlopen(req,timeout=45) as res:
                if start is not None and (res.status!=206 or not res.headers.get("Content-Range","").startswith(f"bytes {start}-{end}/")):
                    raise RuntimeError("Invalid range response")
                raw=res.read()
            path.write_bytes(raw)
        if start is not None and len(raw)!=end-start+1:raise RuntimeError("Invalid byte count")
        if expected and hashlib.sha256(raw).hexdigest()!=expected.get((url,start,end)):
            raise RuntimeError('Input differs from frozen receipt; use a separate campaign')
        receipt["ranges"].append({"url":url,"first_byte":start,"last_byte":end,"bytes":len(raw),"sha256":hashlib.sha256(raw).hexdigest()})
        return raw
    cat=pd.read_csv(BytesIO(fetch(CAT)),sep=r"\s+",comment="#",na_values=["nan","NaN","--",'""'])
    ra0,dec0=53.188672,-27.812554
    dist=np.hypot((cat["sky_centroid.ra"]-ra0)*np.cos(np.deg2rad(dec0)),cat["sky_centroid.dec"]-dec0)*3600
    select=(dist>.8)&(dist<=60)&(~cat.is_extended)&(cat.aper_total_flux/cat.aper_total_flux_err>28)&(cat.nn_dist>5)
    refs=cat.loc[select].copy();refs["distance_arcsec"]=dist[select];refs=refs.sort_values("distance_arcsec").head(8)
    points=[{"label":int(r.label),"ra":float(r["sky_centroid.ra"]),"dec":float(r["sky_centroid.dec"]),"distance_arcsec":float(r.distance_arcsec)} for _,r in refs.iterrows()]
    receipt["selected_controls"]=points
    (out/"transients-goods5324-controls-selection.json").write_text(json.dumps(receipt,indent=2)+"\n")
    gy,gx=np.mgrid[-44:45,-44:45]*.03;rr=np.hypot(gx,gy)
    ap=rr<=.15;ann=(rr>=.27)&(rr<=.42);pix_njy=np.deg2rad(.03/3600)**2*1e15
    results={}
    for epoch,key in PRODUCTS.items():
        stream=BytesIO(fetch(key,0,57599));ph=fits.Header.fromfile(stream);h=fits.Header.fromfile(stream);off=stream.tell()
        assert (h['EXTNAME'],h['BITPIX'],h['BUNIT'])==('SCI',-32,'MJy/sr')
        w=WCS(h,naxis=2);nx,ny=h['NAXIS1'],h['NAXIS2'];radius=int(np.ceil(1.8/(abs(h['CDELT1'])*3600)))
        receipt['products'][epoch]={'utc':ph.get('DATE-BEG'),'shape':[ny,nx],'sci_data_offset':off,
                                   'pipeline_version':ph.get('CAL_VER'),'crds_context':ph.get('CRDS_CTX')}
        for point in points:
            x,y=w.all_world2pix([[point['ra'],point['dec']]],0)[0];xc,yc=int(np.rint(x)),int(np.rint(y))
            record={**point,'epoch':epoch,'status':'no_coverage'}
            if radius<=xc<nx-radius and radius<=yc<ny-radius:
                y0,y1=yc-radius,yc+radius+1
                data=np.frombuffer(fetch(key,off+y0*nx*4,off+y1*nx*4-1),'>f4').reshape(y1-y0,nx)
                xx,yy=w.all_world2pix(point['ra']+gx/(3600*np.cos(np.deg2rad(point['dec']))),point['dec']+gy/3600,0)
                image=map_coordinates(data,[yy-y0,xx],order=1,mode='constant',cval=np.nan)
                if np.isfinite(image[ap]).all() and np.isfinite(image[ann]).mean()>=.9:
                    flux=float(np.sum(image[ap]-np.nanmedian(image[ann]))*pix_njy)
                    record.update(status='measured',aperture_njy=flux)
            results[epoch,point['label']]=record
    rows=[]
    for point in points:
        a,b=results['february',point['label']],results['october',point['label']]
        row={**point,'february_status':a['status'],'october_status':b['status'],
             'february_njy':a.get('aperture_njy'),'october_njy':b.get('aperture_njy')}
        if a.get('aperture_njy',0)>0 and b.get('aperture_njy',0)>0:row['february_over_october']=a['aperture_njy']/b['aperture_njy']
        rows.append(row)
    pd.DataFrame(rows).to_csv(out/'transients-goods5324-controls.csv',index=False)
    ratios=np.array([r['february_over_october'] for r in rows if 'february_over_october'in r])
    receipt['summary']={'selected':len(points),'finite_positive_pairs':len(ratios),
                        'median_february_over_october':float(np.median(ratios)) if len(ratios) else None,
                        'ratio_MAD_scaled':float(1.4826*np.median(abs(ratios-np.median(ratios)))) if len(ratios) else None,
                        'limits':'Post-inspection small-sample comparison. No variability vetting, native-DQ examination, matched-PSF convolution, or full local registration solution.'}
    (out/'transients-goods5324-controls-receipts.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt['summary'],indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--cache',required=True,type=Path);parser.add_argument('--out',required=True,type=Path)
    a=parser.parse_args();main(a.cache,a.out)
