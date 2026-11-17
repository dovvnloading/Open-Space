"""Independent F444W dither subsets at the unresolved Warhol position."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
from urllib.request import Request, urlopen

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from scipy import ndimage as nd
from scipy.optimize import minimize
import macs0416 as obs

HERE=Path(__file__).resolve().parent
# One matched detector in each visit; the central arc is on module B.
VISITS=(("jw01176213001","03101"),("jw06882054001","02101"))
X,Y=345,378
R=70

def native(item):
    epoch,num,sky=item
    visit,group=VISITS[epoch]
    name=f"{visit}_{group}_{num:05d}_nrcblong_i2d.fits"
    url=f"https://stpubdata.s3.amazonaws.com/jwst/public/{visit[:7]}/{visit}/{name}"
    rec={"epoch":epoch,"dither":num,"url":url,"ranges":[]}
    def read(a,n):
        import time
        for attempt in range(4):
            try:
                with urlopen(Request(url,headers={"Range":f"bytes={a}-{a+n-1}"}),timeout=20) as h:
                    if h.status!=206 or not h.headers["Content-Range"].startswith(f"bytes {a}-{a+n-1}/"):
                        raise IOError("Invalid response range")
                    v=h.read()
                    etag=h.headers["ETag"]
                break
            except OSError as exc:
                print(name, a, n, type(exc).__name__, str(exc), flush=True)
                if attempt == 3: raise
                time.sleep(attempt + 1)
        if len(v)!=n:raise IOError("Short FITS window")
        rec["ranges"].append(dict(start=a,bytes=n,sha256=sha256(v).hexdigest(),etag=etag))
        return v
    buf=BytesIO(read(0,57600));primary=fits.Header.fromfile(buf);hdr=fits.Header.fromfile(buf)
    assert primary["FILTER"]=="F444W" and hdr['EXTNAME']=='SCI' and hdr['BITPIX']==-32 and hdr["BUNIT"]=="MJy/sr"
    rec["mjd"]=hdr["MJD-BEG"] if "MJD-BEG" in hdr else primary["MJD-BEG"]
    w=WCS(hdr,naxis=2);xx,yy=w.world_to_pixel_values(*sky)
    width=hdr['NAXIS1']; y0=int(np.floor(yy.min()))-2;y1=int(np.ceil(yy.max()))+3
    assert min(xx.min(),y0)>0 and xx.max()<width-1 and y1<hdr['NAXIS2']
    b=read(buf.tell()+4*y0*width,4*(y1-y0)*width)
    arr=np.frombuffer(b,dtype='>f4').reshape(y1-y0,width)
    im=nd.map_coordinates(arr,[yy-y0,xx],order=1,cval=np.nan)
    return epoch,num,im,rec

def register(im, reference):
    y,x=np.indices(im.shape);rr=np.hypot(x-R,y-R)
    g=(rr>11)&(rr<57)&np.isfinite(reference)&np.isfinite(im)
    # Fit robust translation to extended structure outside the flag center.
    ref=nd.gaussian_filter(reference,1.0); img=nd.gaussian_filter(np.where(np.isfinite(im),im,np.nanmedian(im)),1.0)
    def loss(dxy):
        shifted=nd.shift(img,(dxy[1],dxy[0]),order=1,mode='nearest',prefilter=False)
        delta=shifted[g]-ref[g]
        delta-=np.median(delta)
        return np.mean(np.minimum(delta*delta,.25))
    dx,dy=minimize(loss,[0,0],method='Powell',bounds=[(-8.,8.)]*2).x
    return nd.shift(im,(dy,dx),order=1,mode='nearest',prefilter=False),[float(dx),float(dy)]

def main():
    import sys
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("cache",type=Path);p.add_argument("--base",type=Path,default=obs.LEGACY)
    args=p.parse_args()
    m=obs.module("analyze_elgordo",args.base);m.CENTER=obs.CENTER;m.SIDE=601
    grid=m.getgrid();paths=[args.cache/f"{ep}_nircam_clear-f444w.fits" for ep in obs.EPOCHS]
    d=m.analyze("f444w",paths,grid)
    gy,gx=np.mgrid[Y-R:Y+R+1,X-R:X+R+1];sky=grid.pixel_to_world_values(gx,gy)
    sl=np.s_[Y-R:Y+R+1,X-R:X+R+1]
    data=[[],[]];receipts=[]
    with ThreadPoolExecutor(4) as pool:
        for epoch,num,img,rc in pool.map(native,[(e,n,sky) for e in (0,1) for n in range(1,5)]):
            registered,shift=register(img,d["im"+str(epoch)][sl])
            rc["translation_xy_pixels"]=shift
            data[epoch].append(registered)
            receipts.append(rc)
    rows=[]
    y,x=np.indices(gx.shape); ann=(np.hypot(x-R,y-R)>16)&(np.hypot(x-R,y-R)<45)
    good=d["good"][sl]&ann
    for tag,inds in [("1+3",[0,2]),("2+4",[1,3])]+[(str(i+1),[i]) for i in range(4)]:
        images=[np.nanmean([data[e][i] for i in inds],axis=0) for e in (0,1)]
        images=[np.where(np.isfinite(im),im,np.nanmedian(im)) for im in images]
        diff=m.convolve(images[1],d['psfs'][0])-m.convolve(images[0],d['psfs'][1])
        amp=m.response(diff,d['effective'])
        vals=amp[good]; amp=(amp[R,R]-np.median(vals))*d["to_njy"]/d["resp"]
        noise=m.mad(vals)*d["to_njy"]/d["resp"]
        rows.append(dict(dithers_each_epoch=tag,delta_njy=amp,local_annulus_sigma_njy=noise,signed_score=amp/noise,
                         annulus_valid_pixels=int(good.sum())))
    with (HERE/"warhol_split.csv").open("w",newline="") as h:
        writer=csv.DictWriter(h,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    (HERE/"dither_ranges.json").write_text(json.dumps(receipts,indent=1)+"\n")

if __name__=="__main__":main()
