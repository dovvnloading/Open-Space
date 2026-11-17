"""Forced apertures and shifted-position controls for the C/2025 V1 footprint list."""
import base64
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import sys
from time import sleep
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from astropy import units as u
from astropy.wcs import WCS

BASE=Path(sys.argv[1])
OUT=Path(sys.argv[2])
OUT.mkdir(parents=True,exist_ok=True)
METADATA=BASE/"archive-results.json"
HEADERS={"User-Agent":"Open-Space-archive/1.0 (+https://github.com/dovvnloading/Open-Space)"}


def fetch(url, limit=8_000_000):
    with urlopen(Request(url,headers=HEADERS),timeout=50) as response:
        value=response.read(limit+1)
        if len(value)>limit:
            raise ValueError("Response exceeds requested cutout size")
        return value,{"url":response.url,"bytes":len(value),"sha256":hashlib.sha256(value).hexdigest(),
                      "http_date":response.headers.get("Date")}


def at_times(times):
    answers={}
    receipts=[]
    for i in range(0,len(times),25):
        selection=times[i:i+25]
        params={"format":"json","COMMAND":"'DES=C/2025 V1;'","EPHEM_TYPE":"'OBSERVER'","CENTER":"'I41'",
                "TLIST":"'"+",".join(f"{t:.9f}" for t in selection)+"'",
                "QUANTITIES":"'1,36'","ANG_FORMAT":"'DEG'","CSV_FORMAT":"'YES'",
                "CAL_FORMAT":"'JD'","EXTRA_PREC":"'YES'"}
        url="https://ssd.jpl.nasa.gov/api/horizons.api?"+urlencode(params)
        raw,rec=fetch(url)
        receipts.append(rec)
        s=json.loads(raw)["result"]
        head,data=s.split("$$SOE")
        keys=[k.strip() for k in next(csv.reader([next(l for l in reversed(head.splitlines()) if "R.A." in l)]))]
        ra=next(j for j,k in enumerate(keys) if "R.A." in k)
        dec=next(j for j,k in enumerate(keys) if "DEC" in k.upper() and j>ra)
        for line in csv.reader(data.split("$$EOE")[0].strip().splitlines()):
            jd=float(line[0])
            iref=min(selection,key=lambda t:abs(jd-t))
            if abs(iref-jd)>1e-5:
                raise ValueError("Ephemeris time mismatch")
            answers[iref]=dict(ra=float(line[ra]),dec=float(line[dec]),
                               uncertainties={k:v.strip() for k,v in zip(keys[dec+1:],line[dec+1:]) if k})
        if not all(t in answers for t in selection):
            raise ValueError("Incomplete precise ephemeris")
    return answers,receipts


def image_url(row,suffix):
    tag=str(row["filefracday"])
    filename=f"ztf_{tag}_{int(row['field']):06d}_{row['filtercode']}_c{int(row['ccdid']):02d}_o_q{row['qid']}_{suffix}"
    base="https://irsa.ipac.caltech.edu/ibe/data/ztf/products/sci/"+"/".join([tag[:4],tag[4:8],tag[8:],filename])
    return base+"?"+urlencode({"center":f"{row['pred_ra_deg']:.8f},{row['pred_dec_deg']:.8f}",
                              "size":"160arcsec","gzip":"false"})


def read_stamp(row,suffix):
    uri=image_url(row,suffix)
    attempts=[]
    for attempt,address in enumerate([uri,uri.split("?")[0],uri]):
        if attempt:
            sleep(1.5**attempt)
        try:
            full="?" not in address
            raw,receipt=fetch(address,48_000_000 if full else 8_000_000)
            with fits.open(io.BytesIO(raw)) as hdus:
                h=next(p for p in hdus if p.data is not None and p.data.ndim==2)
                hdr=h.header.copy()
                data=h.data
                if full:
                    pos=SkyCoord(row["pred_ra_deg"],row["pred_dec_deg"],unit="deg")
                    cut=Cutout2D(data,pos,160*u.arcsec,wcs=WCS(hdr).celestial)
                    data=cut.data
                    hdr.update(cut.wcs.to_header(relax=True))
                receipt["earlier_attempts"]=attempts
                return data.astype(float),hdr,receipt
        except Exception as exc:
            attempts.append({"url":address,"error":str(exc)})
    raise OSError(json.dumps(attempts))


def aperture(data,x,y):
    xc,yc=round(x),round(y)
    if xc<17 or yc<17 or xc+17>=data.shape[1] or yc+17>=data.shape[0]:
        return np.nan
    patch=data[yc-17:yc+18,xc-17:xc+18]
    yg,xg=np.ogrid[-17:18,-17:18]
    d2=(xg+xc-x)**2+(yg+yc-y)**2
    ap=patch[d2<=36]
    ring=patch[(d2>=121)&(d2<=225)]
    if not np.isfinite(ap).all() or np.isfinite(ring).sum()<200:
        return np.nan
    return float(ap.sum()-len(ap)*np.nanmedian(ring))


def measure(stamp,point,zp):
    data,hdr,receipt=stamp
    x,y=map(float,WCS(hdr).celestial.all_world2pix(point["ra"],point["dec"],0))
    flux=aperture(data,x,y)
    controls=np.array([aperture(data,x+dx,y+dy) for dx in range(-52,53,14) for dy in range(-52,53,14)
                       if np.hypot(dx,dy)>16])
    controls=controls[np.isfinite(controls)]
    if not np.isfinite(flux) or len(controls)<12:
        return {**receipt,"error":"target or control apertures lack valid pixels"}
    med=float(np.median(controls))
    mad=float(1.4826*np.median(np.abs(controls-med)))
    return {**receipt,"x0":x,"y0":y,"aperture_dn":flux,"control_count":len(controls),
            "control_median_dn":med,"control_mad_sigma_dn":mad,
            "controls_at_least_target":int(np.sum(controls>=flux)),
            "snr_empirical":(flux-med)/mad if mad else None,
            "uncorrected_aperture_mag":float(zp)-2.5*math.log10(flux) if zp is not None and flux>0 else None}


def obtain(row):
    try:
        sci=read_stamp(row,"sciimg.fits")
        hdr=sci[1]
        jd=float(hdr["OBSJD"])+float(hdr["EXPTIME"])/172800
        return row,sci,jd
    except Exception as exc:
        return row,None,{"error":str(exc),"url":image_url(row,"sciimg.fits")}


def main():
    products=json.loads(METADATA.read_text())["ztf"]
    if "ztf_nominal_footprints" not in products:
        raise ValueError(f"Footprint search incomplete: {products}")
    images=products["ztf_nominal_footprints"]
    with ThreadPoolExecutor(max_workers=3) as pool:
        results=list(pool.map(obtain,images))
    times=sorted(set(jd for _,s,jd in results if s))
    points,receipts=at_times(times)
    out=[]
    for row,sci,jd in results:
        entry={"filefracday":row["filefracday"],"field":row["field"],"filter":row["filtercode"],
               "ccdid":row["ccdid"],"qid":row["qid"]}
        out.append(entry)
        if sci is None:
            entry.update(jd)
            continue
        point=points[jd]
        entry.update(mid_jd_utc=jd,position=point,seeing=sci[1].get("SEEING"),infobits=sci[1].get("INFOBITS"),
                     magzp=sci[1].get("MAGZP"))
        try:
            mask,mask_header,rec=read_stamp(row,"mskimg.fits")
            if mask.shape!=sci[0].shape:
                raise ValueError("mask shape mismatch")
            ny,nx=mask.shape
            coords=np.array([[0.,0.],[nx-1.,0.],[0.,ny-1.],[nx-1.,ny-1.],[(nx-1.)/2,(ny-1.)/2]])
            world=WCS(sci[1]).celestial.all_pix2world(coords,0)
            back=WCS(mask_header).celestial.all_world2pix(world,0)
            if not np.isfinite(back).all() or np.max(np.abs(back-coords))>1e-3:
                raise ValueError("mask WCS mismatch")
            entry["mask"]=rec
            sci[0][mask!=0]=np.nan
        except Exception as exc:
            entry["mask_error"]=str(exc)
            entry["science"]={"error":"required science mask unavailable"}
            entry["difference"]={"error":"no validated difference-image mask"}
            continue
        entry["science"]=measure(sci,point,entry["magzp"])
        # Difference products need their own validated quality mask and WCS.
        entry["difference"]={"error":"no validated difference-image mask"}
        print("MEASUREMENT",entry["filefracday"],entry["science"].get("snr_empirical"),entry["difference"].get("snr_empirical"),flush=True)
    result={"footprint_input_sha256":hashlib.sha256(METADATA.read_bytes()).hexdigest(),"horizons_receipts":receipts,"exposures":out}
    dest=OUT/"ztf_forced.json"
    dest.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("RESULT",dest.name,base64.b64encode(dest.read_bytes()).decode())

if __name__=="__main__":
    main()
