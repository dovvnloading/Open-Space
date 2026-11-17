"""Test stationary contamination of the Oct 25 ZTF aperture at its measured sky position."""
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import time
from urllib.parse import urlencode
from urllib.request import Request,urlopen

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np

OUT=Path(__file__).parent/"image-search"/"oct25"
OUT.mkdir(exist_ok=True)
RA,DEC=167.839372253,-11.366983117
MID=2460974.041956011
receipts=[]
def request(name,base,pars=None):
    url=base+("?"+urlencode(pars) if pars else "")
    rec={"file":name,"url":url}
    receipts.append(rec)
    try:
        with urlopen(Request(url,headers={"User-Agent":"Open-Space-static-sky/1.0 (+https://github.com/dovvnloading/Open-Space)"}),timeout=60) as resp:
            raw=resp.read()
            rec.update(status=resp.status,http_date=resp.headers.get("Date"),bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest())
        if not name.endswith(".fits"):
            (OUT/(name+".gz")).write_bytes(gzip.compress(raw,mtime=0))
        return raw
    except Exception as e:
        rec["error"]=str(e)
        return None

def offset(ra,dec):
    east=((float(ra)-RA+180)%360-180)*math.cos(math.radians(DEC))*3600
    north=(float(dec)-DEC)*3600
    return [east,north,math.hypot(east,north)]

ps1=request("panstarrs.csv","https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/mean.csv",
            {"ra":RA,"dec":DEC,"radius":.025,"columns":"[objID,raMean,decMean,nDetections,rMeanPSFMag,rMeanKronMag,iMeanPSFMag,gMeanPSFMag]"})
psources=[]
if ps1:
    for r in csv.DictReader(io.StringIO(ps1.decode())):
        if not r.get("raMean"):
            continue
        p=offset(r["raMean"],r["decMean"])
        psources.append({**r,"offset_ENR_arcsec":p})
psources.sort(key=lambda r:r["offset_ENR_arcsec"][2])
gaia=request("gaia.tsv","https://vizier.cds.unistra.fr/viz-bin/asu-tsv",{"-source":"I/355/gaiadr3","-c":f"{RA} {DEC}","-c.rs":"90","-out":"Source,RA_ICRS,DE_ICRS,Gmag,pmRA,pmDE","-out.max":"5000"})
gsources=[]
if gaia:
    ls=[l for l in gaia.decode().splitlines() if l.strip() and not l.startswith("#")]
    for r in csv.DictReader(ls,delimiter="\t"):
        try:
            ra,dec=float(r["RA_ICRS"]),float(r["DE_ICRS"])
            years=(MID-2457388.5)/365.25
            if r["pmRA"].strip() and r["pmDE"].strip():
                ra+=float(r["pmRA"])*years/(3.6e6*math.cos(math.radians(dec)))
                dec+=float(r["pmDE"])*years/3.6e6
            p=offset(ra,dec)
            gsources.append({**r,"offset_ENR_arcsec":p})
        except (ValueError,KeyError):
            continue
gsources.sort(key=lambda r:r["offset_ENR_arcsec"][2])
print("STATIC CATALOGS",json.dumps([psources[:3],gsources[:3]]),flush=True)

# Apply the same six-pixel aperture at all three sky positions and at the
# first aperture in independent epochs. This tests persistence at fixed sky.
rows=json.loads((OUT.parent/"archive-results.json").read_text())["ztf"]["ztf_nominal_footprints"]
rows=sorted([r for r in rows if r["filefracday"].startswith("20251025")],key=lambda r:r["filefracday"])
jd=[float(r["obsjd"])+15/86400 for r in rows]  # These Palomar science visits are 30 s.
raw=request("horizons.json","https://ssd.jpl.nasa.gov/api/horizons.api",{"format":"json","COMMAND":"'DES=C/2025 V1;'","CENTER":"'I41'","EPHEM_TYPE":"'OBSERVER'","TLIST":"'"+",".join(f"{x:.8f}" for x in jd)+"'","QUANTITIES":"'1'","CSV_FORMAT":"'YES'","CAL_FORMAT":"'JD'","ANG_FORMAT":"'DEG'","EXTRA_PREC":"'YES'"})
coords=[(RA,DEC)]
if raw:
    try:
        before,body=json.loads(raw)["result"].split("$$SOE",1)
        head=next(csv.reader([next(l for l in reversed(before.splitlines()) if "R.A." in l)]))
        ri=next(i for i,k in enumerate(head) if "R.A." in k)
        di=next(i for i,k in enumerate(head) if "DEC" in k)
        coords=[(float(r[ri]),float(r[di])) for r in csv.reader(body.split("$$EOE")[0].strip().splitlines())]
        if math.dist(offset(*coords[0])[:2],[0,0])>.5: raise ValueError("Unexpected original coordinate")
    except Exception as e:
        print("EPOCH ISSUE",str(e),flush=True)
        coords=[(RA,DEC)]
def measure(raw):
    with fits.open(io.BytesIO(raw),memmap=False) as h:
        hdu=next(hd for hd in h if hd.data is not None and hd.data.ndim==2)
        hdr=h[0].header.copy()
        hdr.update(hdu.header)
        data=np.asarray(hdu.data,dtype=float)
        w=WCS(hdr).celestial
        def ap(x,y):
            xi,yi=round(x),round(y)
            if xi-17<0 or yi-17<0 or xi+18>data.shape[1] or yi+18>data.shape[0]:return float("nan")
            a=data[yi-17:yi+18,xi-17:xi+18]
            iy,ix=np.mgrid[yi-17:yi+18,xi-17:xi+18]
            ds=(ix-x)**2+(iy-y)**2
            core=a[ds<=36]
            sky=a[(ds>=121)&(ds<=225)]
            if not np.all(np.isfinite(core)) or not np.isfinite(sky).any():return float("nan")
            return float(np.sum(core-np.nanmedian(sky)))
        phot=[]
        for ra,dec in coords:
            x,y=map(float,w.all_world2pix(ra,dec,0))
            flux=ap(x,y)
            cs=np.array([ap(x+dx*20,y+dy*20) for dx in range(-3,4) for dy in range(-3,4) if dx or dy])
            cs=cs[np.isfinite(cs)]
            if len(cs)<12 or not np.isfinite(flux):
                phot.append({"error":"Insufficient control or target pixels"})
                continue
            md=np.median(cs)
            sigma=1.4826*np.median(abs(cs-md))
            phot.append({"flux_DN":flux,"background_control_median":float(md),"sigma_DN":float(sigma),
                         "control_n":len(cs),"SNR":float((flux-md)/sigma),"controls_geq":int(sum(cs>=flux))})
        return {"coords":phot,"magzp":hdr.get("MAGZP",hdr.get("ZPDIFF")),"seeing":hdr.get("SEEING"),"shape":list(data.shape)}
params={"center":f"{RA},{DEC}","size":"280arcsec","gzip":"false"}
fields={"reference":"https://irsa.ipac.caltech.edu/ibe/data/ztf/products/ref/001/field001360/zr/ccd13/q4/ztf_001360_zr_c13_q4_refimg.fits"}
for r in rows:
    d=r["filefracday"]
    f=f"ztf_{d}_001360_zr_c13_o_q4"
    pref="https://irsa.ipac.caltech.edu/ibe/data/ztf/products/sci/"+"/".join([d[:4],d[4:8],d[8:],f+"_"])
    for suffix in ("sciimg.fits","scimrefdiffimg.fits.fz"):
        fields[d+"_"+suffix]=pref+suffix
measurements={}
try:
    for label,uri in fields.items():
        for pars in (params,None):
            raw=request(label+".fits",uri,pars)
            if raw:
                try:
                    measurements[label]=measure(raw)
                    break
                except Exception as e:
                    receipts[-1]["FITS_error"]=str(e)
        print("STATIC FRAME",label,measurements.get(label),flush=True)
finally:
    result={"initial_aperture":{"ra":RA,"dec":DEC,"mid_jd":MID,"radius_pixels":6},"coords_jd":jd,"coords":coords,
            "PS1_sources_within_90arcsec":psources,"Gaia_sources_within_90arcsec":gsources,
            "frames":measurements,"receipts":receipts}
    (OUT/"static-sky.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
