"""Search each pre-discovery night in ZTF positive difference images along the comet's motion.

Use a fixed local-background-subtracted aperture in every exposure. Scan a 1 arcsec
grid across the Horizons 3-sigma ellipse and compare each night's maximum to the
same maximum measured on nearby, parallel tracks. Report flux in aperture DN at
a common zero point when image headers define that calibration.
"""
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import argparse
import csv
import hashlib
import io
import json
import math
import time

from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import convolve, map_coordinates
import numpy as np

AGENT = "Open-Space-ZTF-stacks/1.0 (+https://github.com/dovvnloading/Open-Space)"
GRID = 1.0
RADIUS = 3.5  # pixels; ZTF sampling ~1 arcsec/pixel


def request(url, timeout=100):
    with urlopen(Request(url,headers={"User-Agent":AGENT}),timeout=timeout) as response:
        raw = response.read()
        return raw, {"url":url,"status":response.status,"date":response.headers.get("Date"),
                     "bytes":len(raw),"sha256":hashlib.sha256(raw).hexdigest()}


def image_url(r,suffix):
    d=r["filefracday"]
    f=f"ztf_{d}_{int(r['field']):06d}_{r['filtercode']}_c{int(r['ccdid']):02d}_o_q{r['qid']}_{suffix}"
    return "https://irsa.ipac.caltech.edu/ibe/data/ztf/products/sci/"+"/".join([d[:4],d[4:8],d[8:],f])


def get_array(raw):
    with fits.open(io.BytesIO(raw),memmap=False) as hdul:
        hdu=next(h for h in hdul if h.data is not None and h.data.ndim==2)
        header=hdul[0].header.copy()
        header.update(hdu.header)
        return np.asarray(hdu.data,dtype=float),header


def prepare(row, receipts):
    approx=(float(row["pred_ra_deg"]),float(row["pred_dec_deg"]))
    params=urlencode({"center":"%.7f,%.7f"%approx,"size":"800arcsec","gzip":"false"})
    uri=image_url(row,"scimrefdiffimg.fits.fz")
    for idx,url in enumerate([uri+"?"+params,uri]):
        try:
            raw,req=request(url)
            data,hdr=get_array(raw)
            receipts.append(req)
            if min(data.shape)<100:
                raise ValueError("Unexpected postage stamp size")
            break
        except Exception as exc:
            receipts.append({"url":url,"error":str(exc)})
            if idx==1:
                raise
    wcs=WCS(hdr).celestial
    x,y=map(float,wcs.all_world2pix(*approx,0))
    x0,y0=max(0,round(x)-405),max(0,round(y)-405)
    data=data[y0:round(y)+406,x0:round(x)+406]
    if min(data.shape)<100:
        raise ValueError("Target region falls outside usable difference-image WCS")
    # Difference products generally carry ZPDIFF. Preserve the actual keyword.
    zpkey=next((k for k in ("ZPDIFF","MAGZPDIFF","MAGZP") if hdr.get(k) is not None),None)
    if hdr.get("OBSMJD") is not None and hdr.get("EXPTIME") is not None:
        mid=float(hdr["OBSMJD"])+2400000.5+float(hdr["EXPTIME"])/172800
        clock="FITS OBSMJD + EXPTIME/2"
    elif hdr.get("OBSJD") is not None and hdr.get("EXPTIME") is not None:
        mid=float(hdr["OBSJD"])+float(hdr["EXPTIME"])/172800
        clock="FITS OBSJD + EXPTIME/2"
    else:
        sciuri=image_url(row,"sciimg.fits")+"?"+urlencode({"center":"%.7f,%.7f"%approx,"size":"20arcsec","gzip":"false"})
        b,rec=request(sciuri)
        receipts.append(rec)
        _,sc=get_array(b)
        mjd = sc.get("OBSMJD")
        jd = float(mjd)+2400000.5 if mjd is not None else float(sc["OBSJD"])
        mid=jd+float(sc["EXPTIME"])/172800
        clock="science FITS shutter epoch + EXPTIME/2"
    yy,xx=np.mgrid[-10:11,-10:11]
    rr=xx*xx+yy*yy
    ap=rr<=RADIUS**2
    bg=(rr>=7**2)&(rr<=10**2)
    kernel=ap.astype(float)-bg*(ap.sum()/bg.sum())
    valid=np.isfinite(data)
    # A gap or bad pixel within the aperture invalidates that trial.
    bad=convolve((~valid).astype(float),(ap|bg).astype(float),mode="constant",cval=1)>0
    values=convolve(np.where(valid,data,0),kernel,mode="constant",cval=0)
    values[bad]=np.nan
    return {"row":row,"map":values,"x0":x0,"y0":y0,"wcs":wcs,"jd":mid,
            "clock":clock,"zp":float(hdr[zpkey]) if zpkey else None,"zpkey":zpkey,
            "seeing":hdr.get("SEEING"),"infobits":hdr.get("INFOBITS"),"shape":list(data.shape),
            "file":row["filefracday"]}


def ephemeris(frames,receipts,out):
    # Horizons TLIST matches supplied UTC Julian dates, rounded only here.
    epochs=sorted(set(round(f["jd"],8) for f in frames))
    args={"format":"json","COMMAND":"'DES=C/2025 V1;'","CENTER":"'I41'","EPHEM_TYPE":"'OBSERVER'",
          "TLIST":"'"+",".join(f"{x:.8f}" for x in epochs)+"'","QUANTITIES":"'1,36,37'",
          "CAL_FORMAT":"'JD'","ANG_FORMAT":"'DEG'","EXTRA_PREC":"'YES'","CSV_FORMAT":"'YES'"}
    url="https://ssd.jpl.nasa.gov/api/horizons.api?"+urlencode(args)
    raw,rec=request(url)
    receipts.append(rec)
    (out/"stack-horizons.json").write_bytes(raw)
    result=json.loads(raw)["result"]
    if "$$SOE" not in result:
        raise ValueError(result[:1000])
    before,body=result.split("$$SOE",1)
    header=next(csv.reader([next(x for x in reversed(before.splitlines()) if "R.A." in x)]))
    def idx(p):return next(i for i,x in enumerate(header) if p.lower() in x.lower())
    fields=[idx(x) for x in ("R.A.","DEC_","SMAA","SMIA","Theta")]
    vals={}
    for row in csv.reader(body.split("$$EOE")[0].strip().splitlines()):
        vals[float(row[0])]=[float(row[k]) for k in fields]
    if len(vals)!=len(epochs):
        raise ValueError("Missing Horizons epochs")
    for f in frames:
        matched=min(vals,key=lambda jd:abs(jd-f["jd"]))
        if abs(matched-f["jd"])>2e-5:
            raise ValueError("Unmatched shutter time")
        f["ra"],f["dec"],f["a"],f["b"],f["theta"]=vals[matched]


def sample(f, east, north, sky=None):
    ra,dec=sky if sky else (f["ra"],f["dec"])
    ra=(ra+east/(3600*np.cos(np.radians(dec))))%360
    dec=dec+north/3600
    x,y=f["wcs"].all_world2pix(ra,dec,0)
    return map_coordinates(f["map"],[y-f["y0"],x-f["x0"]],order=1,mode="constant",cval=np.nan)


def scan(frames):
    frames=sorted(frames,key=lambda f:f["jd"])
    mid=frames[len(frames)//2]
    # All frames in a night span minutes; add one arcsec for WCS/systematics
    # and include an entire pixel of tolerance at the ellipse boundary.
    a=max(f["a"] for f in frames)+1.5
    b=max(f["b"] for f in frames)+1.5
    angle=np.radians(mid["theta"])
    extent=math.ceil(a)
    n,e=np.mgrid[-extent:extent+1:GRID,-extent:extent+1:GRID]
    mask=((np.cos(angle)*e+np.sin(angle)*n)/a)**2+((-np.sin(angle)*e+np.cos(angle)*n)/b)**2<=1
    E,N=e[mask],n[mask]
    # Nonoverlapping sky boxes for the comparison tracks, also keeping each
    # box clear of all other target positions used in this night.
    step=math.ceil(max(2*a+12,32))
    controls=[(x*step,y*step) for y in range(-6,7) for x in range(-6,7)
              if (x or y) and max(abs(x*step),abs(y*step))+a<365]
    shifted=np.array([(0,0)]+controls)
    east=shifted[:,0,None]+E
    north=shifted[:,1,None]+N
    samples=np.stack([sample(f,east,north) for f in frames])
    if all(f["zp"] is not None for f in frames):
        scales=np.array([10**(.4*(25-f["zp"])) for f in frames])
        unit="DN at zero point 25"
    else:
        scales=np.ones(len(frames))
        unit="native per-frame DN"
    samples*=scales[:,None,None]
    center=int(np.argmin(E**2+N**2))
    median=np.nanmedian(samples[:,1:,center],axis=1)
    sig=1.4826*np.nanmedian(abs(samples[:,1:,center]-median[:,None]),axis=1)
    if not np.all(sig>0):
        raise ValueError("Indeterminate empirical frame noise")
    w=1/sig**2
    norm=w.sum()
    stack=((samples-median[:,None,None])*w[:,None,None]).sum(axis=0)/norm
    # Control trial count equals target trial count; discard any control with
    # missing pixels rather than comparing maxima of different-sized regions.
    if not np.all(np.isfinite(stack[0])):
        raise ValueError("Masked target grid; inspection required")
    good=np.isfinite(stack).all(axis=1)
    off=stack[1:][good[1:]]
    if len(off)<20:
        raise ValueError("Fewer than 20 fully sampled control ellipses")
    peak=np.argmax(stack[0])
    offmax=off.max(axis=1)
    empirical_p=(int(np.sum(offmax>=stack[0,peak]))+1)/(len(off)+1)
    noises=1.4826*np.median(abs(off[:,center]-np.median(off[:,center])))
    p95=float(np.quantile(offmax,.95))
    # Compare a point fixed on the sky at the first target position against
    # the co-moving trial. Difference frames already remove static sky.
    ecs=np.array([[E[peak]]]);ncs=np.array([[N[peak]]])
    sky=(frames[0]["ra"],frames[0]["dec"])
    stationary=[float(sample(f,ecs,ncs,sky=sky)[0,0])*s for f,s in zip(frames,scales)]
    sn=(samples[:,0,peak]-median)/sig
    return {"frames":[f["file"] for f in frames],"track_extent_arcsec":[float(a),float(b)],
            "grid_step_arcsec":GRID,"grid_points":len(E),"valid_offset_tracks":len(off),
            "best_east_north_arcsec":[float(E[peak]),float(N[peak])],"flux_units":unit,
            "aperture_radius_pix":RADIUS,"nominal_flux":float(stack[0,center]),
            "best_flux":float(stack[0,peak]),"stack_control_sigma":float(noises),
            "peak_flux_over_sigma":float(stack[0,peak]/noises),"offtrack_p_for_searched_peak":empirical_p,
            "control_peak_95th_percentile":p95,"aperture_5sigma_flux":float(5*noises),
            "aperture_5sigma_mag_zp25":float(25-2.5*np.log10(5*noises)) if unit.endswith("25") else None,
            "peak_frame_snr":sn.tolist(),"stationary_at_first_frame_snr":((stationary-median)/sig).tolist(),
            "frame_flux_sigmas":sig.tolist(),"nominal_frame_snr":((samples[:,0,center]-median)/sig).tolist()}


def main(folder):
    dataset=json.loads((folder/"archive-results.json").read_text())
    rows=dataset["ztf"]["ztf_nominal_footprints"]
    # The first catalogued night serves as a blind pipeline positive control.
    rows=[r for r in rows if r["filefracday"][:8]<="20251029"]
    out=folder/"stack"
    out.mkdir(exist_ok=True)
    receipts,errors,frames=[],[],[]
    for r in rows:
        try:
            f=prepare(r,receipts)
            frames.append(f)
            print("IMAGE OK",f["file"],f["zpkey"],f["zp"],f["clock"],flush=True)
        except Exception as exc:
            errors.append({"file":r["filefracday"],"error":str(exc)})
            print("IMAGE ISSUE",errors[-1],flush=True)
    if frames:
        ephemeris(frames,receipts,out)
    groups=defaultdict(list)
    for f in frames:
        groups[f["file"][:8]].append(f)
    result={}
    for date,group in sorted(groups.items()):
        try:
            if len(group)<2:
                raise ValueError("Insufficient images to check motion")
            result[date]=scan(group)
        except Exception as exc:
            result[date]={"error":str(exc)}
        print("STACK",date,json.dumps(result[date]),flush=True)
    summary={"targets":len(rows),"loaded":len(frames),"errors":errors,"nights":result,
             "frames":[{k:v for k,v in f.items() if k not in {"map","wcs","row"}} for f in frames]}
    (out/"stack-results.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    (out/"stack-receipts.json").write_text(json.dumps(receipts,indent=2)+"\n")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("folder",type=Path)
    main(p.parse_args().folder)
