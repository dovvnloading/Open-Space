#!/usr/bin/env python3
"""Networked IRSA SIA discovery + S3 FITS byte-range cutouts.
No unsupported ?center=...&size=... query string. Only small science slices.
Pixels inspected algorithmically; detector/exposure-level artefact rejection
and human visual validation are NOT implied by obtaining a PNG.
"""
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
from urllib.parse import urlsplit
import hashlib, json, math, io, csv, sys
import numpy as np
import requests
import s3fs
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"output"
PNGS=OUT/"pixel_previews"
sys.path.insert(0,str(ROOT.parent/"euclid-dr1-wide-highz-vetting"))
import image_review
SIA="https://irsa.ipac.caltech.edu/SIA"
BUCKET="nasa-irsa-euclid-dr1"
BANDS=("VIS","Y","J","H")

def sha(b):return hashlib.sha256(b).hexdigest()
def parse(resp):
    if resp.status_code!=200:raise RuntimeError("SIA HTTP "+str(resp.status_code)+" "+resp.text[:100])
    return image_review.csv_rows(resp)

def key_for(row):
    """Use only original public IRSA Euclid MER S3 objects, not arbitrary URLs."""
    access=(row.get("access_url") or "").strip()
    u0=urlsplit(access)
    if u0.scheme!="https" or u0.hostname not in ("irsa.ipac.caltech.edu","irsacloud.ipac.caltech.edu"):
        raise ValueError("SIA access URL is not trusted IRSA host")
    delim="/ibe/data/euclid/dr1/"
    if delim not in u0.path:raise ValueError("Not an IRSA Euclid DR1 IBE FITS resource")
    suffix=u0.path.split(delim,1)[1]
    if not suffix.startswith("MER/") or not suffix.endswith(".fits"):
        raise ValueError("Not a MER mosaic science FITS object")
    return "dr1/"+suffix

def bounded_stamp(key,ra,dec,band,path):
    # Primary header and science region read via HTTP Range behind s3fs.
    # Never touch hdu.data, which would download the entire mosaic.
    src="s3://"+BUCKET+"/"+key
    coord=SkyCoord(ra*u.deg,dec*u.deg,frame="icrs")
    with fits.open(src,mode="readonly",use_fsspec=True,
                   fsspec_kwargs={"anon":True,
                     "client_kwargs":{"region_name":"us-east-1"}}) as hdul:
        viable=[(idx,h) for idx,h in enumerate(hdul) if
                h.header.get("NAXIS")==2 and
                h.header.get("NAXIS1",0)>0 and h.header.get("NAXIS2",0)>0]
        if not viable:raise ValueError("No 2D image HDU found; not a science-pixel null")
        idx,hdu=viable[0]
        wc=WCS(hdu.header).celestial
        x,y=wc.world_to_pixel(coord)
        if not np.isfinite(x) or not np.isfinite(y):raise ValueError("Invalid target WCS pixel")
        nx=int(hdu.header["NAXIS1"]);ny=int(hdu.header["NAXIS2"])
        if x < 0 or x >= nx or y < 0 or y >= ny:
            raise ValueError("SIA mosaic WCS does not contain requested source")
        scale=np.sqrt(np.abs(np.linalg.det(wc.pixel_scale_matrix)))*3600
        if not np.isfinite(scale) or scale<=0:raise ValueError("Invalid WCS scale")
        half=max(12,min(180,int(np.ceil(12.5/scale))))
        x0=max(0,int(round(x))-half); x1=min(nx,int(round(x))+half+1)
        y0=max(0,int(round(y))-half); y1=min(ny,int(round(y))+half+1)
        if (x1-x0)*(y1-y0)>350000:raise ValueError("Cutout pixel budget exceeded")
        # Astropy section indexes FITS array directly (without .data).
        arr=np.asarray(hdu.section[y0:y1,x0:x1],dtype=np.float32)
        preview,metrics=image_review.stamp(arr)
        data_hash=sha(np.ascontiguousarray(arr,dtype="<f4").tobytes())
        preview_hash=image_review.write_png(path,preview)
        finite=np.isfinite(arr)
        central=arr[max(0,int(round(y))-y0-3):min(arr.shape[0],int(round(y))-y0+4),
                    max(0,int(round(x))-x0-3):min(arr.shape[1],int(round(x))-x0+4)]
        return dict(
          band=band,s3_uri=src,s3_key=key,hdu_index=idx,
          source_pixel_xy=[round(float(x),2),round(float(y),2)],
          selected_bounds=[x0,x1,y0,y1],science_pixel_sha256=data_hash,
          science_pixels=int(finite.sum()),finite_fraction=float(finite.mean()),
          center_7x7_finite_fraction=float(np.isfinite(central).mean()) if central.size else None,
          approximate_pixel_scale_arcsec=round(float(scale),4),
          preview_png=path.relative_to(ROOT).as_posix(),
          preview_sha256=preview_hash,metrics=metrics,
          status="REAL_FITS_SCIENCE_PIXELS_RANGE_READ_NOT_HUMAN_CLEARED",
          image_artefacts_cleared=False,detector_persistence_reviewed=False)

def direct_tile_sources(tileid):
    """S3 public prefix listing bypasses unavailable SIA (HTTP 502).
    Uses MER tileid read from the publisher's HATS catalogue row. Never
    constructs an unverified FITS filename or reads the full mosaic.
    """
    tile=str(int(tileid))
    fs=s3fs.S3FileSystem(anon=True,client_kwargs={"region_name":"us-east-1"},
                         connect_timeout=15,read_timeout=60)
    sources={}
    errors={}
    for dirname in ("VIS","NISP"):
        prefix=f"{BUCKET}/dr1/MER/{tile}/{dirname}"
        try:
            names=fs.ls(prefix,detail=False)
            for name in names:
                key=name.split(BUCKET+"/",1)[1] if BUCKET+"/" in name else name
                tail=key.rsplit("/",1)[-1]
                if "BGSUB-MOSAIC" not in tail or not tail.endswith(".fits"):
                    continue
                if dirname=="VIS":
                    band="VIS" if "-VIS_TILE" in tail else None
                else:
                    band=next((b for b in ("Y","J","H")
                           if f"-NIR-{b}_TILE" in tail),None)
                if band in BANDS and band not in sources:
                    sources[band]=key
        except Exception as e:
            errors[dirname]=type(e).__name__+": "+str(e)[:250]
    return sources,errors

def one(row,rank):
    ra,dec=float(row["ra"]),float(row["dec"])
    result={"rank":rank,"id":row["id"],"ra":ra,"dec":dec,
            "tileid":row.get("tileid"),"query_endpoint":SIA,
            "status":"IMAGE_SOURCES_PENDING",
            "image_artefacts_cleared":False,"bands":[]}
    keys={}
    source={}
    try:
        q=requests.get(SIA,params={"COLLECTION":"euclid_dr1_wide_DpdMerBksMosaic",
              "POS":f"CIRCLE {ra:.8f} {dec:.8f} 0.002",
              "FORMAT":"image/fits","MAXREC":50,"RESPONSEFORMAT":"CSV"},timeout=(12,35))
        result["sia_status"]=q.status_code
        result["sia_sha256"]=sha(q.content)
        records=parse(q)
        result["sia_rows"]=len(records)
        for rec in records:
            band=str(rec.get("energy_bandpassname") or "").strip().upper()
            if band not in BANDS or band in keys:continue
            subtype=str(rec.get("dataproduct_subtype") or "").strip().lower()
            if subtype and subtype!="science":continue
            try:
                key=key_for(rec)
                if band=="VIS" and "/VIS/" not in key:continue
                if band!="VIS" and "/NISP/" not in key:continue
                keys[band]=key
                source[band]="IRSA_SIA_ACCESS_URL"
            except Exception:continue
    except Exception as e:
        result["sia_error_not_nondetection"]=type(e).__name__+": "+str(e)[:310]
    try:
        listing,errors=direct_tile_sources(row["tileid"])
        result["direct_public_S3_listing_errors"]=errors
        result["direct_public_S3_tile"]=str(row["tileid"])
        for band,key in listing.items():
            if band not in keys:
                keys[band]=key;source[band]="PUBLISHER_TILEID_DIRECT_S3_LIST"
    except Exception as e:
        result["direct_public_S3_error_not_nondetection"]=type(e).__name__+": "+str(e)[:330]
    for band in BANDS:
        if band not in keys:continue
        key=keys[band]
        image_path=PNGS/f"rank_{rank:02d}_{band}.png"
        try:
            item=bounded_stamp(key,ra,dec,band,image_path)
            item["discovery_route"]=source[band]
        except Exception as e:
            item={"band":band,"s3_key":key,"discovery_route":source[band],
                  "status":"FITS_PIXEL_ERROR_NOT_A_NEGATIVE_OBSERVATION",
                  "error":type(e).__name__+": "+str(e)[:320],
                  "image_artefacts_cleared":False}
        result["bands"].append(item)
    n=sum(v["status"]=="REAL_FITS_SCIENCE_PIXELS_RANGE_READ_NOT_HUMAN_CLEARED"
          for v in result["bands"])
    result["science_bands_with_pixels"]=n
    result["status"]=("SCIENCE_PIXELS_READ_AUTOMATED_METRICS_ONLY" if n else
                      "NO_USABLE_SCIENCE_PIXELS_NOT_CLEARED")
    return result

def main():
    rows=json.loads((OUT/"top_candidates.json").read_text())
    PNGS.mkdir(parents=True,exist_ok=True)
    result=[]
    with ThreadPoolExecutor(max_workers=2) as pool:
        f={pool.submit(one,r,i+1):i for i,r in enumerate(rows[:8])}
        for future in as_completed(f):
            i=f[future]
            try:x=future.result()
            except Exception as e:x={"rank":i+1,"id":rows[i]["id"],"status":"PIPELINE_ERROR",
                                      "error":str(e)[:300],"image_artefacts_cleared":False}
            result.append((i,x))
            print("S3_PIXEL",i+1,x["status"],flush=True)
    result=[x for _,x in sorted(result)]
    (OUT/"image_pixel_receipts.json").write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n")
    enriched=OUT/"top_candidates_with_checks.json"
    if enriched.exists():
        data=json.loads(enriched.read_text())
        lookup={r["id"]:r for r in result}
        for rec in data:
            if rec["id"] in lookup:
                r=lookup[rec["id"]]
                rec["pixel_vetting_status"]=r["status"]
                rec["pixel_science_bands"]=r.get("science_bands_with_pixels",0)
                rec["image_artefacts_cleared"]=False
        enriched.write_text(json.dumps(data,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("PIXEL_SUMMARY",json.dumps({
      "targets":len(result),
      "targets_with_actual_FITS_pixels":sum(any(
          b.get("status")=="REAL_FITS_SCIENCE_PIXELS_RANGE_READ_NOT_HUMAN_CLEARED"
          for b in r.get("bands",[])) for r in result),
      "full_artifact_cleared":0}),flush=True)
if __name__=="__main__":main()
