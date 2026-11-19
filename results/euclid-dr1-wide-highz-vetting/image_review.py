#!/usr/bin/env python3
"""Best-effort, bounded VIS/NISP pixel-cutout inspection receipts for top Euclid rows.

Image metadata alone is NOT an artefact inspection. Attempt genuine <20 arcsec
FITS cutouts only when the SIA access URL supports IBE cutouts; never fetch a
whole multi-GB mosaic. Save auditable PNG postage stamps where possible.
Detector ghosts/persistence need per-exposure human review and remain UNRESOLVED.
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import struct
import zlib
from urllib.parse import urlsplit, urlunsplit, urlencode, parse_qsl

import numpy as np
from astropy.io import fits
import requests

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"output"
PNGS=OUT/"pixel_previews"
SIA="https://irsa.ipac.caltech.edu/SIA"
MAX_BYTES=8_000_000


def sha(data):
    return hashlib.sha256(data).hexdigest()


def write_png(path, pixels):
    """Minimal grayscale PNG writer; no image library / plotting dependency."""
    pixels=np.ascontiguousarray(pixels, dtype=np.uint8)
    height,width=pixels.shape
    def chunk(tag,data):
        return struct.pack(">I",len(data))+tag+data+struct.pack(">I",zlib.crc32(tag+data)&0xffffffff)
    scan=b"".join(b"\x00"+pixels[i,:].tobytes() for i in range(height))
    encoded=(b"\x89PNG\r\n\x1a\n"+chunk(b"IHDR",struct.pack(">2I5B",width,height,8,0,0,0,0))+
             chunk(b"IDAT",zlib.compress(scan,6))+chunk(b"IEND",b""))
    path.write_bytes(encoded)
    return sha(encoded)


def stamp(data):
    arr=np.asarray(data,dtype=np.float64)
    arr=np.squeeze(arr)
    if arr.ndim != 2 or min(arr.shape)<8:
        raise ValueError("not a 2D cutout")
    finite=np.isfinite(arr)
    if finite.sum()<40:
        raise ValueError("few finite science pixels")
    lo,hi=np.percentile(arr[finite],[1,99.5])
    if hi<=lo:
        raise ValueError("flat cutout")
    gray=np.clip((arr-lo)/(hi-lo),0,1)
    gray=(255*np.sqrt(gray)).astype(np.uint8)
    gray[~finite]=0
    h,w=arr.shape
    central=arr[max(0,h//2-2):min(h,h//2+3),max(0,w//2-2):min(w,w//2+3)]
    source=np.nanmedian(central)
    border=np.concatenate((arr[0,:],arr[-1,:],arr[:,0],arr[:,-1]))
    return gray, {"shape":[int(h),int(w)],
        "finite_fraction":round(float(finite.mean()),4),
        "edge_finite_fraction":round(float(np.isfinite(border).mean()),4),
        "center_median":float(source) if math.isfinite(source) else None,
        "global_p01":float(lo),"global_p995":float(hi),
        "image_artifact_verdict":"NOT_CLEARED_REQUIRES_EXPOSURE_LEVEL_AND_VISUAL_INSPECTION"}


def csv_rows(response):
    if response.status_code!=200:
        raise RuntimeError("SIA_HTTP_"+str(response.status_code))
    content=response.text
    if content.lstrip().startswith("<"):
        from astropy.io.votable import parse_single_table
        table=parse_single_table(io.BytesIO(response.content)).to_table()
        rows=[]
        for record in table:
            rows.append({k:(None if np.ma.is_masked(record[k]) else str(record[k]))
                         for k in table.colnames})
    else:
        rows=list(csv.DictReader(content.splitlines()))
    if rows and "access_url" not in rows[0]:
        raise RuntimeError("SIA_NO_ACCESS_URL_COLUMN")
    return rows


def cutout_url(url,ra,dec):
    u=urlsplit(url)
    if u.scheme!="https" or u.hostname not in ("irsa.ipac.caltech.edu","irsacloud.ipac.caltech.edu") or "/ibe/data/" not in u.path:
        return None
    params=dict(parse_qsl(u.query))
    params.update({"center":f"{ra:.7f},{dec:.7f}deg","size":"20arcsec","gzip":"false"})
    return urlunsplit((u.scheme,u.netloc,u.path,urlencode(params),u.fragment))


def read_bounded(url):
    response=requests.get(url,timeout=(8,30),stream=True)
    response.raise_for_status()
    if int(response.headers.get("Content-Length") or 0)>MAX_BYTES:
        response.close()
        raise RuntimeError("CUTOUT_RESPONSE_TOO_LARGE")
    chunks=[];total=0
    try:
        for part in response.iter_content(chunk_size=65536):
            total+=len(part)
            if total>MAX_BYTES:
                raise RuntimeError("CUTOUT_BYTES_LIMIT_EXCEEDED")
            chunks.append(part)
    finally:
        response.close()
    data=b"".join(chunks)
    if not data.startswith(b"SIMPLE") and not data.startswith(b"\x1f\x8b"):
        raise RuntimeError("NOT_FITS_CUTOUT_RESPONSE")
    return data


def review_one(row, rank):
    ra,dec=float(row["ra"]),float(row["dec"])
    receipt={"rank":rank+1,"id":row["id"],"ra":ra,"dec":dec,
        "service":SIA,"status":"QUERY_ERROR",
        "image_artifact_cleared":False,"detector_persistence_reviewed":False,
        "fields":[]}
    try:
        response=requests.get(SIA,params={"COLLECTION":"euclid_dr1_wide_DpdMerBksMosaic",
            "POS":f"CIRCLE {ra:.7f} {dec:.7f} 0.002",
            "FORMAT":"image/fits","MAXREC":32,"RESPONSEFORMAT":"CSV"},
            timeout=(8,20))
        receipt["query_http_status"]=response.status_code
        receipt["sia_response_sha256"]=sha(response.content)
        found=csv_rows(response)
        receipt["sia_records"]=len(found)
        if not found:
            receipt["status"]="NO_SIA_IMAGE_ROWS_NOT_EVIDENCE_OF_NO_COVERAGE"
            return receipt
        seen=set()
        for entry in found:
            subtype=str(entry.get("dataproduct_subtype") or "").strip().lower()
            if subtype and subtype!="science":
                continue  # background models and weight images are not science cutouts
            band=str(entry.get("energy_bandpassname") or entry.get("obs_title") or "unknown").strip()
            if band not in ("VIS","Y","J","H"):
                continue
            if band in seen:
                continue
            seen.add(band)
            url=str(entry.get("access_url") or "").strip()
            field={"band":band,"access_url":url,"status":"NO_SUPPORTED_CUTOUT_SERVICE"}
            receipt["fields"].append(field)
            url2=cutout_url(url,ra,dec)
            if url2 is None:
                continue
            field["cutout_url"]=url2
            try:
                binary=read_bounded(url2)
                field["cutout_sha256"]=sha(binary)
                field["cutout_bytes"]=len(binary)
                with fits.open(io.BytesIO(binary),memmap=False) as f:
                    image=next((hdu.data for hdu in f if hdu.data is not None and
                        np.squeeze(np.asarray(hdu.data)).ndim==2),None)
                    if image is None:
                        raise ValueError("no 2D science FITS HDU")
                    preview,details=stamp(image)
                name=f"rank_{rank+1:02d}_band_{len(receipt['fields']):02d}.png"
                field["preview_path"]=str(PNGS.relative_to(ROOT)/name)
                field["preview_sha256"]=write_png(PNGS/name,preview)
                field.update(details)
                field["status"]="FITS_PIXELS_SAVED_FOR_HUMAN_REVIEW"
            except Exception as err:
                field["status"]="CUTOUT_ERROR_NOT_NEGATIVE_EVIDENCE"
                field["error"]=type(err).__name__+": "+str(err)[:220]
            if len(receipt["fields"])>=8:
                break
        receipt["status"]=("CUTOUT_PIXELS_AVAILABLE_UNCLEARED" if any(
            v["status"]=="FITS_PIXELS_SAVED_FOR_HUMAN_REVIEW" for v in receipt["fields"])
            else "SIA_METADATA_ONLY_NO_PIXEL_ARTIFACT_VERDICT")
    except Exception as err:
        receipt["error"]=type(err).__name__+": "+str(err)[:240]
    return receipt


def main():
    rows=json.loads((OUT/"top_candidates_with_checks.json").read_text())
    PNGS.mkdir(exist_ok=True)
    results=[]
    with ThreadPoolExecutor(max_workers=3) as pool:
        futures={pool.submit(review_one,row,i):i for i,row in enumerate(rows[:8])}
        for f in as_completed(futures):
            i=futures[f]
            try:
                results.append((i,f.result()))
            except Exception as err:
                results.append((i,{"id":rows[i]["id"],"status":"PIPELINE_ERROR",
                   "error":type(err).__name__+": "+str(err)[:240],
                   "image_artifact_cleared":False}))
    results=[x for _,x in sorted(results)]
    (OUT/"image_pixel_receipts.json").write_text(
        json.dumps(results,indent=2,sort_keys=True,allow_nan=False)+"\n")
    for i,r in enumerate(results):
        rows[i]["pixel_vetting_status"]=r["status"]
        rows[i]["image_level_artifact_status"]="NOT_CLEARED_REQUIRES_EXPOSURE_LEVEL_AND_VISUAL_REVIEW"
    (OUT/"top_candidates_with_checks.json").write_text(
        json.dumps(rows,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("PIXEL_REVIEW",json.dumps({"attempted":len(results),
        "with_pixels":sum(any(v.get("status")=="FITS_PIXELS_SAVED_FOR_HUMAN_REVIEW"
                     for v in r.get("fields",[])) for r in results),
        "cleared":0}),flush=True)


if __name__=="__main__":
    main()
