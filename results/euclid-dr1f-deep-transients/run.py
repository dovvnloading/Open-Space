#!/usr/bin/env python3
"""Frozen three-field DR1-F visit-pair difference survey. Network runs on Actions.

Publisher files are read as FITS byte-range sections, never as whole mosaics.
An unavailable service is an error receipt, not a nondetection. Results are
conditional on visited/masked pixels, not the nominal 53 deg2 EDF footprint.
"""
import csv, hashlib, io, json, math, os, re, sys, time
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlsplit
from xml.etree import ElementTree as ET

import numpy as np
import requests
from scipy import ndimage as ndi
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.path import Path as Polygon

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"output"
FROZEN=json.loads((ROOT/"freeze.json").read_text())
BANDS=("VIS","Y","J","H")
BUCKET="nasa-irsa-euclid-dr1"
SIA="https://irsa.ipac.caltech.edu/SIA"
PIXEL_ARCSEC=0.1
NP=FROZEN["patch_pixels"]
SESSION=requests.Session()
SESSION.headers.update({"User-Agent":"Open-Space-Euclid-DR1F-deep-time-domain/1.0"})
S3FS=None

def utc():
    return datetime.now(timezone.utc).isoformat()
def sha(data):
    return hashlib.sha256(data).hexdigest()
def dump(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2,sort_keys=True,allow_nan=False)+"\n")
def safe(e):
    return type(e).__name__+": "+str(e)[:320]
def fetch(url,**kwargs):
    for k in range(3):
        try:
            x=SESSION.get(url,timeout=(12,55),**kwargs)
            x.raise_for_status()
            return x
        except Exception:
            if k==2: raise
            time.sleep(2**k)
def key_from_url(url):
    """SIA/IRSA object identity, not a guessed mosaic filename."""
    v=urlsplit(url)
    if v.hostname not in ("irsacloud.ipac.caltech.edu","irsa.ipac.caltech.edu"):
        raise ValueError("Untrusted image origin")
    marker="/ibe/data/euclid/dr1/"
    if marker not in v.path: raise ValueError("Not Euclid DR1 IRSA IBE")
    suffix=v.path.split(marker,1)[1]
    if not suffix.startswith("MER/") or not suffix.endswith(".fits"):
        raise ValueError("Not a MER FITS file")
    return "dr1/"+suffix
def https_key(key):
    return "https://nasa-irsa-euclid-dr1.s3.us-east-1.amazonaws.com/"+key

def table_rows(content):
    if content.lstrip().startswith(b"<"):
        doc=ET.fromstring(content)
        statuses=[x.attrib.get("value") for x in doc.findall(".//{*}INFO") if x.attrib.get("name")=="QUERY_STATUS"]
        if statuses and statuses[-1]!="OK":
            raise ValueError("IRSA SIA QUERY_STATUS="+str(statuses))
        table=doc.find(".//{*}TABLE")
        if table is None: raise ValueError("No VOTable TABLE")
        names=[x.attrib["name"] for x in table.findall("./{*}FIELD")]
        rows=[]
        for tr in table.findall(".//{*}TR"):
            rows.append(dict(zip(names,[(x.text or "") for x in tr.findall("./{*}TD")])))
        return rows
    text=content.decode("utf-8-sig")
    return list(csv.DictReader(io.StringIO(text)))

def discover(field):
    name,ra,dec=field["name"],field["ra_deg"],field["dec_deg"]
    receipt={"field":name,"utc":utc(),"endpoint":SIA,"collection":"euclid_dr1_deep_DpdMerBksMosaic",
             "center":[ra,dec],"source_rows":0,"status":"NOT_QUERIED"}
    groups={}
    try:
        response=fetch(SIA,params={"COLLECTION":"euclid_dr1_deep_DpdMerBksMosaic",
            "POS":f"CIRCLE {ra:.7f} {dec:.7f} 0.0025",
            "MAXREC":"10000","RESPONSEFORMAT":"CSV"})
        receipt.update(status="SIA_OK",http_status=response.status_code,
                       response_sha256=sha(response.content),response_bytes=len(response.content))
        rows=table_rows(response.content)
        receipt["source_rows"]=len(rows)
        for row in rows:
            obs=row.get("obs_id","").strip()
            key=(row.get("energy_bandpassname","") or "").strip().upper().replace("NIR-","")
            if key not in BANDS or not obs:continue
            url=row.get("access_url","").strip()
            try: objkey=key_from_url(url)
            except Exception: continue
            tail=objkey.rsplit("/",1)[-1]
            kind=next((kind for kind,needle in (("science","BGSUB-MOSAIC"),("rms","-RMS_"),
                         ("flag","-FLAG_"),("layer","DETECTOR-LAYERING"))
                         if needle in tail),None)
            if kind is None:continue
            if (row.get("dataproduct_subtype") or "").strip().lower()=="science" and kind!="science":continue
            groups.setdefault((key,obs),{})[kind]=objkey
    except Exception as e:
        receipt.update(status="SIA_ERROR_NOT_NO_COVERAGE",error=safe(e))
    visits=[]
    for band in BANDS:
        eligible=[(obs,paths) for (b,obs),paths in groups.items() if b==band and "science" in paths]
        eligible.sort()
        receipt[band+"_science_obs_total"]=len(eligible)
        for obs,paths in eligible[:FROZEN["max_visits_per_band"]]:
            visits.append({"band":band,"obs_id":obs,"products":paths,
                           "selection":"FIRST_LEXICAL_SCIENCE_OBS_IDS_BEFORE_PIXELS"})
    receipt["selected_visits"]=len(visits)
    return receipt,visits

def layer_receipt(v,ra,dec):
    out={"obs_id":v["obs_id"],"band":v["band"],
         "status":"MISSING_LAYER_NOT_INDEPENDENTLY_VERIFIED","exposure_ids":[]}
    key=v["products"].get("layer")
    if not key:return out
    try:
        resp=fetch(https_key(key))
        with fits.open(io.BytesIO(resp.content),memmap=False) as h:
            tbl=h[1].data
            xy=(math.radians(ra),math.radians(dec))
            sub=[r for r in tbl if Polygon(np.asarray(r["POLYGON"]).reshape(-1,2)).contains_point(xy)]
            exposures=[{"id":str(r["CALIB_IMAGE"]).strip(),"extension":str(r["EXTNAME"]).strip(),
                        "mjd_float32":float(r["MJD-OBS"])} for r in sub]
        out.update(status="LAYER_POINT_COVERAGE" if exposures else "LAYER_POINT_NOT_COVERED",
                   source_key=key,source_bytes=len(resp.content),source_sha256=sha(resp.content),
                   exposure_ids=sorted(set(x["id"] for x in exposures)),
                   exposures=exposures,layer_total_rows=len(tbl),
                   center_covering_rows=len(sub),
                   median_mjd=float(np.median([x["mjd_float32"] for x in exposures])) if exposures else None)
    except Exception as e:
        out["status"]="LAYER_READ_ERROR_NOT_CLEARED";out["error"]=safe(e)
    return out

def inventory():
    manifest={"freeze_file_sha256":sha((ROOT/"freeze.json").read_bytes()),
        "method_version":"2026-11-19-psf-diff-v1","created_utc":utc(),
        "release":"IRSA Euclid DR1-Foundation (2026-11-12)",
        "fields":[],"visits":[],"receipts":[]}
    for field in FROZEN["fields"]:
        receipt,visits=discover(field)
        manifest["receipts"].append(receipt)
        for v in visits:
            v.update(field=field["name"],center=[field["ra_deg"],field["dec_deg"]])
            l=layer_receipt(v,*v["center"])
            v["layer"]=l
            manifest["visits"].append(v)
        manifest["fields"].append(field)
    for field in FROZEN["fields"]:
        for band in BANDS:
            vv=[v for v in manifest["visits"] if v["field"]==field["name"] and v["band"]==band]
            for i in range(len(vv)):
                for j in range(i+1,len(vv)):
                    a,b=vv[i]["layer"],vv[j]["layer"]
                    independent=(a.get("status")=="LAYER_POINT_COVERAGE" and
                                 b.get("status")=="LAYER_POINT_COVERAGE" and
                                 not set(a["exposure_ids"]).intersection(b["exposure_ids"]) and
                                 abs(a["median_mjd"]-b["median_mjd"])>1.0)
                    vv[j].setdefault("pair_independence",{})[vv[i]["obs_id"]]=bool(independent)
    dump(OUT/"frozen_manifest.json",manifest)
    print("INVENTORY",[(r["field"],r["status"],r["selected_visits"]) for r in manifest["receipts"]],flush=True)

def s3():
    global S3FS
    if S3FS is None:
        import s3fs
        S3FS=s3fs.S3FileSystem(anon=True,client_kwargs={"region_name":"us-east-1"})
    return S3FS

def target_wcs(ra,dec):
    w=WCS(naxis=2);w.wcs.crval=[ra,dec];w.wcs.crpix=[(NP+1)/2,(NP+1)/2]
    w.wcs.cdelt=[-PIXEL_ARCSEC/3600,PIXEL_ARCSEC/3600];w.wcs.ctype=["RA---TAN","DEC--TAN"]
    return w

def image(key,wout,order=1):
    """Range-read a rectangle then reproject to a fixed tangent WCS.
    No access to full FITS .data, no full mosaic download.
    """
    uri="s3://"+BUCKET+"/"+key
    with fits.open(uri,use_fsspec=True,
                   fsspec_kwargs={"anon":True,"client_kwargs":{"region_name":"us-east-1"}},
                   memmap=False) as h:
        hd=next((x for x in h if x.header.get("NAXIS")==2 and x.header.get("NAXIS1",0)>100),None)
        if hd is None:raise ValueError("No 2D FITS HDU")
        w=WCS(hd.header).celestial
        ys,xs=np.mgrid[0:NP,0:NP]
        ra,dec=wout.all_pix2world(xs.astype(float),ys.astype(float),0)
        xx,yy=w.all_world2pix(ra,dec,0)
        finite=np.isfinite(xx)&np.isfinite(yy)
        if finite.sum()==0:raise ValueError("Target WCS not finite")
        x0=max(0,int(np.floor(np.min(xx[finite])))-4); x1=min(int(hd.header["NAXIS1"]),int(np.ceil(np.max(xx[finite])))+5)
        y0=max(0,int(np.floor(np.min(yy[finite])))-4); y1=min(int(hd.header["NAXIS2"]),int(np.ceil(np.max(yy[finite])))+5)
        if x1<=x0 or y1<=y0 or (x1-x0)*(y1-y0)>4*NP*NP:
            raise ValueError("Input pixel bounds invalid")
        block=np.asarray(hd.section[y0:y1,x0:x1])
        # For bit flags use float64 with nearest-neighbour: float32 loses low
        # bits of large integer masks, silently corrupting the INVALID bit.
        pixel_dtype=np.float64 if order==0 else np.float32
        sampled=ndi.map_coordinates(block.astype(pixel_dtype),[yy-y0,xx-x0],order=order,
                     mode="constant",cval=np.nan if order==1 else 2147483647,
                     prefilter=False).astype(pixel_dtype)
        hash_dtype="<f8" if order==0 else "<f4"
        receipt={"key":key,"bounds_xyxy":[x0,y0,x1,y1],"block_shape":list(block.shape),
                 "sampled_sha256":sha(np.ascontiguousarray(sampled,dtype=hash_dtype).tobytes()),
                 "valid_coordinate_fraction":float(np.isfinite(sampled).mean()),
                 "header_wcs_crval":[float(a) for a in w.wcs.crval],
                 "interpolated":order==1}
        return sampled,receipt

def read_visit(v,field,w):
    arrays={};receipts=[];errors={}
    for kind,order in (("science",1),("rms",1),("flag",0)):
        key=v["products"].get(kind)
        if not key:
            errors[kind]="PUBLISHER_PRODUCT_MISSING"
            continue
        try:
            arr,record=image(key,w,order=order)
            arrays[kind]=arr
            record.update(kind=kind)
            receipts.append(record)
        except Exception as e:errors[kind]=safe(e)
    return arrays,{"obs_id":v["obs_id"],"band":v["band"],"pixels":receipts,"errors":errors}

def good_mask(a):
    if not all(k in a for k in ("science","rms","flag")):
        return np.zeros((NP,NP),bool)
    flag=a["flag"].astype(np.int64)
    # Euclid DR1 DPDD VIS/NIR calibrated-frame definitions: bit 0 is the
    # INVALID bit, activated for critical defects (hot/cold/saturated/cosmic).
    # MER propagates flags as OR; requiring every *informational* bit to be
    # zero over many stacked inputs destroys actual sky coverage.
    # For VIS additionally exclude direct SAT, COSMIC and GHOST flags.
    # Additional warning bits and exposure-level artefacts still require a
    # separate detector-level review before any discovery is validated.
    fatal=1 | (8|16|32 if a["band"]=="VIS" else 0)
    return (np.isfinite(a["science"]) & np.isfinite(a["rms"]) &
        (a["rms"]>0)&np.isfinite(a["flag"]) & ((flag&fatal)==0))

def robust(x):
    arr=np.asarray(x);arr=arr[np.isfinite(arr)]
    if len(arr)<40:return float("nan")
    m=float(np.median(arr))
    return max(float(1.4826*np.median(np.abs(arr-m))),1e-12)

def sources(im,mask,rms):
    """Prior North v5 positive local-peak/RMS rule, not a global image MAD.
    Empirical compact-source PSF proxy; stellar identities unverified.
    """
    good=mask & np.isfinite(rms) & (rms>0)
    if good.sum()<100:return []
    sky=float(np.median(im[good]))
    base=np.where(good,im,sky)
    background=ndi.gaussian_filter(base,6)
    signal=base-background
    maxmap=ndi.maximum_filter(signal,size=5)
    # Match the pre-existing V5 four-pixel aperture validity rule, not an
    # 11-pixel erosion that removes all stars in sparse flagged mosaics.
    m=good & (signal>8*rms) & (signal==maxmap) & ndi.binary_erosion(good,iterations=4)
    yy,xx=np.where(m)
    found=[]
    for x,y in zip(xx,yy):
        if x<16 or y<16 or x>=NP-16 or y>=NP-16:continue
        sub=im[y-9:y+10,x-9:x+10]
        Y,X=np.mgrid[-9:10,-9:10]
        ring=(X*X+Y*Y>=40)
        local_sky=float(np.median(sub[ring]))
        stamp=np.maximum(sub-local_sky,0)
        ss=stamp.sum()
        if ss<=0:continue
        cx=float((stamp*X).sum()/ss)
        cy=float((stamp*Y).sum()/ss)
        r2=float((stamp*((X-cx)**2+(Y-cy)**2)).sum()/ss)
        sigma=math.sqrt(max(r2/2,.25))
        if .8<sigma<10 and abs(cx)<2 and abs(cy)<2:
            found.append((int(x),int(y),float(sigma)))
    return found

def psf_calibrate(a,b,mask,rms_a,rms_b):
    sa=sources(a,mask,rms_a);sb=sources(b,mask,rms_b)
    psfa=float(np.median([p[2] for p in sa])) if len(sa)>=6 else None
    psfb=float(np.median([p[2] for p in sb])) if len(sb)>=6 else None
    # Never invent a PSF if too few usable calibration sources exist.
    if psfa is None or psfb is None:
        return None,{"status":"INSUFFICIENT_PSF_STARS","sources_a":len(sa),"sources_b":len(sb)}
    pairs=[];r=max(5,int(math.ceil(2.5*max(psfa,psfb))))
    yy,xx=np.mgrid[-r:r+1,-r:r+1];ap=(xx*xx+yy*yy<=r*r)
    for x,y,_ in sa:
        near=[(xx,yy) for xx,yy,_ in sb if (xx-x)**2+(yy-y)**2<=2.5**2]
        if len(near)!=1 or x<r+3 or y<r+3 or x>=NP-r-3 or y>=NP-r-3:continue
        x2,y2=near[0]
        if x2<r+3 or y2<r+3 or x2>=NP-r-3 or y2>=NP-r-3:continue
        f1=float(np.sum(a[y-r:y+r+1,x-r:x+r+1][ap]))
        f2=float(np.sum(b[y2-r:y2+r+1,x2-r:x2+r+1][ap]))
        if f1>0 and f2>0:pairs.append(f2/f1)
    if len(pairs)<6:
        return None,{"status":"INSUFFICIENT_MATCHED_STARS","matched":len(pairs),"sources_a":len(sa),"sources_b":len(sb)}
    rat=float(np.median(pairs))
    spread=robust(pairs)
    if not (.2<=rat<=5 and np.isfinite(spread) and spread/rat<.25):
        return None,{"status":"UNSTABLE_PHOTOMETRIC_SCALE","n":len(pairs),"scale":rat,"scatter":spread}
    return (psfa,psfb,rat),{"status":"EMPIRICAL_COMPACT_SOURCE_PSF_AND_FLUX_PROXY","n_matched":len(pairs),
                          "psf_sigmas_px":[psfa,psfb],"scale_b_over_a":rat,
                          "robust_fractional_scale_scatter":spread/rat}

def matched_diff(a,b,mask,cal):
    s1,s2,ratio=cal
    common=max(s1,s2)
    # Match the narrower PSF to the broader PSF, symmetrically per pair.
    c1=ndi.gaussian_filter(np.where(mask,a,0),math.sqrt(max(0,common**2-s1**2))) if common>s1+.05 else np.where(mask,a,0)
    c2=ndi.gaussian_filter(np.where(mask,b,0),math.sqrt(max(0,common**2-s2**2))) if common>s2+.05 else np.where(mask,b,0)
    difference=c2-ratio*c1
    difference-=float(np.median(difference[mask]))
    clean=ndi.binary_erosion(mask,iterations=max(12,int(3*common)))
    filt=ndi.gaussian_filter(difference,common)
    scatter=robust(filt[clean])
    return difference,filt,clean,scatter,common

def peaks(filt,mask,scatter,width):
    if not np.isfinite(scatter) or scatter<=0:return []
    candidates=[]
    for sign in (-1,1):
        image=filt*sign
        high=(image==ndi.maximum_filter(image,size=max(13,int(3*width)*2+1)))
        threshold=FROZEN["threshold_sigma"]*scatter
        yy,xx=np.where(high & (image>threshold) & mask)
        for x,y in zip(xx,yy):
            if x<17 or y<17 or x>=NP-17 or y>=NP-17:continue
            cube=np.maximum(sign*filt[y-12:y+13,x-12:x+13],0)
            Y,X=np.mgrid[-12:13,-12:13]
            ss=float(cube.sum())
            if ss<=0:continue
            xm=float((cube*X).sum()/ss);ym=float((cube*Y).sum()/ss)
            qxx=float((cube*(X-xm)**2).sum()/ss)
            qyy=float((cube*(Y-ym)**2).sum()/ss)
            qxy=float((cube*(X-xm)*(Y-ym)).sum()/ss)
            eig=np.linalg.eigvalsh([[qxx,qxy],[qxy,qyy]])
            elongation=float(math.sqrt(max(eig[1],.01)/max(eig[0],.01)))
            if elongation>1.7 or abs(xm)>3 or abs(ym)>3:continue
            candidates.append({"x":int(x),"y":int(y),"sign":int(sign),
                "sigma_empirical":float(image[y,x]/scatter),
                "shape_axis_ratio":elongation})
    return sorted(candidates,key=lambda x:-x["sigma_empirical"])[:60]

def injections(filt,clean,scatter,width):
    """Add PSF-shaped synthetic transients to a held-fixed null difference.
    Re-run the identical peak finder, not a known-location threshold oracle.
    """
    rng=np.random.default_rng(FROZEN["random_seed"])
    ys,xs=np.where(ndi.binary_erosion(clean,iterations=18))
    if len(xs)<100:return {"status":"NO_VALID_BLANK_POSITIONS"}
    trials=[]
    # Define flux on the PSF-filtered DIFFERENCE map in control-noise units.
    # Inject into the original difference prior to filtering.
    sigma=max(.8,float(width))
    Y,X=np.mgrid[-int(5*sigma):int(5*sigma)+1,-int(5*sigma):int(5*sigma)+1]
    g=np.exp(-(X*X+Y*Y)/(2*sigma*sigma))
    response=float(ndi.gaussian_filter(g,sigma)[len(g)//2,len(g)//2])
    for strength in FROZEN["injection_levels_sigma"]:
        recovered=0;valid=0
        order=rng.permutation(len(xs))
        for k in order[:min(40,len(order))]:
            x=int(xs[k]);y=int(ys[k])
            d=int(5*sigma)
            if min(x,y,NP-x,NP-y)<d+15:continue
            # Spatially local application of the same matched-PSF filter.
            # Existing observed filter map supplies the undisturbed baseline.
            injected=filt.copy()
            h=ndi.gaussian_filter(g,sigma)*float(strength)*scatter/response
            injected[y-d:y+d+1,x-d:x+d+1]+=h
            detected=peaks(injected,clean,scatter,sigma)
            valid+=1
            if any(q["sign"]==1 and abs(q["x"]-x)<=3 and abs(q["y"]-y)<=3 for q in detected):
                recovered+=1
        trials.append({"input_filtered_peak_sigma":strength,"placed":valid,"recovered":recovered,
                       "efficiency":recovered/valid if valid else None})
    return {"status":"CONDITIONAL_DIFFERENCE_MAP_INJECTIONS",
            "trials":trials,"nuisance_model":"empirical control-noise, fixed photometric scale and PSF",
            "not_included":"raw exposure persistence, time-dependent host subtraction, missed mosaic coverage"}

def analyse_field(field,manifest):
    n=field["name"];ra,dec=field["ra_deg"],field["dec_deg"]
    visits=[x for x in manifest["visits"] if x["field"]==n]
    output={"field":n,"utc":utc(),"center":[ra,dec],"patch_pixels":NP,
            "nominal_square_deg":float((NP*PIXEL_ARCSEC/3600)**2),
            "band_results":[],"pixel_receipts":[],"raw_triggers":[],"status":"INCOMPLETE"}
    w=target_wcs(ra,dec)
    for band in BANDS:
        selected=[v for v in visits if v["band"]==band]
        if len(selected)<2:
            output["band_results"].append({"band":band,"status":"FEWER_THAN_TWO_SIA_VISITS"})
            continue
        cache=[]
        for visit in selected:
            arrays,receipt=read_visit(visit,field,w)
            arrays["band"]=band
            output["pixel_receipts"].append(receipt)
            cache.append((visit,arrays))
        for i in range(len(cache)):
            for j in range(i+1,len(cache)):
                v1,a=cache[i];v2,b=cache[j]
                pair={"band":band,"obs_ids":[v1["obs_id"],v2["obs_id"]],
                      "status":"NOT_SCREENED","validated_discovery":False}
                good=bool(v2.get("pair_independence",{}).get(v1["obs_id"]))
                pair["disjoint_detector_layers_center"]=good
                if not good:
                    pair["status"]="INDEPENDENCE_NOT_VERIFIED";output["band_results"].append(pair);continue
                ma,mb=good_mask(a),good_mask(b)
                mask=ma & mb
                pair["unmasked_pixels"]=int(mask.sum())
                pair["mask_provenance"]="DR1 DPDD bit0=INVALID; VIS also SAT/COSMIC/GHOST excluded; uint32 bit integrity retained in float64 sampling"
                pair["individual_good_pixels"]=[int(ma.sum()),int(mb.sum())]
                pair["joint_area_deg2"]=float(mask.sum())*(PIXEL_ARCSEC/3600)**2
                if mask.sum()<NP*NP*.15:
                    pair["status"]="INSUFFICIENT_JOINT_UNMASKED_PIXELS"
                    output["band_results"].append(pair);continue
                noise=robust(np.concatenate((a["science"][mask][::50],b["science"][mask][::50])))
                if not np.isfinite(noise):
                    pair["status"]="INVALID_SCIENCE_NOISE";output["band_results"].append(pair);continue
                cal,calrec=psf_calibrate(a["science"],b["science"],mask,a["rms"],b["rms"])
                pair["calibration"]=calrec
                pair["calibration"]["source_detection"]="reused prior Euclid North v5 signal>8*RMS local-source screen; no independent stellar IDs"
                if cal is None:
                    pair["status"]="NO_STABLE_PSF_OR_FLUX_CALIBRATION";output["band_results"].append(pair);continue
                diff,filt,clean,scatter,width=matched_diff(a["science"],b["science"],mask,cal)
                pair["control_scatter"]=float(scatter) if np.isfinite(scatter) else None
                pair["detection_area_deg2"]=float(clean.sum())*(PIXEL_ARCSEC/3600)**2
                pair["difference_pixel_sha256"]=sha(np.ascontiguousarray(diff,dtype="<f4").tobytes())
                hits=peaks(filt,clean,scatter,width)
                skyra,skydec=w.all_pix2world(
                    np.array([p["x"] for p in hits]),np.array([p["y"] for p in hits]),0)
                for q,x,y in zip(hits,skyra,skydec):
                    q.update(field=n,band=band,ra_deg=float(x),dec_deg=float(y),
                             pair_obs_ids=pair["obs_ids"],status="UNVETTED_DIFFERENCE_TRIGGER")
                pair["raw_trigger_count"]=len(hits)
                pair["top_triggers"]=hits[:8]
                output["raw_triggers"]+=hits
                pair["injection"]=injections(filt,clean,scatter,width)
                pair["status"]="EMPIRICALLY_SCREENED_UNCONFIRMED"
                output["band_results"].append(pair)
                print("PAIR",n,band,pair["obs_ids"],"area",round(pair["detection_area_deg2"],6),
                      "raw",len(hits),flush=True)
    output["status"]="FIELDS_SCREENED_WITH_LIMITS" if any(
        b["status"]=="EMPIRICALLY_SCREENED_UNCONFIRMED" for b in output["band_results"]) else "NO_CALIBRATED_PAIR"
    dump(OUT/(n+".json"),output)

def external_one(hit):
    ra,dec=hit["ra_deg"],hit["dec_deg"]
    rec={"field":hit["field"],"ra_deg":ra,"dec_deg":dec,"band":hit["band"],
         "discovery_status":"NOT_VALIDATED","services":{}}
    queries={
      "gaia_dr3":("https://gea.esac.esa.int/tap-server/tap/sync",{"REQUEST":"doQuery","LANG":"ADQL",
        "FORMAT":"csv","QUERY":("SELECT TOP 12 source_id,ra,dec,pmra,pmdec,phot_variable_flag "
          "FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec),"
          f"CIRCLE('ICRS',{ra:.8f},{dec:.8f},0.00084))")}),
      "euclid_phz":("https://irsa.ipac.caltech.edu/TAP/sync",{"REQUEST":"doQuery","LANG":"ADQL",
        "FORMAT":"csv","QUERY":("SELECT TOP 8 * FROM euclid_dr1_phz_photo_z_deep_survey WHERE "
          "1=CONTAINS(POINT('ICRS',ra,dec),"
          f"CIRCLE('ICRS',{ra:.8f},{dec:.8f},0.00084))")}),
      "tns_public":("https://www.wis-tns.org/search",{"ra":ra,"decl":dec,"radius":3,"units":"arcsec"}),
      "jpl_small_body_ident":("https://ssd-api.jpl.nasa.gov/sb_ident.api",{
           "obs-time":"2025-01-01","lat":0,"lon":0,"radius":"0.01","mag-limit":24})
    }
    for name,(url,params) in queries.items():
        try:
            r=SESSION.get(url,params=params,timeout=(8,19))
            payload=r.content[:150000]
            rec["services"][name]={"http_status":r.status_code,"response_sha256":sha(payload),
                                  "bytes_sampled":len(payload),"query_url":r.url.split("?")[0],
                                  "status":"RESPONSE_REQUIRES_HUMAN_CROSSMATCH"}
            if r.status_code!=200:rec["services"][name]["status"]="QUERY_ERROR_NOT_NONDTECTION"
            if name=="jpl_small_body_ident":
                rec["services"][name]["status"]="NOT_POSITION_EPOCH_VALIDATED"
        except Exception as e:
            rec["services"][name]={"status":"QUERY_ERROR_NOT_NONDTECTION","error":safe(e)}
    rec["literature_status"]="TNS_RESPONSE_NOT_OBJECT_LEVEL_LITERATURE_CLEARANCE"
    rec["moving_object_status"]="NOT_CLEARED_MPC_HORIZONS_EPHEMERIS_REQUIRED"
    rec["host_and_redshift_status"]="NOT_VALIDATED"
    return rec

def summary(manifest):
    fields=[json.loads((OUT/(f["name"]+".json")).read_text()) for f in FROZEN["fields"]]
    screened=[p for f in fields for p in f["band_results"] if p["status"]=="EMPIRICALLY_SCREENED_UNCONFIRMED"]
    raw=[c for f in fields for c in f["raw_triggers"]]
    # Non-independent bands or pairs are not counted as a single discovery.
    ranked=sorted(raw,key=lambda x:-x["sigma_empirical"])
    follow=[external_one(x) for x in ranked[:min(12,len(ranked))]]
    dump(OUT/"external_receipts.json",follow)
    area=sum(float(f["nominal_square_deg"]) for f in fields)
    effective=sum(p["detection_area_deg2"] for p in screened)
    s={"status":"COMPLETE_LIMITED_FIELD_SEARCH","fields_requested":len(fields),
       "fields_with_calibrated_pairs":sum(any(p["status"]=="EMPIRICALLY_SCREENED_UNCONFIRMED"
           for p in f["band_results"]) for f in fields),
       "selected_visit_band_observations":len(manifest["visits"]),
       "independent_calibrated_visit_pairs":len(screened),
       "nominal_patch_area_deg2":area,
       "sum_band_pair_areas_deg2_NOT_UNION":effective,
       "raw_difference_triggers":len(raw),"triggers_queried_for_external_checks":len(follow),
       "validated_astrophysical_candidates":0,
       "candidate_disposition":"ALL_UNVALIDATED_UNTIL_EXPOSURE_ARTIFACT_MOTION_HOST_AND_TNS_CHECKS",
       "rate_status":"NO_EMPIRICAL_TRANSIENT_RATE_LIMIT_WITHOUT_COMPLETE_FOLLOWUP_AND_CONTROL_TIME",
       "method":"PSF-matched Gaussian-convolved visit-pair differencing on S3 FITS subrectangles",
       "caveat":"VIS/NISP band-pair areas overlap; do not sum them as sky area or use nominal EDF footprints",
       "created_utc":utc()}
    dump(OUT/"summary.json",s)
    print("SUMMARY",json.dumps(s),flush=True)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    manifest_path=OUT/"frozen_manifest.json"
    if not manifest_path.exists():
        inventory();return
    manifest=json.loads(manifest_path.read_text())
    if manifest["freeze_file_sha256"]!=sha((ROOT/"freeze.json").read_bytes()):
        raise RuntimeError("Freeze changed after inventory; stop, do not silently reselect")
    for field in FROZEN["fields"]:
        dest=OUT/(field["name"]+".json")
        if not dest.exists():
            analyse_field(field,manifest)
            return
    if not (OUT/"summary.json").exists():summary(manifest)
    else:print("COMPLETE (no reprocessing)")
if __name__=="__main__":
    main()
