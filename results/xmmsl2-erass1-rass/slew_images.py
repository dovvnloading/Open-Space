"""Re-examine the release-named soft image in three individual ESA slew packages."""
import argparse
import gzip
import hashlib
import json
import re
import shutil
import tarfile
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

from full_only_review import XSA, digest, plain

TARGETS = ("XMMSL2-total-25537", "XMMSL2-total-10038", "XMMSL2-total-21146")


def unpack(obs, image, work):
    params = {"obsno": obs, "level": "PPS", "name": image}
    receipt = {"url": XSA, "params": params, "selected_files": [], "members_seen": 0,
               "similar_images": [], "event_files_seen": [], "complete_listing": False}
    matched = {}
    try:
        with requests.get(XSA, params=params, stream=True, timeout=(15,110)) as r:
            receipt.update(final_url=r.url, http=r.status_code, content_length=r.headers.get("Content-Length"),
                           content_type=r.headers.get("Content-Type"), content_disposition=r.headers.get("Content-Disposition"))
            r.raise_for_status()
            r.raw.decode_content = True
            with tarfile.open(fileobj=r.raw,mode="r|*",bufsize=1024*256) as archive:
                for m in archive:
                    if not m.isfile():
                        continue
                    name = Path(m.name).name
                    canonical = re.sub(r"\.gz$", "", name,flags=re.I).upper()
                    receipt["members_seen"] += 1
                    if ("IMAGE_6" in canonical or "IMAGE6" in canonical) and len(receipt["similar_images"]) < 180:
                        receipt["similar_images"].append(name)
                    event = bool(re.search(r"PIEVLI|EVENLI|EVLIST|EVLI0",canonical) and "PN" in canonical and obs in canonical)
                    if event and len(receipt["event_files_seen"])<100:
                        receipt["event_files_seen"].append(name)
                    role = ("image" if canonical==image.upper() else "events" if event else "")
                    if not role or role in matched:
                        continue
                    dest = work/role
                    stream = archive.extractfile(m)
                    if stream is None:
                        continue
                    with dest.open("wb") as f:
                        shutil.copyfileobj(stream, f, length=1024*1024)
                    if dest.stat().st_size != m.size:
                        raise ValueError("Incomplete tar member "+name)
                    matched[role] = dest
                    receipt["selected_files"].append({"role":role,"member":m.name,"bytes":m.size,"sha256":digest(dest)})
                    if set(matched)=={"image","events"}:
                        break
                else:
                    receipt["complete_listing"]=True
        receipt["status"]="exact release member extracted" if "image" in matched else "archive enumerated; exact release member missing"
    except Exception as exc:
        receipt.update(status="incomplete ESA archive transport or packaging",error=str(exc)[:440])
    return matched,receipt


def offsets(origin, grid):
    x,y=origin.spherical_offsets_to(grid.icrs)
    return x.arcsec,y.arcsec


def describe_image(path, pos, rivals):
    with fits.open(path,memmap=False) as hdus:
        h=next(h for h in hdus if h.header.get("NAXIS")==2 and WCS(h.header).has_celestial)
        grid=WCS(h.header).celestial
        scales=proj_plane_pixel_scales(grid)*3600
        xc,yc=map(float,grid.world_to_pixel(pos))
        if not np.all(np.isfinite(scales)) or np.any((scales<.1)|(scales>30)):
            raise ValueError("Unexpected native pixel size")
        im=h.data
        sy,sx=im.shape
        if not (-.5<=xc<sx-.5 and -.5<=yc<sy-.5):
            raise ValueError("Release source outside supplied WCS bounds")
        r=int(np.ceil(100/min(scales)))
        x0,y0=int(round(xc)),int(round(yc))
        y,x=np.mgrid[max(0,y0-r):min(sy,y0+r+1),max(0,x0-r):min(sx,x0+r+1)]
        weights=np.asarray(im[y,x],dtype=float)
        valid=np.isfinite(weights)&(weights>=0)
        if not valid.any() or np.any(np.abs(weights[valid]-np.rint(weights[valid]))>1e-4):
            raise ValueError("Pixel values do not have integer-count semantics")
        dx,dy=offsets(pos,grid.pixel_to_world(x,y))
        d=np.hypot(dx,dy)
        bg=(d>=45)&(d<90)&valid
        pixarea=float(np.prod(scales))
        bgarea=float(bg.sum())*pixarea
        st={"size_pixels":[sx,sy],"source_xy_zero_based":[xc,yc],"pixel_scale_arcsec":list(scales),
            "bunit":str(h.header.get("BUNIT","")),"exposure_header_s":h.header.get("EXPOSURE"),
            "annulus_45_90_counts":int(weights[bg].sum()),"annulus_fraction_of_planar_area":bgarea/(np.pi*(90**2-45**2)),
            "apertures":[],"interpretation":"Geometric area on available image pixels only. A valid image pixel can have partial or zero exposure. Net counts and sparse shape measures have no exposure correction or selection-error correction."}
        for radius in (15.,25.,30.,40.):
            mask=(d<=radius)&valid
            if not mask.any():
                continue
            w=weights[mask]
            n=int(w.sum())
            aperture={"radius_arcsec":radius,"gross_counts":n,"occupied_pixels":int((w>0).sum()),
                      "brightest_pixel_counts":int(w.max()),"fraction_of_planar_area":float(mask.sum()*pixarea/(np.pi*radius**2))}
            if bgarea>0:
                b=float(weights[bg].sum()*mask.sum()/bg.sum())
                aperture.update(local_background_counts=b,geometric_net_counts=n-b)
            if n>0:
                a=np.column_stack((dx[mask],dy[mask]))
                centre=np.average(a,axis=0,weights=w)
                major,minor=np.sqrt(np.maximum(0,np.linalg.eigvalsh(((a-centre).T*w)@(a-centre)/n)))[::-1]
                p=pos.spherical_offsets_by(centre[0]*u.arcsec,centre[1]*u.arcsec)
                aperture.update(centroid_ra=p.ra.deg,centroid_dec=p.dec.deg,offset_east_north_arcsec=centre.tolist(),
                                rms_major_minor_arcsec=[major,minor],optical_centroid_separations=[
                    {"host":r["host"],"sep_arcsec":p.separation(spot(r)).arcsec} for r in rivals
                    if r.get("host_ra") is not None and r.get("host_dec") is not None])
            st["apertures"].append(aperture)
        return st


def spot(r):
    return SkyCoord(float(r["host_ra"])*u.deg,float(r["host_dec"])*u.deg)


def event_positions(h, good):
    names=h.columns.names
    if {"RA","DEC"}.issubset(names):
        units=[str(h.columns[k].unit).lower() for k in ("RA","DEC")]
        if not all(s in {"deg","degree","degrees"} for s in units):
            raise ValueError("Unknown event RA/DEC units")
        return SkyCoord(h.data["RA"][good],h.data["DEC"][good],unit="deg")
    w=WCS(naxis=2)
    for i,k in enumerate(("X","Y")):
        idx=names.index(k)+1
        dim="XY"[i]
        for key,prefix,attr in (("TCRPX","CRPX","crpix"),("TCRVL","CRVL","crval"),("TCDLT","CDLT","cdelt"),("TCTYP","CTYP","ctype")):
            val=h.header.get(f"{key}{idx}",h.header.get(f"REF{dim}{prefix}"))
            if val is None:
                raise ValueError("Unavailable event list celestial WCS "+key+dim)
            getattr(w.wcs,attr)[i]=val
        if h.header.get(f"TCROT{idx}",0.)!=0 or any(re.match(rf"T[PC]{idx}_\d+",t) for t in h.header):
            raise ValueError("Rotated event table WCS requires a separate calibration")
    if not w.has_celestial:
        raise ValueError("No recognized event-sky coordinate pair")
    world=w.all_pix2world(h.data["X"][good],h.data["Y"][good],1)
    return SkyCoord(world[0],world[1],unit="deg")


def analyze_events(path,pos):
    with fits.open(path,memmap=False) as hdus:
        h=next(h for h in hdus if h.name=="EVENTS")
        data=h.data
        names=h.columns.names
        good=np.isfinite(data["TIME"]) & (data["PI"]>=200) & (data["PI"]<=2000)
        pi_unit=str(h.columns["PI"].unit).lower()
        if pi_unit not in ("ev","channel","chan","adu","none"):
            raise ValueError("Unexpected PI unit "+pi_unit)
        cuts=["200 <= PI <= 2000 (nominal pn PI in eV)"]
        if "FLAG" in names:
            good &= data["FLAG"]==0
            cuts.append("FLAG=0")
        if "PATTERN" in names:
            good &= data["PATTERN"]<=4
            cuts.append("PATTERN<=4")
        sky=event_positions(h,good)
        dx,dy=offsets(pos,sky)
        mask=np.hypot(dx,dy)<=30.
        x,y=dx[mask],dy[mask]
        t=np.asarray(data["TIME"][good][mask],dtype=float)
        ix=np.argsort(t)
        temporal=[]
        for label,idx in zip(("first_half","second_half"),np.array_split(ix,2)):
            if len(idx):
                temporal.append({"part":label,"counts":len(idx),"time_min":float(t[idx].min()),"time_max":float(t[idx].max()),
                                 "centroid_east_north_arcsec":[x[idx].mean(),y[idx].mean()]})
        return {"status":"native event-list positions and time order inspected","cuts":cuts,"pi_unit":pi_unit,
                "time_unit":str(h.columns["TIME"].unit),"counts_in_30as":len(t),"time_halves":temporal,
                "events_in_30as":[[float(t[i]),float(x[i]),float(y[i])] for i in ix],
                "scope":"Ordered photons are from this same observation. Small count differences between time halves or image apertures do not establish motion."}


def rows(path):
    op=gzip.open if path.suffix==".gz" else open
    with op(path,"rt") as f:
        return [json.loads(line) for line in f if line.strip()]


def coordinate_alias(name):
    m=re.search(r"J(\d\d)(\d\d)(\d\d(?:\.\d+)?)([+-])(\d\d)(\d\d)(\d\d(?:\.\d+)?)",name)
    if m:
        hh,mm,ss,sign,dd,dm,ds=m.groups()
        return SkyCoord(f"{hh}h{mm}m{ss}s",f"{sign}{dd}d{dm}m{ds}s")


def aliases_and_history(out,sid,gal,raw,field,obspos):
    rivals=[{"host":r["host"],"host_ra":r["host_ra"],"host_dec":r["host_dec"],"local_pcc":r["local_pcc"],
             "sep_high_arcsec":r["host_high_sep_arcsec"],"origin":r["galaxy_support"]} for r in gal]
    near=raw.get("alternatives_within_35arcsec",{}).get("NED",[])
    original="GALEXASC J200452.71-745401.0"
    uv=next((x for x in near if original.replace(" ","")==str(x.get("preferred_name","")).replace(" ","")),None)
    uvpos=SkyCoord(float(uv["ra"])*u.deg,float(uv["dec"])*u.deg) if uv else coordinate_alias(original)
    if uvpos is None:
        raise ValueError("Invalid GALEX coordinate designation")
    distances=[{"host":r["host"],"GALEX_sep_arcsec":uvpos.separation(spot(r)).arcsec} for r in rivals]
    rivals.append({"host":original,"host_ra":uvpos.ra.deg,"host_dec":uvpos.dec.deg,
                   "sep_high_arcsec":obspos.separation(uvpos).arcsec,"origin":"NED returned coordinates" if uv else "coordinates decoded from release-listed GALEXASC designation"})
    history={}
    for file in ("casebook.jsonl","historical_casebook.jsonl"):
        cases=rows(out/"casework"/file)
        match=next((c for c in cases if c.get("id")==sid),None)
        if match:
            history[file]={"aliases":match.get("alias"),"checks":match.get("checks"),"evidence":match.get("evidence")}
    return rivals,{"GALEX_to_optical":distances,"NED_cone":near,"Gaia_high_field":field.get("high_gaia"),"Gaia_status":field.get("high_gaia_status"),
                   "prior_archive_dossiers":history}


def literature(names):
    records=[]
    for name in names:
        params={"query.bibliographic":name,"rows":20,"select":"DOI,title,abstract,published"}
        key=re.sub(r"[^a-z0-9]","",name.lower())
        receipt={"source":"Crossref","query":name}
        try:
            r=requests.get("https://api.crossref.org/works",params=params,timeout=(10,25))
            receipt.update(url=r.url,http=r.status_code,sha256=hashlib.sha256(r.content).hexdigest())
            r.raise_for_status()
            works=r.json()["message"]["items"]
            hits=[]
            for item in works:
                words=" ".join(item.get("title",[]))+" "+item.get("abstract","")
                if key in re.sub(r"[^a-z0-9]","",words.lower()):
                    hits.append({"doi":item.get("DOI"),"title":item.get("title")})
            receipt.update(status="ok",results_examined=len(works),exact_names_in_titles_or_abstracts=hits)
        except Exception as exc:
            receipt.update(status="unavailable",error=str(exc)[:240])
        records.append(receipt)
    return records


def run_one(sid,src,previous,out,field,raw):
    obspos=SkyCoord(src.ra*u.deg,src.dec*u.deg)
    row=previous["instrument_checks"]["release_row"]
    galaxies=field.get("galaxy_hypotheses",[])
    rivals=[{"host":g["host"],"host_ra":g["host_ra"],"host_dec":g["host_dec"]} for g in galaxies]
    details={}
    if sid==TARGETS[0]:
        rivals,details=aliases_and_history(out,sid,galaxies,raw,field,obspos)
        names=[src["name"],row["IDENT"]]+[g["host"] for g in galaxies if "identified bright" in g["galaxy_support"]]
        details["named_paper_queries"]=literature(list(dict.fromkeys(names)))
    with tempfile.TemporaryDirectory(prefix="slew-pps-") as work:
        products,receipt=unpack(row["OBSID"],row["XIMNAME_B6"],Path(work))
        result={"detection_uid":sid,"name":src["name"],"native_release_row":row,"xmm_ra":src.ra,"xmm_dec":src.dec,
                "native_catalogue_soft_net_counts":src.xmm_counts,"native_catalogue_likelihood":src.xmm_ml,
                "archive":receipt,"local_rivals":rivals,"supporting_evidence":details,
                "high_state_field":{k:field.get(k) for k in ("high_field_radius_arcsec","high_gaia","high_gaia_status",
                    "other_eRASS1_sources","high_ero_status","same_OBSID_other_slew")}}
        for role,fn in (("image",lambda p:describe_image(p,obspos,rivals)),("events",lambda p:analyze_events(p,obspos))):
            if role not in products:
                result[role]={"status":"product unavailable in retrieved observation package"}
                continue
            try:
                result[role]=fn(products[role])
                result[role].setdefault("status","exact native soft image inspected")
            except Exception as exc:
                result[role]={"status":"retrieved product unmeasurable","error":str(exc)[:400]}
    return result


def main(out):
    past=json.loads((out/"review_independent_epochs.json").read_text())
    imgs={x["detection_uid"]:x for x in past["slew_images"]}
    fields={x["detection_uid"]:x for x in rows(out/"host_adjudication_fields.jsonl.gz")}
    raw={x["id"]:x for x in rows(out/"vetting/counterpart_evidence.jsonl.gz")}
    src=pd.read_csv(out/"catalog_sources.csv.gz",low_memory=False).set_index("detection_uid")
    with ThreadPoolExecutor(max_workers=2) as pool:
        jobs=[pool.submit(run_one,s,src.loc[s],imgs[s],out,fields[s],raw[s]) for s in TARGETS]
        results=[job.result() for job in jobs]
    record={"sources":results,"scope":"Exact release B6 image filenames are matched by tar-member basename. Retain even failed transports, missing images, or uncalibrated event astrometry. Crossref checks cover exact strings in up to 20 index results per name and cannot establish novelty."}
    (out/"slew_image_followup.json").write_text(json.dumps(plain(record),indent=2,allow_nan=False)+"\n")
    for r in results:
        print(r["name"],r["archive"]["status"],r["image"]["status"],r["events"]["status"],flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("out",type=Path)
    main(p.parse_args().out)
