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
from astropy.time import Time
from scipy.spatial import ConvexHull, QhullError
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

from full_only_review import XSA, digest, plain

TARGETS = ("XMMSL2-total-25537", "XMMSL2-total-10038", "XMMSL2-total-21146")


def unpack(obs, image, work):
    params = {"obsno": obs, "level": "PPS", "name": image}
    receipt = {"url": XSA, "params": params, "selected_files": [], "members_seen": 0,
               "similar_images": [], "event_files_seen": [], "complete_listing": False}
    matched = {"images": [], "events": []}
    seen=set()
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
                    # The current pipeline's tile serials may differ from the
                    # catalogue even when the OBSID and instrument agree.
                    im=re.fullmatch(r"P"+obs+r"PN[A-Z0-9]{4}IMAGE_6([A-Z0-9]{3})\.(?:FTZ|FITS?)",canonical)
                    ev=re.fullmatch(r"P"+obs+r"PN[A-Z0-9]{4}(?:PIEVLI|EVENLI)([A-Z0-9]{4})\.(?:FTZ|FITS?)",canonical)
                    if not im and not ev:
                        continue
                    role="images" if im else "events"
                    if canonical in seen:
                        raise ValueError("Duplicated product basename "+canonical)
                    seen.add(canonical)
                    tile=(im or ev).group(1).lstrip("0") or "0"
                    dest=work/canonical
                    with dest.open("wb") as f:
                        shutil.copyfileobj(archive.extractfile(m),f,length=1024*1024)
                    if dest.stat().st_size!=m.size:
                        raise ValueError("Truncated member "+canonical)
                    item={"name":canonical,"tile":tile,"role":role,"bytes":m.size,"sha256":digest(dest),
                          "release_image_basename":canonical==image.upper(),"tar_member":m.name}
                    receipt["selected_files"].append(item)
                    matched[role].append((item,dest))
                receipt["complete_listing"]=True
        receipt["status"]="complete observation package enumerated"
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
        return SkyCoord(h.data["RA"][good],h.data["DEC"][good],unit="deg"),"RA,DEC degree columns"
    if not {"X","Y"}.issubset(names):
        raise ValueError("No RA/DEC or SAS sky X/Y event columns")
    col=[names.index(k)+1 for k in ("X","Y")]
    x,y=h.data["X"][good],h.data["Y"][good]
    try:
        pixel=WCS(h.header,keysel=["pixel"],colsel=col,relax=True)
        if pixel.pixel_n_dim!=2 or not pixel.has_celestial or list(pixel.wcs.colax)!=col:
            raise ValueError("Inapplicable pixel list WCS")
        world=pixel.all_pix2world(x,y,1)
        return SkyCoord(world[pixel.wcs.lng],world[pixel.wcs.lat],unit="deg"),"FITS pixel-list WCS on SAS X,Y"
    except (ValueError,TypeError,RuntimeError,AttributeError):
        pass
    w=WCS(naxis=2)
    for i,idx in enumerate(col):
        dim="XY"[i]
        for key,prefix,attr in (("TCRPX","CRPX","crpix"),("TCRVL","CRVL","crval"),("TCDLT","CDLT","cdelt"),("TCTYP","CTYP","ctype")):
            val=h.header.get(f"{key}{idx}",h.header.get(f"REF{dim}{prefix}"))
            if val is None:
                raise ValueError("Unavailable event list celestial WCS "+key+dim)
            getattr(w.wcs,attr)[i]=val
        if h.header.get(f"TCROT{idx}",0.)!=0 or any(re.fullmatch(rf"(?:TC|TP){idx}_\d+",t) for t in h.header):
            raise ValueError("Rotated event table WCS requires a separate calibration")
    if not w.has_celestial:
        raise ValueError("No recognized event-sky coordinate pair")
    world=w.all_pix2world(x,y,1)
    return SkyCoord(world[w.wcs.lng],world[w.wcs.lat],unit="deg"),"SAS sky X/Y table WCS (TCRPX, TCRVL, TCDLT, TCTYP or REFX/REFY)"


def time_span(header, observed=None):
    def val(k):
        return float(header[k]) if k in header else None
    scale=str(header.get("TIMESYS","")).strip().lower()
    ref=val("MJDREF")
    if ref is None and "MJDREFI" in header:
        ref=val("MJDREFI")+(val("MJDREFF") or 0)
    raw=(None,None) if observed is None else observed
    if observed is None:
        raw=val("TSTART"),val("TSTOP")
    if ref is not None and all(a is not None and np.isfinite(a) for a in raw):
        unit=u.Unit(header.get("TIMEUNIT","s"))
        factor=unit.to(u.day)
        mjd=np.asarray(raw)*factor+ref+float(header.get("TIMEZERO",0))*factor
        utc=Time(mjd,format="mjd",scale=scale).utc.mjd if scale in Time.SCALES else None
        return {"raw_scale":scale or "unspecified","native_mjd":list(mjd),"utc_mjd":list(utc) if utc is not None else None,
                "origin":"observed event TIME" if observed is not None else "TSTART/TSTOP"}
    if "DATE-OBS" in header and "DATE-END" in header:
        try:
            bounds=[Time(str(header[k]),scale="utc").mjd for k in ("DATE-OBS","DATE-END")]
            return {"raw_scale":"utc","utc_mjd":bounds,"native_mjd":bounds,"origin":"DATE-OBS/DATE-END header"}
        except ValueError:
            pass
    return {"raw_scale":scale,"origin":"no recognized time bounds"}


def time_match(span,row):
    a,b=float(row["MJD_START"]),float(row["MJD_STOP"])
    out={}
    for key in ("utc_mjd","native_mjd"):
        bounds=span.get(key)
        if bounds is None:
            continue
        lo,hi=bounds
        gap=max(0.,a-hi,lo-b)*86400
        overlap=max(0.,min(b,hi)-max(a,lo))*86400
        out[key]={"gap_s":gap,"overlap_s":overlap,"midpoint_offset_s":abs((lo+hi-a-b)/2)*86400}
    return out


def timing_rank(t):
    # Release MJD scale is not explicit; report tests on both clocks. Allow
    # 2 s of telemetry rounding. A time-only match does not imply sky coverage.
    good=[r for r in t.values() if r["gap_s"]<=2.]
    if not good:
        return None
    return min(good,key=lambda r:(-r["overlap_s"],r["midpoint_offset_s"]))


def image_footprint(path,pos,obs,row):
    with fits.open(path,memmap=False) as hdus:
        h=next(h for h in hdus if h.header.get("NAXIS")==2 and WCS(h.header).has_celestial)
        hdr=hdus[0].header.copy()
        hdr.update(h.header)
        grid=WCS(h.header).celestial
        sx,sy=h.header["NAXIS1"],h.header["NAXIS2"]
        xc,yc=map(float,grid.world_to_pixel(pos))
        if not np.all(np.isfinite([xc,yc])):
            raise ValueError("Nonfinite target pixel")
        cx,cy=(sx-1)/2,(sy-1)/2
        center=grid.pixel_to_world(cx,cy).icrs
        corners=grid.pixel_to_world(np.array([-.5,sx-.5,sx-.5,-.5]),np.array([-.5,-.5,sy-.5,sy-.5])).icrs
        inside=-.5<=xc<=sx-.5 and -.5<=yc<=sy-.5
        nearest=grid.pixel_to_world(np.clip(xc,-.5,sx-.5),np.clip(yc,-.5,sy-.5)).icrs
        window=time_span(hdr)
        detected=str(hdr.get("OBS_ID",hdr.get("OBSID",""))).strip()
        return {"obsid_header":detected,"obsid_compatible":not detected or detected==obs,
                "center_ra_dec":[center.ra.deg,center.dec.deg],"corners_ra_dec":list(zip(corners.ra.deg,corners.dec.deg)),
                "size_pixels":[sx,sy],"target_pixel":[xc,yc],"inside":inside,"distance_to_footprint_arcsec":pos.separation(nearest).arcsec,
                "center_distance_arcsec":pos.separation(center).arcsec,"window":window,"time_test":time_match(window,row),
                "edge_margin_pixels":min(xc+.5,yc+.5,sx-.5-xc,sy-.5-yc)}


def event_footprint(path,pos,obs,row):
    with fits.open(path,memmap=False) as hdus:
        h=next(h for h in hdus if h.name=="EVENTS")
        hdr=hdus[0].header.copy()
        hdr.update(h.header)
        times=np.asarray(h.data["TIME"],dtype=float)
        indices=np.flatnonzero(np.isfinite(times))
        if not len(indices):
            raise ValueError("Empty event-time array")
        sky,method=event_positions(h,indices)
        dx,dy=offsets(pos,sky)
        good=np.isfinite(dx)&np.isfinite(dy)
        if not good.any():
            raise ValueError("No finite calibrated event coordinates")
        xy=np.column_stack((dx[good],dy[good]))
        covered=False
        if len(xy)>=3:
            try:
                hull=ConvexHull(xy)
                covered=bool((hull.equations[:,-1]<=1e-5).all())
            except QhullError:
                pass
        detected=str(hdr.get("OBS_ID",hdr.get("OBSID",""))).strip()
        spans=time_span(hdr,[times[indices].min(),times[indices].max()])
        headspan=time_span(hdr)
        return {"obsid_header":detected,"obsid_compatible":not detected or detected==obs,
                "event_coord_calibration":method,"finite_events_all_energies":len(xy),"inside_all_event_convex_hull":covered,
                "minimum_observed_event_sep_arcsec":float(np.hypot(xy[:,0],xy[:,1]).min()),
                "east_north_min_arcsec":xy.min(axis=0).tolist(),"east_north_max_arcsec":xy.max(axis=0).tolist(),
                "event_time_window":spans,"header_time_window":headspan,"time_test":time_match(spans,row)}


def analyze_events(path,pos,release):
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
        sky,method=event_positions(h,good)
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
        obs_header=hdus[0].header.copy()
        obs_header.update(h.header)
        window={}
        if len(t):
            # Convert each selected photon with the same reference carried by
            # the product. The catalogue does not explicitly label its time
            # scale here, so publish the two conditional interval counts.
            for key in ("utc_mjd","native_mjd"):
                converted=[time_span(obs_header,[v,v]).get(key) for v in t]
                if all(v is not None for v in converted):
                    mjd=np.array([v[0] for v in converted])
                    mask=(mjd>=float(release["MJD_START"]))&(mjd<=float(release["MJD_STOP"]))
                    window[key]={"counts_in_release_interval":int(mask.sum()),"min_mjd":float(mjd.min()),"max_mjd":float(mjd.max())}
        return {"status":"native event-list positions and time order inspected","cuts":cuts,"pi_unit":pi_unit,"event_coordinate_calibration":method,
                "time_unit":str(h.columns["TIME"].unit),"counts_in_30as":len(t),"time_halves":temporal,"target_aperture_time_comparison":window,
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


def select_tiles(products,obspos,release):
    obs=release["OBSID"]
    trials={"images":[],"events":[]}
    paths={}
    for role,inspector in (("images",image_footprint),("events",event_footprint)):
        for product,path in products[role]:
            record={k:product[k] for k in ("name","tile","release_image_basename","sha256")}
            try:
                record.update(inspector(path,obspos,obs,release))
                record["timing_rank"]=timing_rank(record["time_test"])
                record["status"]="examined"
            except Exception as exc:
                record.update(status="unmeasurable",error=str(exc)[:280])
            paths[product["name"]]=path
            trials[role].append(record)
    inside=[v for v in trials["images"] if v.get("inside") and v.get("obsid_compatible")]
    timed=[v for v in inside if v.get("timing_rank") is not None]
    suitable=[v for v in trials["events"] if v.get("obsid_compatible") and v.get("timing_rank") is not None]
    def order(entry):
        rank=entry["timing_rank"]
        return (-rank["overlap_s"],rank["midpoint_offset_s"])
    # Pair one segment to one image. Association relies on the observation,
    # common segment suffix, source in the image, and event TIME overlap.
    # Match rank never depends on photons inside the source aperture.
    pairs=[(im,ev) for im in inside for ev in suitable if im["tile"]==ev["tile"] and
           (im in timed or not im.get("time_test"))]
    pairs.sort(key=lambda pair:(order(pair[1]),not pair[0]["release_image_basename"],pair[0]["center_distance_arcsec"],pair[0]["name"]))
    if pairs:
        image,event=pairs[0]
        rule="co-indexed slew segment: overlapping event time, containing image WCS, matched observation and segment"
    else:
        timed.sort(key=lambda im:(order(im),not im["release_image_basename"],im["center_distance_arcsec"],im["name"]))
        event_hull=[e for e in suitable if e.get("inside_all_event_convex_hull")]
        event_hull.sort(key=lambda e:(order(e),e["name"]))
        image=timed[0] if timed else None
        event=event_hull[0] if event_hull else None
        rule="independent single-image WCS/time and single-event-list footprint/time selections; no co-indexed segment"
    choice={"selection_rule":rule,"number_sky_containing_images":len(inside),
            "number_time_and_sky_images":len(timed),"number_time_compatible_events":len(suitable),
            "number_same_tile_spacetime_pairs":len(pairs),"release_mjd_interval":[release["MJD_START"],release["MJD_STOP"]],
            "closest_images_arcsec":[(v["name"],v["distance_to_footprint_arcsec"]) for v in
                sorted((v for v in trials["images"] if "distance_to_footprint_arcsec" in v),key=lambda v:v["distance_to_footprint_arcsec"])[:3]],
            "closest_observed_event_fields_arcsec":[(v["name"],v["minimum_observed_event_sep_arcsec"]) for v in
                sorted((v for v in trials["events"] if "minimum_observed_event_sep_arcsec" in v),key=lambda v:v["minimum_observed_event_sep_arcsec"])[:3]],
            "tried":trials,"coverage_note":"One image and at most one event list are measured. Overlapping scan products are never added. An image can include unexposed sky. Event convex hulls describe sampled extent only; an empty aperture is not a nondetection limit."}
    result={}
    for role,selected in (("image",image),("events",event)):
        choice["selected_"+role]=None if selected is None else selected["name"]
        if selected is not None:
            result[role]=paths[selected["name"]]
    return result,choice


def run_one(sid,src,previous,out,field,raw,baseline):
    obspos=SkyCoord(src.ra*u.deg,src.dec*u.deg)
    row=previous["instrument_checks"]["release_row"]
    galaxies=field.get("galaxy_hypotheses",[])
    rivals=[{"host":g["host"],"host_ra":g["host_ra"],"host_dec":g["host_dec"]} for g in galaxies]
    details={}
    if sid in baseline:
        rivals=baseline[sid]["local_rivals"]
        details=baseline[sid].get("supporting_evidence",{})
    elif sid==TARGETS[0]:
        rivals,details=aliases_and_history(out,sid,galaxies,raw,field,obspos)
        names=[src["name"],row["IDENT"]]+[g["host"] for g in galaxies if "identified bright" in g["galaxy_support"]]
        details["named_paper_queries"]=literature(list(dict.fromkeys(names)))
    with tempfile.TemporaryDirectory(prefix="slew-pps-") as work:
        downloaded,receipt=unpack(row["OBSID"],row["XIMNAME_B6"],Path(work))
        products,selection=select_tiles(downloaded,obspos,row)
        result={"detection_uid":sid,"name":src["name"],"native_release_row":row,"xmm_ra":src.ra,"xmm_dec":src.dec,
                "native_catalogue_soft_net_counts":src.xmm_counts,"native_catalogue_likelihood":src.xmm_ml,
                "archive":receipt,"product_selection":selection,"local_rivals":rivals,"supporting_evidence":details,
                "high_state_field":{k:field.get(k) for k in ("high_field_radius_arcsec","high_gaia","high_gaia_status",
                    "other_eRASS1_sources","high_ero_status","same_OBSID_other_slew")}}
        for role,fn in (("image",lambda p:describe_image(p,obspos,rivals)),("events",lambda p:analyze_events(p,obspos,row))):
            if role not in products:
                result[role]={"status":"no product verified for this sky position and epoch; no inference about source counts"}
                continue
            try:
                result[role]=fn(products[role])
                result[role].setdefault("status","native soft image with verified sky and epoch inspected")
            except Exception as exc:
                result[role]={"status":"retrieved product unmeasurable","error":str(exc)[:400]}
    return result


def main(out):
    past=json.loads((out/"review_independent_epochs.json").read_text())
    imgs={x["detection_uid"]:x for x in past["slew_images"]}
    fields={x["detection_uid"]:x for x in rows(out/"host_adjudication_fields.jsonl.gz")}
    raw={x["id"]:x for x in rows(out/"vetting/counterpart_evidence.jsonl.gz")}
    src=pd.read_csv(out/"catalog_sources.csv.gz",low_memory=False).set_index("detection_uid")
    saved=out/"slew_image_followup.json"
    baseline={s["detection_uid"]:s for s in json.loads(saved.read_text())["sources"]} if saved.is_file() else {}
    with ThreadPoolExecutor(max_workers=2) as pool:
        jobs=[pool.submit(run_one,s,src.loc[s],imgs[s],out,fields[s],raw[s],baseline) for s in TARGETS]
        results=[job.result() for job in jobs]
    record={"sources":results,"scope":"The native release-named B6 image and all alternate current pn soft tiles are compared by sky WCS, time and observation. Every archive member and tried sky/time footprint is recorded. Shape estimates reuse the high-state photons. Prior optical/literature evidence is reused on targeted imaging replay. Crossref checks cover exact strings in up to 20 index results per name and cannot establish novelty."}
    (out/"slew_image_followup.json").write_text(json.dumps(plain(record),indent=2,allow_nan=False)+"\n")
    for r in results:
        print(r["name"],r["archive"]["status"],r["image"]["status"],r["events"]["status"],flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("out",type=Path)
    main(p.parse_args().out)
