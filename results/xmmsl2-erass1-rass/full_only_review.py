"""Inspect the individually strong, uncorroborated full-only slew detections."""
import argparse
import hashlib
import io
import json
import re
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from bs4 import BeautifulSoup

from local_hosts import measure
from prepare import truth

HEA="https://heasarc.gsfc.nasa.gov/FTP/xmm/data/"
PINNED=HEA+"catalogues/xmmsl2_total.fits.gz"
XSA="https://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio"


def digest(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f,"sha256").hexdigest()


def get(url,work,manifest,params=None,maxbytes=95_000_000):
    identity=json.dumps([url,params],sort_keys=True)
    path=work/hashlib.sha256(identity.encode()).hexdigest()
    if not path.exists():
        with requests.get(url,params=params,stream=True,timeout=(10,55)) as r:
            r.raise_for_status()
            if int(r.headers.get("Content-Length",0))>maxbytes:
                raise ValueError("Product exceeds the targeted single-product size budget")
            size=0
            tmp=path.with_suffix(".part")
            with tmp.open("wb") as stream:
                for block in r.iter_content(128*1024):
                    size+=len(block)
                    if size>maxbytes:
                        raise ValueError("Archive returned more than a targeted single product")
                    stream.write(block)
            tmp.replace(path)
    manifest.append(dict(url=url,params=params,bytes=path.stat().st_size,sha256=digest(path)))
    return path


def release_row(src,work,manifest,catalog_inputs):
    local=work.parent/"downloads/xmm_full.fits.gz"
    if local.is_file():
        manifest.append(dict(url=PINNED,bytes=local.stat().st_size,sha256=digest(local)))
    else:
        local=get(PINNED,work,manifest)
    frozen=next(v["sha256"] for v in catalog_inputs["resources"].values() if v.get("url")==PINNED)
    if manifest[-1]["sha256"]!=frozen:
        raise ValueError("Pinned release differs from the primary catalogue input")
    row=int(src.detection_uid.removeprefix("XMMSL2-total-"))-1
    with fits.open(local,memmap=False) as hdus:
        f=hdus[1].data[row]
        if str(f["UNIQUE_SRCNAME"]).strip()!=src["name"] or abs(float(f["RA"])-src.ra)>.0000003:
            raise ValueError("Exact row and coordinates did not agree with frozen release")
        return {k:str(f[k]).strip("\0 ") for k in ("OBSID","SOURCENUM","XIMNAME_B6","SRCNAME_B6",
                "VAL_FLAG","FLAG_COMMENT","MODE_ID","IDENT","ALTIDENT","ID_RESOURCE","ID_CATEGORY",
                "FILT_CNTS","UNFILT_CNTS","MJD_START","MJD_STOP")}


def native_flux_documentation(work,manifest):
    url="https://www.cosmos.esa.int/web/xmm-newton/xmmsl2-ug"
    try:
        path=get(url,work,manifest,maxbytes=8_000_000)
        doc=BeautifulSoup(path.read_bytes(),"html.parser")
        for t in doc.select("sup"):
            t.replace_with("^"+t.get_text(strip=True))
        for t in doc.select("sub"):
            t.replace_with("_"+t.get_text(strip=True))
        page=" ".join(doc.stripped_strings)
        hits=list(re.finditer(r"\b11\s*[.):-]?\s*Flux\s+conversion",page,re.I))
        sections=[]
        for hit in hits:
            text=page[hit.start():hit.start()+2800]
            stop=re.search(r"\b12\s*[.):-]?\s*(?:Differences|Changes)",text,re.I)
            if stop:
                text=text[:stop.start()]
            score=sum(bool(re.search(pattern,text,re.I)) for pattern in
                      (r"power.{0,2}law",r"photon|spectral index|slope",r"column|absorp",r"erg",r"6[. ]"))
            sections.append((score,text))
        if not sections or max(x[0] for x in sections)<2:
            raise ValueError("Could not verify substantive release-guide section 11")
        section=max(sections,key=lambda x:x[0])[1]
        gamma=re.search(r"(?:photon index|slope|spectral index|power[ -]?law index|gamma)\s*(?:(?:of|=|is|:)\s*)?([0-9]+(?:\.[0-9]+)?)",section,re.I)
        # Preserve the source's display order and superscripts so an ambiguous
        # coefficient, sign, or unit can be checked from the original text.
        tokens=section.split()
        nh=re.search(r"(?:(?:column(?: density)?|absorption)[^.;]{0,90}|N[_ ]?H.{0,30})",section,re.I)
        return {"url":url,"sha256":digest(path),"status":"ESA release guide section 11 located",
                "section_title":"Flux conversion","photon_index_if_explicit":float(gamma[1]) if gamma else None,
                "column_or_absorption_context":nh[0].strip() if nh else None,
                "section_excerpt":" ".join(tokens[:145]),"section_words":len(tokens),
                "ecf_unit_check":"The catalogue reader independently checks each release's median FLUX_B6/RATE_B6 after scaling."}
    except Exception as exc:
        return {"url":url,"status":"unverified","error":str(exc)[:380]}


def extract_image(path,filename):
    try:
        with fits.open(path) as hdus:
            hdus.verify("exception")
        return path
    except (OSError,ValueError):
        # ESA can package even a single requested product as a tar archive.
        with tarfile.open(path) as t:
            files=[m for m in t if m.isfile() and Path(m.name).name.lower()==filename.lower()]
            if len(files)!=1 or files[0].size>95_000_000:
                raise ValueError("Archive did not contain exactly one requested subimage")
            return io.BytesIO(t.extractfile(files[0]).read())


def image_stats(path,src):
    pos=SkyCoord(src.ra*u.deg,src.dec*u.deg)
    with fits.open(path,memmap=False) as hdus:
        h=next(h for h in hdus if h.header.get("NAXIS")==2 and WCS(h.header).has_celestial)
        grid=WCS(h.header).celestial
        xc,yc=grid.world_to_pixel(pos)
        scale=proj_plane_pixel_scales(grid)*3600
        if np.any(~np.isfinite(scale)) or min(scale)<=0 or not np.all((scale>.1)&(scale<30)):
            raise ValueError("Unexpected image pixel scale")
        im=np.asarray(h.data)
        radius=int(np.ceil(92/min(scale)))
        xx,yy=int(round(float(xc))),int(round(float(yc)))
        if not (radius<=xx<im.shape[1]-radius and radius<=yy<im.shape[0]-radius):
            raise ValueError("No complete local background annulus at the source position")
        y,x=np.mgrid[yy-radius:yy+radius+1,xx-radius:xx+radius+1]
        cnt=im[y,x].astype(float)
        if not np.all(np.isfinite(cnt)&(cnt>=0)&(abs(cnt-np.rint(cnt))<1e-5)):
            raise ValueError("Requested image is not an integer count image")
        sky=grid.pixel_to_world(x,y).icrs
        dx,dy=pos.spherical_offsets_to(sky)
        dx,dy=dx.arcsec,dy.arcsec
        sep=np.hypot(dx,dy)
        aperture,background=sep<=30.,(sep>=45.)&(sep<90.)
        n,m=int(cnt[aperture].sum()),int(cnt[background].sum())
        event_pixels=int(np.count_nonzero(cnt[aperture]))
        result={"aperture_arcsec":30.,"annulus_arcsec":[45.,90.],"gross_counts":n,
                "annulus_counts":m,"geometric_background_counts":float(m*aperture.sum()/background.sum()),
                "occupied_aperture_pixels":event_pixels,"pixel_scale_arcsec":list(scale),
                "brightest_pixel_counts":int(cnt[aperture].max()),"bunit":str(h.header.get("BUNIT",""))}
        if n:
            xy=np.stack([dx[aperture],dy[aperture]],axis=1)
            weights=cnt[aperture]
            center=np.average(xy,axis=0,weights=weights)
            moments=((xy-center).T*weights)@(xy-center)/weights.sum()
            eig=np.maximum(0.,np.linalg.eigvalsh(moments))
            result.update(centroid_offset_arcsec=float(np.linalg.norm(center)),rms_minor_major_arcsec=list(np.sqrt(eig)))
        result["interpretation"]="The event pixels and image use the same slew photons as the catalogue. Sparse moments and a geometric annulus do not test the calibrated, moving PSF or exposure-dependent background."
        return result


def inspect(src,work,manifest,catalog_inputs):
    attempts=[]
    try:
        f=release_row(src,work,manifest,catalog_inputs)
        filename=f["XIMNAME_B6"]
        obs=f["OBSID"]
        if not re.fullmatch(r"[A-Za-z0-9_.+-]+",filename) or "IMAGE_6" not in filename.upper() or not obs.isdigit():
            raise ValueError("Release lacks an exact native band-6 slew subimage identifier")
        routes=[(HEA+f"rev0/{obs}/PPS/{filename}",None),
                (XSA,{"obsno":obs,"level":"PPS","name":filename})]
        for url,params in routes:
            try:
                path=get(url,work,manifest,params)
                data=image_stats(extract_image(path,filename),src)
                return dict(status="official native soft-band image inspected",release_row=f,
                            product_url=url,parameters=params,image_measurements=data,failed_routes=attempts)
            except Exception as exc:
                attempts.append({"url":url,"params":params,"error":str(exc)[:320]})
        return dict(status="unresolved; official subimage unavailable or unsuitable",release_row=f,failed_routes=attempts)
    except Exception as exc:
        return dict(status="unresolved; native release reference unavailable",error=str(exc)[:320],failed_routes=attempts)


def plain(obj):
    # JSON null represents missing, including unavailable fluxes and absent IDs.
    if isinstance(obj,dict):
        return {str(k):plain(v) for k,v in obj.items()}
    if isinstance(obj,(tuple,list,np.ndarray)):
        return [plain(v) for v in obj]
    if pd.isna(obj):
        return None
    return obj.item() if isinstance(obj,np.generic) else obj


def flag(v):
    return str(v).lower() in ("true","1","1.0","t","yes")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",required=True,type=Path)
    p.add_argument("--work",required=True,type=Path)
    args=p.parse_args()
    out=args.out
    work=args.work/"full-only-products"
    work.mkdir(parents=True,exist_ok=True)
    cat=pd.read_csv(out/"catalog_sources.csv.gz",dtype={"xmm_obs":str},low_memory=False)
    chosen=cat.loc[~truth(cat.clean)&truth(cat.xmm_native_quality)&~truth(cat.xmm_full_only_corroborated)&
                   cat.xmm_counts.ge(5)&(cat.xmm_flux/cat.xmm_err).ge(3)].copy()
    pairs=pd.read_csv(out/"comparisons.csv.gz",low_memory=False)
    lm=pd.read_csv(out/"limits.csv",low_memory=False)
    inputs=json.loads((out/"catalog_inputs.json").read_text())
    summary=json.loads((out/"catalog_statistics.json").read_text())["full_only_audit"]
    expected=summary["strong_soft_photometry"]-summary["corroborated_rows"]
    if len(chosen)!=expected:
        raise ValueError(f"Full-only selection mismatch: expected {expected}, got {len(chosen)}")
    result=[]
    lines=[]
    manifest=[]
    for _,src in chosen.iterrows():
        sid=src.detection_uid
        comp=pairs.loc[pairs.detection_uid.eq(sid)].copy()
        comparable=[]
        for _,q in comp.iterrows():
            h,l=str(q.high),str(q.low)
            legit=all(m=="xmm" or not flag(src[m+"_ambiguous"]) for m in (h,l))
            legit &= h=="xmm" or flag(src[h+"_quality"])
            legit &= l=="xmm" or q.low_kind=="upper_limit" or flag(src[l+"_quality"])
            if (h,l) in [("ero","rass"),("rass","ero")] and q.low_kind=="detection":
                legit &= flag(src.ero_rass_consistent)
            vals=q[["high","low","low_kind","high_flux","high_error","low_flux_or_limit","low_error",
                    "ratio","guarded_ratio","model_floor_ratio"]].to_dict()
            vals["all_photometric_gates_except_slew_corroboration"]=bool(legit)
            comparable.append(vals)
        good=[q for q in comparable if q["all_photometric_gates_except_slew_corroboration"]]
        nominal=any(q["ratio"]>=10 for q in good)
        guarded=any(q["guarded_ratio"]>=10 for q in good)
        stress=any(q["model_floor_ratio"]>=10 for q in good)
        # Inspect even a nominal-only tenfold crossing before interpreting it.
        imaging=inspect(src,work,manifest,inputs) if nominal else {"status":"not triggered; no eligible nominal tenfold pair"}
        host,hits=measure(dict(detection_uid=sid,ra=src.ra,dec=src.dec,poserr_arcsec=src.xmm_poserr))
        host.pop("control_scores",None)
        host.update(legacy_alternatives=hits,release_identification=src.xmm_counterpart,
                    release_alternative=src.xmm_altident,release_class=src.xmm_xtyp,release_id_separation_arcsec=src.xmm_iddist)
        epochs={}
        for m in ("xmm","ero","rass"):
            epochs[m]={k:src.get(m+"_"+k) for k in ("flux","err","ml","time","mjd","mjd_min","mjd_max","sep","ambiguous","quality","poserr")}
            epochs[m]["name"]=src["name"] if m=="xmm" else src.get(m+"_name")
            if m!="xmm":
                limits=lm.loc[lm.detection_uid.eq(sid)&lm.mission.eq(m)]
                if len(limits)!=1:
                    raise ValueError(f"Missing or ambiguous map record for {sid} {m}")
                epochs[m]["limit"]=limits.iloc[0][["status","limit_flux_cgs","confidence","exp_s","method","map_ref"]].to_dict()
        warningcols=["xmm_issue","xmm_veropt","xmm_verext","xmm_verhalo","xmm_verbg","xmm_verfalse","xmm_veredge"]
        boundary=("Guarded and model-tested tenfold pair; instrument authenticity and external event/host vetting remain open." if stress else
                  "Nominal/guarded crossing fails the spectral stress test." if guarded else
                  "Nominal crossing alone; the uncertainty-guarded tenfold cut fails." if nominal else
                  "No admissible nominal tenfold comparison in available epochs/limits.")
        record=dict(detection_uid=sid,name=src["name"],ra=src.ra,dec=src.dec,obsid=src.xmm_obs,
                 source_number=src.xmm_id,net_soft_counts=src.xmm_counts,position_sigma_arcsec=src.xmm_poserr,
                 exposure_s=src.xmm_exp,snr=src.xmm_flux/src.xmm_err,extent_pixels=src.xmm_ext,
                 warning_flags={c:src.get(c) for c in warningcols},warning_fields=src.xmm_warning_fields,
                 release_validation_flag=src.get("val"),image_background_rate_cps=src.get("xmm_background_rate"),
                 native_B6_flux_per_cps=src.xmm_flux/src.xmm_rate if src.xmm_rate>0 else None,
                 full_only_reasons=src.xmm_full_only_reasons,other_strong_obsids=src.xmm_other_strong_slew_obsids,
                 close_secure_ero=src.xmm_close_secure_ero,close_secure_rass=src.xmm_close_secure_rass,
                 epochs=epochs,photometric_comparisons=comparable,optical=host,slew_product=imaging,
                 candidate=False,eligible_for_source_authenticity_followup=stress,boundary=boundary)
        result.append(record)
        lines.append(dict(detection_uid=sid,name=src["name"],ra=src.ra,dec=src.dec,obsid=src.xmm_obs,
                     epoch=src.xmm_time,net_soft_counts=src.xmm_counts,snr=src.xmm_flux/src.xmm_err,
                     xmm_flux=src.xmm_flux,ero_flux=src.ero_flux,rass_flux=src.rass_flux,
                     ero_limit=epochs["ero"]["limit"]["limit_flux_cgs"],rass_limit=epochs["rass"]["limit"]["limit_flux_cgs"],
                     max_nominal=max((q["ratio"] for q in good),default=None),
                     max_guarded=max((q["guarded_ratio"] for q in good),default=None),
                     max_model_floor=max((q["model_floor_ratio"] for q in good),default=None),
                     release_identification=src.xmm_counterpart,archive_image=imaging["status"],boundary=boundary))
    native_doc=native_flux_documentation(work,manifest)
    both=cat.xmm_flux.gt(0)&cat.xmm_rate.gt(0)
    native_doc["measured_median_B6_flux_CGS_per_cps"]=float((cat.loc[both,"xmm_flux"]/cat.loc[both,"xmm_rate"]).median())
    native_doc["positive_B6_flux_and_rate_rows"]=int(both.sum())
    dossier=dict(sources=result,external_inputs=manifest,native_xmm_release_calibration=native_doc,
                 selection="Full-only, native instrument & point-source quality, ≥5 net soft counts, flux/error≥3, absent independent corroboration.",
                 common_flux_units="observed 0.2–2.0 keV erg cm^-2 s^-1",
                 source_followup="Independent repeat or secure external positions are still required for promotion. No selected-epoch stacking.")
    (out/"full_only_review.json").write_text(json.dumps(plain(dossier),allow_nan=False,indent=2)+"\n")
    pd.DataFrame(lines).to_csv(out/"full_only_review.csv",index=False,float_format="%.10g")
    print(f"Audited {len(lines)} strong uncorroborated full-only slew rows",flush=True)


if __name__=="__main__":
    main()
