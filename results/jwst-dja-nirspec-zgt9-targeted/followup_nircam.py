#!/usr/bin/env python3
"""Bounded NIRCam catalogue crossmatch and exact-filter SED *diagnostic*.

Not a native image forced-photometry measurement. Never converts a missing
catalogue match into an upper limit, and never reports a physical z posterior.
Published JADES GOODS-S/N and GLASS A2744 catalogues only. Other fields are
explicitly unobserved by this test (not inferred non-detections).
"""
import argparse
from io import BytesIO
import json
import math
from pathlib import Path
import re
import sys
from urllib.error import HTTPError
from urllib.request import urlopen, Request

import numpy as np
from astropy.io import fits, votable
from astropy.coordinates import SkyCoord
from astropy import units as u

HERE=Path(__file__).resolve().parent
LEDGER=HERE/"literature_referee_2026-11-19.json"
URLS={
"goods-n":[
"https://archive.stsci.edu/hlsps/jades/jwst/nircam/goods-n/hlsp_jades_jwst_nircam_goods-n_photometry_v1.0_catalog.fits",
"https://archive.stsci.edu/pub/hlsp/jades/jwst/nircam/goods-n/hlsp_jades_jwst_nircam_goods-n_photometry_v1.0_catalog.fits"],
"goods-s":[
"https://archive.stsci.edu/hlsps/jades/jwst/nircam/goods-s-deep/hlsp_jades_jwst_nircam_goods-s-deep_photometry_v2.0_catalog.fits",
"https://archive.stsci.edu/pub/hlsp/jades/jwst/nircam/goods-s-deep/hlsp_jades_jwst_nircam_goods-s-deep_photometry_v2.0_catalog.fits"],
"abell2744":[
"https://archive.stsci.edu/hlsps/glass-jwst/jwst/nircam/abell2744/hlsp_glass-jwst_jwst_nircam_abell2744_multiband_v1.0_cat.fits",
"https://archive.stsci.edu/pub/hlsp/glass-jwst/jwst/nircam/abell2744/hlsp_glass-jwst_jwst_nircam_abell2744_multiband_v1.0_cat.fits"]}
BANDS=["F090W","F115W","F150W","F200W","F277W","F356W","F444W"]
CAT_CAP=760*1024*1024

def safe(v):
    if isinstance(v,np.generic):return safe(v.item())
    if isinstance(v,bytes):return v.decode("utf-8","replace")
    if isinstance(v,(float,int,str,bool)) or v is None:
        if isinstance(v,float) and not np.isfinite(v):return None
        return v
    if isinstance(v,list):return [safe(t) for t in v]
    if isinstance(v,dict):return {str(k):safe(t) for k,t in v.items()}
    if isinstance(v,np.ndarray):return safe(v.tolist())
    return str(v)

def dump(path,r):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(safe(r),indent=2,sort_keys=True,allow_nan=False)+"\n")

def get_file(urls,path):
    errors=[]
    if path.exists() and path.stat().st_size>2880:
        return {"url":"ON_DISK_FROM_THIS_RUN","bytes":path.stat().st_size}
    for url in urls:
        tmp=path.with_suffix(".fits.partial")
        n=0
        try:
            req=Request(url,headers={"User-Agent":"Open-Space-PR197-NIRCam-referee/1.0"})
            with urlopen(req,timeout=32) as response, tmp.open("wb") as out:
                if response.status!=200:raise ValueError("Unexpected HTTP "+str(response.status))
                while True:
                    buf=response.read(1<<20)
                    if not buf:break
                    n+=len(buf)
                    if n>CAT_CAP:raise ValueError("CAT_CAP_EXCEEDED")
                    out.write(buf)
            if n<2880:raise ValueError("TOO_SHORT")
            tmp.replace(path)
            return {"url":url,"bytes":n}
        except Exception as e:
            tmp.unlink(missing_ok=True)
            errors.append({"url":url,"reason":type(e).__name__+": "+str(e)[:180]})
    return {"status":"CATALOG_DOWNLOAD_FAILED", "attempts":errors}

def table_data(path):
    with fits.open(path,memmap=True) as h:
        for i,unit in enumerate(h):
            if not isinstance(unit,(fits.BinTableHDU,fits.TableHDU)):
                continue
            names=unit.columns.names
            norm={n.lower():n for n in names}
            ra=next((norm[k] for k in ("ra","ra_deg","alpha_j2000","ra_2000") if k in norm),None)
            dec=next((norm[k] for k in ("dec","dec_deg","delta_j2000","dec_2000") if k in norm),None)
            if ra is not None and dec is not None and len(unit.data)>0:
                return i,ra,dec,names
    raise ValueError("NO_KNOWN_RA_DEC_TABLE")

def extract_matches(path,targets):
    idx,racol,deccol,names=table_data(path)
    matches={}
    with fits.open(path,memmap=True) as h:
        t=h[idx].data
        ra=np.asarray(t[racol],float)
        dec=np.asarray(t[deccol],float)
        finite=np.isfinite(ra)&np.isfinite(dec)
        for g in targets:
            key=str(g["group_id"])
            center=SkyCoord(g["coord_deg"]["ra"]*u.deg,g["coord_deg"]["dec"]*u.deg)
            # Fast, exact SkyCoord angular separation only for points within 10 arcsec.
            r0,d0=g["coord_deg"]["ra"],g["coord_deg"]["dec"]
            around=np.flatnonzero(finite&(np.abs(dec-d0)<.003)&(np.abs((ra-r0)*np.cos(np.deg2rad(d0)))<.003))
            close=[]
            if len(around)>0:
                c=SkyCoord(ra[around]*u.deg,dec[around]*u.deg)
                order=np.argsort(center.separation(c).arcsec)[:4]
                for q in order:
                    row=int(around[q])
                    sep=float(center.separation(c[q]).arcsec)
                    close.append({"row":row,"sep_arcsec":sep,"ra":float(ra[row]),"dec":float(dec[row])})
            trusted=close[0] if close and close[0]["sep_arcsec"]<=.2 else None
            matches[key]={"status":"MATCH_0P2ARCSEC" if trusted else "NO_CLOSE_CATALOG_MATCH_NOT_UPPER_LIMIT",
                          "top4":close,
                          "source_row":trusted["row"] if trusted else None}
            if trusted:
                rec=t[trusted["row"]]
                meta={}
                for keycol in names:
                    low=keycol.lower()
                    if low in ("id","number","objid","source_id","photid","ra","dec","x","y") or "flag" in low and len(meta)<8:
                        try: meta[keycol]=safe(rec[keycol])
                        except Exception:pass
                matches[key]["identifiers"]=meta
                fluxes={}
                for band in BANDS:
                    namesb=[v for v in names if band.lower() in v.lower()]
                    # Prefer un-convolved aperture and a measured flux/uncertainty pair
                    # (the field nomenclature can vary by catalogue release).
                    chosen=None
                    for k in namesb:
                        lk=k.lower()
                        if "flux" in lk and ("circ1" in lk or "aper" in lk):
                            errors=[v for v in namesb if v.lower() in (lk+"_e",lk+"_err",lk+"_error",lk.replace("flux","fluxerr"),lk.replace("flux","err"))]
                            if errors: chosen=(k,errors[0]);break
                    if chosen is None:
                        for k in namesb:
                            lk=k.lower()
                            if "circ1" in lk and not any(w in lk for w in ("err","_e","error","snr","flag")):
                                errors=[v for v in namesb if v.lower() in (lk+"_e",lk+"_err",lk+"_error",lk.replace("circ1","circ1_err"))]
                                if errors:chosen=(k,errors[0]);break
                    if chosen:
                        f,e=chosen
                        flux=rec[f];err=rec[e]
                        fluxes[band]={"flux":safe(flux),"error":safe(err),"flux_column":f,"error_column":e,
                                     "flux_unit":str(h[idx].columns[f].unit or "UNDECLARED"),
                                     "error_unit":str(h[idx].columns[e].unit or "UNDECLARED")}
                matches[key]["measured_catalog_fluxes"]=fluxes
                matches[key]["available_band_columns"]={b:[n for n in names if b.lower() in n.lower()][:30] for b in BANDS}
    return matches,{"hdu":idx,"ra_column":racol,"dec_column":deccol,"table_columns":len(names)}

def throughput(band):
    url="https://svo2.cab.inta-csic.es/theory/fps/fps.php?ID=JWST/NIRCam."+band
    with urlopen(Request(url,headers={"User-Agent":"Open-Space-PR197-referee/1.0"}),timeout=30) as resp:
        xml=resp.read(4_000_000)
    tb=votable.parse_single_table(BytesIO(xml)).to_table()
    w=np.asarray(tb["Wavelength"],float)/10000.  # SVO FPS wavelengths in Angstrom; microns
    t=np.asarray(tb["Transmission"],float)
    ok=np.isfinite(w)&np.isfinite(t)&(w>.3)&(t>0)
    if ok.sum()<15:raise ValueError("invalid filter "+band)
    w,t=w[ok],t[ok];o=np.argsort(w)
    return w[o],t[o],{"url":url,"profile_bytes":len(xml),"points":len(w)}

def sed(w,kind,z,slope):
    w=np.asarray(w,float)
    base=(w/2)**slope
    if kind=="highz":
        return base*np.where(w<.121567*(1+z),.008,1.)
    if kind=="balmer":
        return base*np.exp(-.5*(2./w)**.7)*np.where(w<.4*(1+z),.10,1.)
    if kind=="line":
        out=base*.08
        for l,amp in ((.3727,.4),(.4861,.2),(.5007,1.),(.6563,.5)):
            cen=l*(1+z)
            out+=amp*np.exp(-.5*((w-cen)/.016)**2)
        return out
    if kind=="dwarf":
        out=base
        for c,d in ((1.4,.6),(1.9,.7),(2.7,.6),(3.3,.5)):
            out*=1-d*np.exp(-.5*((w-c)/.1)**2)
        return out
    raise ValueError(kind)

def fit_photo(matches,profiles):
    # Only physically interpretable photometry and throughput in the same
    # catalogue units is allowed. Merely having a catalogue match is not enough.
    data=[]
    for b,v in matches.get("measured_catalog_fluxes",{}).items():
        try:
            y,e=float(v["flux"]),float(v["error"])
            if not np.isfinite(y+e) or e<=0:continue
            if v["flux_unit"]=="UNDECLARED" or v["error_unit"]!=v["flux_unit"]:continue
            if b in profiles:data.append((b,y,e,v["flux_unit"]))
        except (KeyError,ValueError,TypeError):pass
    if len(data)<4 or len({u for _,_,_,u in data})!=1:
        return {"status":"INSUFFICIENT_VERIFIED_SAME_UNIT_FLUXES_AND_FILTERS",
                "bands":len(data),"not_a_physical_posterior":True}
    obs=np.array([y for b,y,e,u in data]);sigma=np.array([e for b,y,e,u in data])
    chis={}
    grid={"highz":np.arange(12.1,18.01,.3),
          "balmer":np.arange(1.5,5.5,.2),"line":np.arange(.6,8,.4),
          "dwarf":[0]}
    for fam,zvals in grid.items():
        best={"chi2":float("inf")}
        for z in zvals:
            for slope in [-1.,0.,1.]:
                f=np.array([
                    np.trapezoid(sed(profiles[b][0],fam if fam!="highz" else "highz",float(z),slope)*
                                 profiles[b][1]/profiles[b][0],profiles[b][0])/
                    np.trapezoid(profiles[b][1]/profiles[b][0],profiles[b][0])
                    for b,*_ in data])
                norm=np.sum(obs*f/sigma**2)/np.sum((f/sigma)**2)
                norm=max(0.,norm)
                chi=float(np.sum(((obs-norm*f)/sigma)**2))
                if chi<best["chi2"]:
                    best={"chi2":chi,"z_grid":float(z),"slope":slope,"amplitude":float(norm)}
        chis[fam]=best
    return {"status":"LIMITED_SVO_FILTER_CONVOLVED_SED_CHI2_NOT_POSTERIOR",
            "filters_used":[b for b,*_ in data],
            "photometry_unit":data[0][-1],"n_measured_fluxes":len(data),
            "best_fit_by_template_family":chis,"not_a_physical_posterior":True}

def run(out):
    ledger=json.loads(LEDGER.read_text())
    groups=[g for g in ledger["entries"] if g["referee_class"]=="UNRESOLVED_NOT_NEW"]
    where={}
    for g in groups:
        ra,dec=g["coord_deg"]["ra"],g["coord_deg"]["dec"]
        field=("abell2744" if ra<10 else
               "goods-n" if dec>30 else
               "goods-s" if 50<ra<55 else None)
        where[str(g["group_id"])]=field
    work=Path("/tmp/dja-pr197-photcatalog");work.mkdir(parents=True,exist_ok=True)
    report={"method":"Existing published NIRCam catalogues at exact sky positions; no forced photometry or uncertainty calibration",
            "status":"CATALOG_BOUNDED_PROBE_NO_REDSHIFT_POSTERIOR","groups":{},
            "catalogues":{},"throughput":{}}
    for g in groups:
        key=str(g["group_id"])
        report["groups"][key]={"coord_deg":g["coord_deg"],
                               "catalogue_field":where[key],
                               "status":"NOT_CHECKED_YET"}
    dump(out/"followup_nircam.json",report)
    for field in ["goods-n","goods-s","abell2744"]:
        targets=[g for g in groups if where[str(g["group_id"])]==field]
        if not targets:continue
        path=work/(field+".fits")
        dl=get_file(URLS[field],path)
        report["catalogues"][field]=dl
        if "status" in dl:
            for g in targets:report["groups"][str(g["group_id"])]["status"]="CATALOG_DOWNLOAD_FAILED_NO_UPPER_LIMIT"
            dump(out/"followup_nircam.json",report);continue
        try:
            matches,meta=extract_matches(path,targets)
            report["catalogues"][field]["parsed_table"]=meta
            for g in targets:
                key=str(g["group_id"])
                report["groups"][key].update(matches[key])
                report["groups"][key]["flux_at_slit_position"]=bool(matches[key].get("measured_catalog_fluxes"))
        except Exception as e:
            report["catalogues"][field]["parse_error"]=type(e).__name__+": "+str(e)[:250]
            for g in targets:report["groups"][str(g["group_id"])]["status"]="CATALOG_PARSE_FAILED_NO_UPPER_LIMIT"
        dump(out/"followup_nircam.json",report)
    for g in groups:
        key=str(g["group_id"])
        if where[key] is None:report["groups"][key]["status"]="NO_COVERED_CATALOGUE_IN_THIS_PROBE_NOT_NONDETECTION"
    if any(len(g.get("measured_catalog_fluxes",{}))>=4 for g in report["groups"].values()):
        responses={}
        for band in BANDS:
            try: responses[band]=throughput(band);report["throughput"][band]=responses[band][2]
            except Exception as e:
                report["throughput"][band]={"status":"UNAVAILABLE","error":type(e).__name__+": "+str(e)[:180]}
        for key,item in report["groups"].items():
            if "measured_catalog_fluxes" in item:
                item["convolved_template_comparison"]=fit_photo(item,responses)
    dump(out/"followup_nircam.json",report)
    print(json.dumps({"NIRCam_status":report["status"],
                      "catalogues":{k:v.get("status","received") for k,v in report["catalogues"].items()},
                      "groups":{k:v["status"] for k,v in report["groups"].items()}}))

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    run(args.out)
