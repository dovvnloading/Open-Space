"""Archival identification of the two-shot HETDEX bridges that pass both knot tests."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import csv
import copy
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import re

import numpy as np
import requests
from astropy.io import fits
from astropy.io.votable import parse_single_table

from desi_lines import air_to_vac, desi_lines, fits_for_group, measure_arrays
from desi_static import search as desi_search
from evidence_summary import adjudications, frozen_digest, ledger, stable_key
from focus import dump
from repeat_blobs import clean_bridge, swapped_view
from vetting import (AGN_GAP, AGN_PASSED, Ancillary, Cache, atom, chord,
                     hsc_photometry, number, ordinary_agn_disfavored, separation, tap_batch)


REST={"Ly_alpha":1215.6701,"NV":1240.81,"CIV":1549.48,"HeII":1640.42,"CIII]":1908.73,
      "MgII":2799.94,"[OII]":3728.4835,"[NeIII]":3869.86,"Hdelta":4102.89,
      "Hgamma":4341.68,"Hbeta":4862.683,"[OIII]4959":4960.295,"[OIII]5007":5008.240,"Halpha":6564.614}
IS_ACTIVE={"QSO","AGN","SY1","SY2","SEYFERT","BLL","BLAZAR","BLA"}
CLASSICAL=("[OII]","[NeIII]","Hdelta","Hgamma","Hbeta","[OIII]4959","[OIII]5007","Halpha")
C=299792.458


def sha(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f,"sha256").hexdigest()


def write_json(path,obj):
    path.write_text(json.dumps(atom(obj),indent=2,allow_nan=False)+"\n")


def safe(fn,*args):
    try:
        return fn(*args)
    except Exception as exc:
        return {"status":"unavailable","error":str(exc)[:450]}


def expected(p):
    vac=air_to_vac(p["wave"])
    return {k:vac/l-1 for k,l in REST.items() if vac/l-1>=-.005}


def classical_test(w,f,e,p,z):
    """Test known foreground assignments only at independently anchored systemic z."""
    observed=air_to_vac(p["wave"])
    close=[k for k in CLASSICAL if z>=-.005 and abs(observed/((1+z)*REST[k])-1)*C<400]
    measurements={}
    for k in CLASSICAL if close else ():
        status,fits=fits_for_group(w,f,e,[REST[k]*(1+z)],(0,),(100,))
        if fits:
            _,_,v,cov,_=fits[0]
            se=float(np.sqrt(cov[0,0]))
            measurements[k]={"flux":float(v[0]),"err":se,"snr":float(v[0])/se,
                             "upper3":max(float(v[0]),0)+3*se}
        else:
            measurements[k]={"coverage":status}
    # The integrated [O II] primary is fitted independently as a resolved
    # doublet in desi_lines; this optional one-Gaussian test is conservative.
    supported=[k for k in close if number(measurements.get(k,{}).get("snr"))>=4 and
               any(other!=k and abs(REST[other]-REST[k])>8 and number(m.get("snr"))>=5
                   for other,m in measurements.items())]
    return {"redshift":z,"possible_primary":close,"fixed_lines":measurements,"supported_primary":supported}


def catalog_fallback(service,p):
    # Use a separate public coordinate interface when a TAP service fails.
    ra,dec=p["ra"],p["dec"]
    if service=="NED":
        url="https://ned.ipac.caltech.edu/cgi-bin/objsearch"
        params={"search_type":"Near Position Search","in_csys":"Equatorial","in_equinox":"J2000.0",
                "lon":f"{ra:.9f}d","lat":f"{dec:+.9f}d","radius":"0.05","of":"xml_main"}
    else:
        url="https://simbad.cds.unistra.fr/simbad/sim-coo"
        params={"Coord":f"{ra:.9f} {dec:+.9f}","CooFrame":"ICRS","CooEpoch":"2000",
                "Radius":"3","Radius.unit":"arcsec","output.format":"VOTable"}
    with requests.get(url,params=params,timeout=(6,22)) as response:
        response.raise_for_status()
        raw=response.content
        receipt={"query_url":response.url,"response_sha256":hashlib.sha256(raw).hexdigest()}
    if b'QUERY_STATUS' in raw and re.search(br'QUERY_STATUS["\'][^>]*value\s*=\s*["\']ERROR',raw,re.I):
        raise ValueError("Archive reported an unsuccessful cone search")
    table=parse_single_table(io.BytesIO(raw)).to_table(use_names_over_ids=True)
    names={re.sub("[^a-z0-9]","",k.lower()):k for k in table.colnames}
    def column(*choices):
        return next((names[k] for k in choices if k in names),None)
    rac,dc=column("radeg","ra","rad"),column("decdeg","dec","decd")
    if not rac or not dc:
        raise ValueError("Cone response lacks degree coordinate fields")
    matches=[]
    for row in table:
        r,d=number(row[rac]),number(row[dc])
        if not np.isfinite(r+d):
            raise ValueError("Cone coordinates are missing or require an unverified unit conversion")
        ds=float(separation(r,d,ra,dec))
        if ds<=3.001:
            m={"ra":r,"dec":d,"sep_arcsec":ds}
            for key,keys in {"name":("objectname","mainid"),"pretype":("type","otype"),
                             "z":("redshift","rvzredshift"),"refcode":("rvzbibcode","refcode")}.items():
                col=column(*keys)
                if col and not np.ma.is_masked(row[col]):
                    m[key]=atom(row[col])
            matches.append(m)
    return dict(service=service,status="queried",matches=matches,**receipt)


def archive_cones(points,prior,cache):
    from legacy_bricks import legacy_brick_cones
    with ThreadPoolExecutor(max_workers=3) as pool:
        queries=dict(zip(("NED","SIMBAD","Legacy"),pool.map(lambda s:tap_batch(s,points),("NED","SIMBAD","Legacy"))))
    out=[]
    for i,p in enumerate(points):
        cones={}
        for s,qs in queries.items():
            q=qs[i]
            if q.get("status")!="queried" and prior[i].get(s,{}).get("status")=="queried":
                q=dict(prior[i][s],refresh_failed=q)
            elif q.get("status")!="queried" and s in ("NED","SIMBAD"):
                q=dict(safe(catalog_fallback,s,p),failed_tap=q)
            earlier=prior[i].get(s) or {}
            if s in ("NED","SIMBAD") and q.get("status")=="queried" and earlier.get("status")=="queried":
                def source_key(m):
                    return (str(m.get("ra")),str(m.get("dec")),str(m.get("z",m.get("rvz_redshift"))),
                            str(m.get("main_id",m.get("prefname",m.get("name","")))))
                seen={source_key(m) for m in q.get("matches",[])}
                retained=[dict(m,trace="previous successful cone") for m in earlier.get("matches",[]) if source_key(m) not in seen]
                if retained:
                    q=dict(q,matches=list(q.get("matches",[]))+retained,
                           earlier_response_sha256=earlier.get("response_sha256"),earlier_query_url=earlier.get("query_url"))
            cones[s]=q
        h=prior[i].get("HSC",{})
        cones["HSC"]=h if h.get("status")=="measured" else safe(hsc_photometry,p)
        out.append(cones)
    # A successful empty TAP cone does not establish a usable optical exposure.
    retry=[i for i,c in enumerate(out) if queries["Legacy"][i].get("status")!="queried" or
           not c["Legacy"].get("matches") and c["HSC"].get("status")!="measured"]
    if retry:
        try:
            native=legacy_brick_cones([points[i] for i in retry],cache)
            if len(native)!=len(retry):
                raise ValueError("Native Legacy response multiplicity differs")
        except Exception as exc:
            native=[{"status":"unavailable","error":str(exc)} for _ in retry]
        for i,brick in zip(retry,native,strict=True):
            tap=out[i]["Legacy"]
            if brick.get("status") in ("queried","outside_footprint") and (brick.get("matches") or not tap.get("matches")):
                out[i]["Legacy"]=dict(brick,TAP_attempt=tap)
            else:
                out[i]["Legacy"]=dict(tap,native_attempt=brick)
    for p,c in zip(points,out,strict=True):
        c["optical_coverage"]=optical_coverage(c)
    return out


def optical_coverage(c):
    q=c["Legacy"]
    native=q.get("native_attempt") or q
    image=native.get("r_exposure") or {}
    if c["HSC"].get("status")=="measured":
        return dict(status="measured",basis="HSC central r aperture")
    if number(image.get("nexp_r"))>0:
        return dict(status="measured",basis="native Legacy r exposure at exact WCS pixel",nexp_r=image["nexp_r"])
    for m in q.get("matches",[]):
        if number(m.get("sep_arcsec"),math.inf)<=1.25 and number(m.get("nobs_r"))>0:
            return dict(status="measured",basis="Legacy Tractor counterpart with r exposures",nobs_r=m["nobs_r"],
                        sep_arcsec=m["sep_arcsec"])
    if native.get("status")=="outside_footprint" or number(image.get("nexp_r"))==0:
        return dict(status="outside_usable_r_coverage",basis="no published Legacy brick or zero on-position r exposures")
    return dict(status="unverified",basis="no calibrated aperture, exact on-position exposure or exposed r counterpart")


def wise_assessment(q):
    pool=q.get("matches") or []
    if len(pool)!=1:return None
    m=pool[0]
    if number(m.get("sep_arcsec"))>1.25 or number(m.get("maskbits"))!=0 or not all(
            0<=number(m.get("fracflux_"+b))<.1 for b in ("w1","w2")):
        return None
    f1,f2,i1,i2=(number(m.get(k)) for k in ("flux_w1","flux_w2","flux_ivar_w1","flux_ivar_w2"))
    if min(f1,f2,i1,i2)<=0:return None
    sn1,sn2=f1*math.sqrt(i1),f2*math.sqrt(i2)
    w1,w2=22.5-2.5*math.log10(f1)-2.699,22.5-2.5*math.log10(f2)-3.339
    ce=2.5/math.log(10)*math.hypot(1/sn1,1/sn2)
    return dict(w1_w2_vega=w1-w2,color_err=ce,w2_vega=w2,red_agn_flag=min(sn1,sn2)>=5 and w2<=15 and w1-w2-2*ce>=.8,
                ls_id=m.get("ls_id") or m.get("tractor_id"))


def lan_info(aux,p,ids):
    if "LAN" not in aux.tables:
        return {"status":"unavailable","error":aux.errors.get("LAN")}
    _,tab,ix,tree=aux.tables["LAN"]
    xyz=chord([p["ra"],p["secondary_ra"]],[p["dec"],p["secondary_dec"]])
    mid=xyz.sum(axis=0); mid/=np.linalg.norm(mid)
    probes=np.vstack((xyz,mid))
    rr=set()
    for vec in probes:
        rr.update(int(ix[j]) for j in tree.query_ball_point(vec,2*np.sin(np.deg2rad(40/3600)/2)))
    names={k.lower():k for k in tab.names}
    rr.update(map(int,np.flatnonzero(np.isin(tab[names["detectid"]],list(ids)))))
    # The published duplicate list may contain detections whose model centroid moved.
    duplicate=names.get("dups_detectid")
    if duplicate:
        for j,val in enumerate(tab[duplicate]):
            if ids.intersection(map(int,re.findall(r'\d{9,}',str(val)))):
                rr.add(j)
    matches=[]
    for j in sorted(rr):
        r=tab[j]
        rd=chord(r[names["ra"]],r[names["dec"]])
        ds=np.rad2deg(2*np.arcsin(np.minimum(1,np.linalg.norm(probes-rd,axis=1)/2)))*3600
        area=number(r[names["area_r_iso_circ"]],0)
        radius=math.sqrt(max(0.,area)/math.pi)
        cid=int(r[names["detectid"]]); dup=set(map(int,re.findall(r'\d{9,}',str(r[duplicate])))) if duplicate else set()
        z=number(r[names["z_hetdex"]],None)
        exact=bool(ids & ({cid}|dup))
        if exact or ds.min()<=max(5,radius):
            matches.append(dict(detectid=cid,exact_detection_or_duplicate=exact,z=z,
                z_agrees=z is not None and abs(z-p["z_lya"])<.015,flag_resolved=number(r[names["flag_resolved"]],None),
                knot1_knot2_midpoint_offsets_arcsec=ds,isophotal_radius_arcsec=radius))
    return {"status":"queried","matches":matches}


def vac_info(aux,p):
    if "DESI" not in aux.tables:
        return {"status":"unavailable","error":aux.errors.get("DESI")}
    h,t,_,_=aux.tables["DESI"]
    cols={x.lower():x for x in t.names}
    sep=separation(t[cols["target_ra"]],t[cols["target_dec"]],p["ra"],p["dec"])
    res=[]
    for j in np.flatnonzero(sep<=1.25):
        r=t[int(j)]
        m={k:atom(r[cols[k]]) for k in ("targetid","target_ra","target_dec","detectid","vi_z","vi_quality","dex_flag","source_type") if k in cols}
        hetsep=float(separation(r[cols["ra_hetdex"]],r[cols["dec_hetdex"]],p["ra"],p["dec"]))
        m.update(target_sep_arcsec=float(sep[j]),hetdex_sep_arcsec=hetsep)
        m["secure_visual_association"]=bool(hetsep<=1.25 and number(m.get("dex_flag"))==1 and
                                           number(m.get("vi_quality"))>=3 and number(m.get("vi_z"))>=0)
        if hetsep<=1.25 and number(m.get("dex_flag"))==1:
            m["diagnostics"]=safe(desi_lines,h,int(j),p["wave"],"blob")
            if m["secure_visual_association"]:
                try:
                    arrays=[]
                    for name in ("DESI_WAVE_VACUUM","DESI_SPEC","DESI_SPEC_ERR"):
                        arr=h[name].data
                        if arr.ndim==1 or arr.shape.count(1):
                            arrays.append(arr.ravel())
                        elif arr.shape[0]==len(t):
                            arrays.append(arr[j])
                        elif arr.shape[1]==len(t):
                            arrays.append(arr[:,j])
                        else:
                            raise ValueError("unidentified DESI spectral row axis")
                    m["foreground_at_external_z"]=classical_test(*arrays,p,number(m["vi_z"]))
                except Exception as exc:
                    m["additional_diagnostics_error"]=str(exc)
        res.append(m)
    return {"status":"queried","matches":res}


def direct_spectra(points,old,cache):
    from focus import trace_spectra
    coords=[(p["ra"],p["dec"]) for p in points]
    try:
        sky,receipts=desi_search(coords)
    except Exception as exc:
        sky={c:{"status":"unavailable","matches":[],"error":str(exc)} for c in coords}
        receipts={"error":str(exc)}
    todo=[dict(peak=p,first={},same_visit_two_line_repeats=[]) for p in points]
    mapping={c:{"static":dict(status=x["status"],matches=list(x["matches"]))} for c,x in sky.items()}
    spectra=cache/"focus-spectra"
    try:
        spec_receipts=trace_spectra(todo,mapping,spectra)
    except Exception as exc:
        spec_receipts={"error":str(exc)}
    result=[]
    for p,coord,t,prior in zip(points,coords,todo,old,strict=True):
        cone=sky[coord]
        new={m.get("spectrum_key") or m.get("sparcl_id") or str(m.get("specid")):copy.deepcopy(m)
             for m in prior.get("matches",[])}
        for m in t.get("external",{}).get("matches",[]):
            key=m.get("spectrum_key") or m.get("sparcl_id") or str(m.get("specid"))
            if key not in new or m.get("diagnostics"):
                new[key]=m
        for key,m in new.items():
            raw=spectra/(hashlib.sha256(key.encode()).hexdigest()+".npz")
            if raw.is_file() and number(m.get("redshift_warning"))==0 and number(m.get("redshift"))>=-.005:
                try:
                    with np.load(raw) as a:
                        m["foreground_at_external_z"]=classical_test(a["w"],a["f"],a["e"],p,number(m["redshift"]))
                        # N V supplements the existing forced UV diagnostics. Limits
                        # require full velocity/width coverage in the same external fibre.
                        zs=1+p["z_lya"]
                        _,fits=fits_for_group(a["w"],a["f"],a["e"],np.array([1238.821,1242.804])*zs,
                                              (-300,-150,0,150,300),(65,100,175))
                        if len(fits)==15:
                            m.setdefault("supplementary_uv",{})["nv_upper3"]=max(
                                max(0.,float(v.sum()))+3*math.sqrt(max(0.,float(cov.sum()))) for _,_,v,cov,_ in fits)
                except Exception as exc:
                    m["additional_diagnostics_error"]=str(exc)
        result.append(dict(status="matched" if new else cone["status"],matches=list(new.values()),
            initial_status=prior.get("status"),DESI_static={k:v for k,v in cone.items() if k!="matches"}))
    return result,dict(positions=coords,static=receipts,spectra=spec_receipts)


def properties(p,ext,vac):
    hits=[]
    for m in ext.get("matches",[]):
        if number(m.get("sep_arcsec"),math.inf)>1.25 or number(m.get("coadd_fiberstatus"),0)!=0 or number(
                m.get("mean_fibre_separation_arcsec"),0)>1.25:
            continue
        hits.append(dict(source=m.get("data_release"),id=m.get("specid"),sep_arcsec=m.get("sep_arcsec"),
            ra=number(m.get("ra"),None),dec=number(m.get("dec"),None),z=number(m.get("redshift"),None),secure=number(m.get("redshift_warning"))==0,
            spectype=str(m.get("spectype","")).strip().upper(),lines=(m.get("diagnostics") or {}).get("lines",{}),
            foreground=m.get("foreground_at_external_z"),supplementary_uv=m.get("supplementary_uv",{})))
    for m in vac.get("matches",[]):
        if number(m.get("dex_flag"))==1 and number(m.get("hetdex_sep_arcsec"),math.inf)<=1.25:
            hits.append(dict(source="HETDEX-DESI visual VAC",id=m.get("targetid"),sep_arcsec=m.get("target_sep_arcsec"),
                ra=number(m.get("target_ra"),None),dec=number(m.get("target_dec"),None),z=number(m.get("vi_z"),None),secure=m["secure_visual_association"],spectype=None,
                lines=(m.get("diagnostics") or {}).get("lines",{}),foreground=m.get("foreground_at_external_z"),supplementary_uv={}))
    yes=[]; competing=[]; foreground=[]; active=[]; ordinary=False
    narrow=(0<2.355*number(p.get("sigma_A"),math.inf)/p["wave"]*C<900)
    for h in hits:
        ls=h["lines"]; prim=ls.get("primary",{})
        primary=number(prim.get("snr"))>=4
        uvpair=number(ls.get("lya",{}).get("snr"))>=5 and number(ls.get("heii",{}).get("snr"))>=4
        z=h["z"]
        if h["secure"] and z is not None and abs(z-p["z_lya"])<.006:
            yes.append("secure independent redshift agrees")
        if uvpair:
            yes.append("fixed Ly-alpha and He II detected in one independent spectrum")
        if h["secure"] and z is not None and abs(z-p["z_lya"])>.02:
            competing.append(dict(source=h["source"],id=h["id"],z=z))
        if h["spectype"] in ("QSO","STAR") and h["secure"]:
            active.append(dict(source=h["source"],id=h["id"],classification=h["spectype"],z=z,
                same_line_and_lya_redshift=bool(primary and z is not None and abs(z-p["z_lya"])<.006)))
        if (h.get("foreground") or {}).get("supported_primary"):
            foreground.append({"source":h["source"],"id":h["id"],"z":z,"lines":h["foreground"]["supported_primary"],
                               "basis":"catalog spectrum redshift and an independently measured companion"})
        if number(ls.get("alt_oii",{}).get("snr"))>=4:
            companions=[k for k in ("alt_halpha","alt_hbeta","alt_oiii4959","alt_oiii5007","alt_neiii")
                        if number(ls.get(k,{}).get("snr"))>=5]
            if companions:
                foreground.append(dict(source=h["source"],id=h["id"],z=air_to_vac(p["wave"])/3728.4835-1,
                                       lines=["[OII]"]+companions,basis="fixed [O II] doublet fit and optical companion"))
        if primary and narrow and prim.get("flux",0)>0:
            nv=number(h["supplementary_uv"].get("nv_upper3"),math.inf)/prim["flux"]
            civ=number(ls.get("civ",{}).get("max_upper3sigma"),math.inf)/prim["flux"]
            galaxy=h["spectype"]=="GALAXY" and h["secure"] and z is not None and abs(z-p["z_lya"])<.006
            ordinary|=max(nv,civ)<.15 or (galaxy and min(nv,civ)<.15)
    return dict(spectra=hits,independent_lya_evidence=sorted(set(yes)),competing_redshifts=competing,
                multi_line_foreground=foreground,external_active_or_star=active,external_uv_agn_disfavored=ordinary)


def concise_cones(cones):
    return {s:{"status":q.get("status"),"nearby":[{k:m[k] for k in ("name","objname","main_id","pretype","otype",
            "ra","dec","prefname","z","rvz_redshift","zrefcode","z_refcode","refcode","rvz_bibcode","ls_id","type",
            "nobs_r","quality_r","flux_r","flux_ivar_r","sep_arcsec") if k in m}
            for m in q.get("matches",[])],**{k:q[k] for k in ("flux","flux_err","unit","layer","selected_release","outcome","r_exposure","basis","nexp_r","nobs_r") if k in q}}
            for s,q in cones.items()}


def cone_flags(p,cones):
    types=[]; redshifts=[]
    for service in ("NED","SIMBAD"):
        for m in cones[service].get("matches",[]):
            if number(m.get("sep_arcsec"),math.inf)>1.25:
                continue
            typ=str(m.get("pretype",m.get("otype",""))).strip().upper()
            if typ in IS_ACTIVE:
                types.append(dict(service=service,type=typ,sep_arcsec=m["sep_arcsec"]))
            z=number(m.get("z",m.get("rvz_redshift")))
            if z>=0 and abs(z-p["z_lya"])>.02:
                redshifts.append(literature_lines(p,service,m,z))
    return types,redshifts


def literature_lines(p,service,m,z):
    wave=p["wave"]
    pred=3728.4835*(1+z)
    predicted_air=pred*wave/air_to_vac(wave)
    ref=next((m[k] for k in ("zrefcode","z_refcode","rvz_bibcode","redshift_refcode","refcode") if m.get(k)),None)
    return dict(service=service,z=z,sep_arcsec=m.get("sep_arcsec"),reference=ref,
                published_object=m.get("prefname",m.get("main_id",m.get("name"))),
                predicted_OII_air_A=predicted_air,observed_minus_OII_A=wave-predicted_air,
                OII_offset_kms=C*(air_to_vac(wave)/pred-1),within_400_kms=abs(C*(air_to_vac(wave)/pred-1))<400)


def assess(row,group,line,path,first,aux,externals,cones,center_exams):
    p=row["peak"]; prior=row["vet"]; v=dict(prior)
    ps=[p,swapped_view(p)]
    for k in range(2):
        ps[k]=dict(ps[k],z_lya=air_to_vac(ps[k]["wave"])/REST["Ly_alpha"]-1)
    # sigma_A is the primary knot's fitted line width, never infer a partner width.
    ps[1].pop("sigma_A",None)
    source_knots=prior.get("confirmed_knot_catalogs") or []
    ids={int(k["detectid"]) for k in source_knots if k.get("detectid") is not None}
    for obs in group["observations"]:
        ids.update(int(obs["peak"][k]) for k in ("knot_detectid","secondary_detectid"))
    lan=safe(lan_info,aux,p,ids)
    if lan.get("status")!="queried":
        previous=[o.get("LAN_spatial_matches") for o in group["observations"]]
        if previous and all(isinstance(x,list) for x in previous) and not group.get("past_LAN_unavailable",True):
            lan=dict(status="queried",matches=[m for l in previous for m in l],source="frozen geometry-stage positional LAN tests",
                     refresh=lan)
    before=set(prior.get("reasons",[])); gaps=set(before); hard=[]; published=[]
    passed=list(prior.get("passed",[])); knots=[]; props=[]
    for k,(q,c,e) in enumerate(zip(ps,cones,externals,strict=True)):
        vac=safe(vac_info,aux,q)
        if vac.get("status")!="queried" and k==0 and prior.get("DESI"):
            m=dict(prior["DESI"],diagnostics=prior.get("DESI_lines") or {})
            m.setdefault("hetdex_sep_arcsec",m.get("sep_arcsec"))
            m["secure_visual_association"]=(number(m.get("dex_flag"))==1 and number(m.get("target_sep_arcsec"),math.inf)<=1.25 and number(m.get("vi_quality"))>=3)
            vac=dict(status="frozen_match",matches=[m],refresh=vac)
        prop=properties(q,e,vac)
        props.append(prop)
        detections=[source_knots[j] for j in (k,k+2) if j<len(source_knots)]
        catalog_alternatives=[]
        for m in detections:
            zcat=number(m.get("z_hetdex"))
            if zcat>=-.005 and abs(zcat-q["z_lya"])>.02:
                fits=[name for name in CLASSICAL if abs(air_to_vac(q["wave"])/(REST[name]*(1+zcat))-1)*C<400]
                if fits:
                    catalog_alternatives.append(dict(detectid=m.get("detectid"),candidate_ids=fits,z=zcat,
                        source_type=m.get("source_type"),z_method=m.get("z_hetdex_src"),confidence=m.get("z_hetdex_conf")))
        agn,zlit=cone_flags(q,c)
        loc="primary" if k==0 else "secondary"
        wise=wise_assessment(c["Legacy"])
        if wise and wise["red_agn_flag"]:
            gaps.add(loc+" knot has an isolated bright red WISE counterpart consistent with an active nucleus")
        if catalog_alternatives and not any("fixed Ly-alpha" in msg for msg in prop["independent_lya_evidence"]):
            gaps.add(loc+" knot wavelength agrees with a foreground transition at its HETDEX catalog redshift")
        if prop["multi_line_foreground"]:
            hard.append(loc+" knot has an independently measured foreground line and companion")
        if any(a["classification"]=="QSO" and a["same_line_and_lya_redshift"] for a in prop["external_active_or_star"]):
            hard.append(loc+" knot has a spectroscopic QSO counterpart at the Ly-alpha redshift with the same observed line")
        elif prop["external_active_or_star"] or agn:
            gaps.add(loc+" knot lies near a stellar or active counterpart; source-line association unresolved")
        if prop["competing_redshifts"] or zlit:
            gaps.add(loc+" knot has a competing counterpart redshift; association unresolved")
        for service in ("NED","SIMBAD","Legacy"):
            if c[service].get("status")=="queried" or service=="Legacy" and c[service].get("status")=="outside_footprint":
                gaps.discard(service+(" search unavailable" if k==0 else " query unavailable at secondary knot"))
        if c["optical_coverage"]["status"]=="measured":
            if k==0:gaps.discard("direct imaging checks unavailable")
            passed.append(loc+" knot has measured optical coverage")
        elif k==0:
            gaps.add("direct imaging checks unavailable")
        knots.append(dict(knot=loc,ra=q["ra"],dec=q["dec"],wave_air_A=q["wave"],trial_redshifts=expected(q),
             catalog_detections=[dict(m,flux_cgs=number(m.get("flux"),0)*1e-17 if m.get("flux") is not None else None,
                  flux_err_cgs=number(m.get("flux_err"),0)*1e-17 if m.get("flux_err") is not None else None) for m in detections],
             rest_EW_A=p.get("spectral_EW0_A") if k==0 else None,
             rest_EW_lower_A=p.get("spectral_EW0_lower_A") if k==0 else None,cone_results=concise_cones(c),
             published_redshifts=zlit,published_active_types=agn,WISE=wise,HETDEX_foreground_alternatives=catalog_alternatives,HETDEX_DESI_VAC=vac,independent_spectroscopy=prop,
             external_cone_status=e.get("status"),DESI_static=e.get("DESI_static")))
    if lan.get("status")=="queried":
        gaps.discard("published HETDEX LAN properties unavailable")
    else:
        gaps.add("published HETDEX LAN properties unavailable")
    for m in lan.get("matches",[]):
        if m.get("z_agrees") and m.get("flag_resolved")==1:
            published.append("previously resolved HETDEX LAN at matching redshift")
        elif m.get("z_agrees"):
            gaps.add("previous HETDEX LAN entry at matching redshift; novelty of extension requires review")
        elif m.get("exact_detection_or_duplicate"):
            gaps.add("published LAN association uses a discrepant redshift")
    if props[0]["independent_lya_evidence"] and not props[0]["competing_redshifts"]:
        gaps.discard("independent redshift unavailable")
        if not prior.get("low_z_families") or any("He II" in msg for msg in props[0]["independent_lya_evidence"]):
            gaps.discard("foreground interpretation unresolved")
        passed.extend(props[0]["independent_lya_evidence"])
    if props[1]["independent_lya_evidence"]:
        passed.extend("secondary knot: "+s for s in props[1]["independent_lya_evidence"])
    central_lya=[hit for check in center_exams for hit in check["spectral_tests"] if hit["same_bridge_lya"]]
    central_fg=[hit for check in center_exams for hit in check["spectral_tests"] if hit["same_bridge_foreground"]]
    if central_lya and not central_fg:
        gaps.discard("independent redshift unavailable")
        passed.append("independent spectrum inside both bridge masks detects the nominated line and supports Ly-alpha")
    if central_fg:
        hard.append("spectrum inside both bridge masks detects the nominated line and a corroborated foreground family")
    if any(hit["same_bridge_lya"] and hit["spectype"]=="QSO" for hit in central_lya):
        hard.append("spectroscopically identified QSO at the bridge wavelength and redshift")
    v["SPARCL"]=externals[0]
    agn_pass=ordinary_agn_disfavored(p,v,p["z_lya"]) or props[0]["external_uv_agn_disfavored"]
    if not hard and not published and (not gaps or gaps=={AGN_GAP}):
        if agn_pass:
            gaps.discard(AGN_GAP); passed.append(AGN_PASSED)
        else:
            gaps.add(AGN_GAP)
    if prior["disposition"]=="excluded":
        hard.extend(before or ["source excluded in earlier vetting"])
    if prior["disposition"]=="published":
        published.extend(before or ["source published in earlier vetting"])
    disposition="published" if published else "excluded" if hard else "unresolved" if gaps else "candidate"
    obs={o["first_pass_line"]:o for o in group["observations"]}
    edge=row.get("geometry_edge_lines") or [first,prior["matched_second_shot_original_line"]]
    if len(set(edge))!=2 or len({obs[j]["peak"]["shotid"] for j in edge})!=2 or not all(
            clean_bridge(obs[j][key]) for j in edge for key in ("cube","secondary_centered_cube")):
        raise ValueError("Source lacks two verified independent two-knot geometry measurements")
    bridge=[]
    for j in edge:
        o=obs[j]
        c=o["cube"]; c2=o["secondary_centered_cube"]
        bridge.append(dict(first_pass_line=j,shotid=o["peak"]["shotid"],cube=c.get("cube"),
            seeing_arcsec=c.get("seeing_arcsec"),knot_separation_arcsec=o["peak"].get("knot_separation_arcsec"),
            primary_line_window_air_A=c.get("line_window_air_A"),secondary_line_window_air_A=c2.get("line_window_air_A"),
            primary_central=c["central_aperture"],primary_bridge=c["bridge"],secondary_central=c2["central_aperture"],
            secondary_bridge=c2["bridge"],unit=c.get("aperture_unit")))
    summary=dict(release=p["release"],group=group["group"],first_pass_line=first,prior_disposition=prior["disposition"],
        disposition=disposition,prior_reasons=sorted(before),passed=list(dict.fromkeys(passed)),reasons=sorted(set(hard+published)|gaps),
        archived_two_shot_bridge=bridge,knots=knots,between_knots=center_exams,LAN=lan,
        ordinary_AGN_screen="width and UV limits disfavor ordinary AGN; weak-line and obscured activity can remain"
        if agn_pass else "narrow width and/or UV constraints inadequate to exclude ordinary AGN")
    v.update(disposition=disposition,reasons=summary["reasons"],passed=summary["passed"],NED=cones[0]["NED"],
             SIMBAD=cones[0]["SIMBAD"],Legacy=cones[0]["Legacy"],HSC=cones[0]["HSC"],
             secondary_archive_cones=cones[1],secondary_spectroscopy=externals[1],targeted_spatial_source=summary)
    merged=dict(row,peak=p,vet=v,first_pass_line=first,key=stable_key(p),previous_ledger=path,
        previous_ledger_line=line,previous_disposition=prior["disposition"],previous_snapshot_sha256=frozen_digest(p,prior).hex())
    return merged,summary,dict(archive_cones=cones)


def condensed(s):
    result={k:s[k] for k in ("release","group","first_pass_line","prior_disposition","disposition")}
    for role,k in zip(("A","B"),s["knots"],strict=True):
        meas=k["catalog_detections"]
        k1=meas[0] if meas else {}
        props=k["independent_spectroscopy"]
        vals={"ra_deg":k["ra"],"dec_deg":k["dec"],"air_A":k["wave_air_A"],"trial_z_Lya":k["trial_redshifts"]["Ly_alpha"],
             "trial_z_alternatives":{n:z for n,z in k["trial_redshifts"].items() if n!="Ly_alpha"},
             "flux_cgs":k1.get("flux_cgs"),"flux_err_cgs":k1.get("flux_err_cgs"),"rest_EW_A":k["rest_EW_A"],
             "rest_EW_lower_A":k["rest_EW_lower_A"],
             "HETDEX_classes_z":[[m.get("source_type"),m.get("z_hetdex"),m.get("z_hetdex_src")] for m in meas],
             "independent_zs":[[m["source"],m["id"],m["z"],m["secure"]] for m in props["spectra"]],
             "foreground":props["multi_line_foreground"],"HETDEX_foreground_alternatives":k["HETDEX_foreground_alternatives"],"high_z_evidence":props["independent_lya_evidence"],
             "external_line_flux_or_upper3_cgs":[{name:{v:h.get(v) for v in ("coverage","flux","err","upper3sigma")}
                 for name,h in m["lines"].items() if name in ("primary","heii","alt_oii","alt_hbeta","alt_oiii5007")} for m in props["spectra"]],
             "NED_SIMBAD_zs":k["published_redshifts"],
             "cone_status":{n:[c["status"],len(c["nearby"])] for n,c in k["cone_results"].items()},
             "external_spectra_status":k["external_cone_status"],"HETDEX_DESI_VAC_status":k["HETDEX_DESI_VAC"]["status"]}
        result.update({f"{role}_{n}":v for n,v in vals.items()})
    result.update(knot_separation_arcsec=s["archived_two_shot_bridge"][0]["knot_separation_arcsec"],
        shots=[m["shotid"] for m in s["archived_two_shot_bridge"]],
        bridge_empirical_snr=[[m["primary_bridge"]["sn"],m["secondary_bridge"]["sn"]] for m in s["archived_two_shot_bridge"]],
        between_knot_positions=[{k:c[k] for k in ("ra","dec","origins","published_counterparts","spectral_tests","external_status")} for c in s["between_knots"]],
        bridge_flux_cgs=[[m["primary_bridge"]["value"]*1e-17,m["primary_bridge"]["scatter"]*1e-17] for m in s["archived_two_shot_bridge"]],
        LAN_status=s["LAN"]["status"],LAN_matches=s["LAN"].get("matches",[]),ordinary_AGN_screen=s["ordinary_AGN_screen"],
        prior_reasons=s["prior_reasons"],passed=s["passed"],remaining_reasons=s["reasons"])
    return {k:dump(v) if isinstance(v,(list,dict)) else v for k,v in result.items()}


def bridge_dist(point,obs):
    p=obs["peak"]
    a,b=chord([p["ra"],p["secondary_ra"]],[p["dec"],p["secondary_dec"]])
    z=chord(*point)
    fac=206264.806
    direction=b-a; length=np.linalg.norm(direction)
    along=np.dot(z-a,direction)/length
    across=np.linalg.norm(z-a-along*direction/length)*fac
    d1,d2=np.linalg.norm(z-a)*fac,np.linalg.norm(z-b)*fac
    seeing=number(obs["cube"].get("seeing_arcsec"),None)
    return dict(shotid=p["shotid"],along_arcsec=along*fac,cross_arcsec=across,
        distances_arcsec=[d1,d2],seeing_arcsec=seeing,inside=seeing is not None and
        0<along<length and across<.75*seeing and min(d1,d2)>1.3*seeing)


def bridge_targets(p,cones,visits):
    a,b=chord([p["ra"],p["secondary_ra"]],[p["dec"],p["secondary_dec"]])
    mid=(a+b)/np.linalg.norm(a+b)
    midpoint=(float(np.rad2deg(np.arctan2(mid[1],mid[0]))%360),float(np.rad2deg(np.arcsin(mid[2]))))
    targets={midpoint:{"origins":["geometric midpoint"],"published_counterparts":[]}}
    for cone in cones:
        for service in ("NED","SIMBAD","Legacy"):
            for r in cone.get(service,{}).get("matches",[]):
                coord=(number(r.get("ra")),number(r.get("dec")))
                if not np.isfinite(coord).all():
                    continue
                spans=[bridge_dist(coord,o) for o in visits]
                if not any(0<m["along_arcsec"]<o["peak"]["knot_separation_arcsec"] and
                     m["cross_arcsec"]<=max(1.25,m["seeing_arcsec"] or 0) for m,o in zip(spans,visits,strict=True)):
                    continue
                loc=targets.setdefault(coord,{"origins":[],"published_counterparts":[]})
                if service not in loc["origins"]:loc["origins"].append(service)
                zz=number(r.get("z",r.get("rvz_redshift")))
                if zz>=0:
                    lit=literature_lines(p,service,dict(r,sep_arcsec=0.),zz)
                    if lit not in loc["published_counterparts"]:
                        loc["published_counterparts"].append(lit)
    return targets


def sdss_positions(points):
    """Native DR17 spectra at bridge positions; cone status and actual fibre fits are distinct."""
    if not points:return [],{}
    rects=[]
    for p in points:
        rad=1.25/3600
        dr=rad/max(.001,math.cos(math.radians(p["dec"])))
        rects.append(f"(ra between {p['ra']-dr:.9f} and {p['ra']+dr:.9f} and dec between {p['dec']-rad:.9f} and {p['dec']+rad:.9f})")
    q=("SELECT TOP 257 specObjID,ra,dec,z,zWarning,class,plate,mjd,fiberID,run2d,sciencePrimary FROM SpecObjAll WHERE "+
       " OR ".join(rects))
    trace={"query":q,"files":[]}
    out=[{"status":"unavailable","matches":[]} for _ in points]
    try:
        r=requests.get("https://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch",
                       params={"cmd":q,"format":"csv"},timeout=(5,18))
        r.raise_for_status()
        lines=[l for l in r.content.decode("utf-8-sig").splitlines() if l and not l.startswith("#")]
        table=csv.DictReader(lines)
        if not {"specobjid","ra","dec","zwarning","plate","fiberid"}<={x.lower() for x in table.fieldnames or []}:
            raise ValueError("DR17 returned no recognizable spectroscopic catalog schema")
        rows=[{k.lower():v for k,v in row.items()} for row in table]
        if len(rows)>=257:raise ValueError("SDSS result limit reached")
        trace.update(url=r.url,sha256=hashlib.sha256(r.content).hexdigest(),rows=len(rows))
    except Exception as exc:
        return out,dict(trace,error=str(exc)[:350])
    for x in out:x["status"]="queried"
    for row in rows:
        ra,dec=number(row.get("ra")),number(row.get("dec"))
        hits=[(j,float(separation(ra,dec,p["ra"],p["dec"]))) for j,p in enumerate(points)]
        hits=[(j,ds) for j,ds in hits if ds<=1.25]
        if not hits:continue
        run2d=row["run2d"].strip()
        if not re.fullmatch(r"[0-9A-Za-z_]+",run2d):continue
        plate,mjd,fiber=[int(row[k]) for k in ("plate","mjd","fiberid")]
        base=dict(ra=ra,dec=dec,redshift=row["z"],redshift_warning=row["zwarning"],spectype=row["class"],
                  specid=row["specobjid"],data_release="SDSS-DR17",plate=plate,mjd=mjd,fiberid=fiber,run2d=run2d)
        tracefile=dict(specid=row["specobjid"],attempts=[])
        arrays=None
        for survey in (("eboss","sdss") if not run2d.isdigit() else ("sdss","eboss")):
            url=f"https://data.sdss.org/sas/dr17/{survey}/spectro/redux/{run2d}/spectra/lite/{plate:04d}/spec-{plate:04d}-{mjd}-{fiber:04d}.fits"
            try:
                with requests.get(url,timeout=(5,18),stream=True) as rr:
                    rr.raise_for_status()
                    if int(rr.headers.get("Content-Length","0"))>16000000:raise ValueError("unexpected spectrum file size")
                    buf=bytearray()
                    for block in rr.iter_content(65536):
                        buf.extend(block)
                        if len(buf)>16000000:raise ValueError("spectrum exceeds single-target byte limit")
                with fits.open(io.BytesIO(buf),memmap=False) as h:
                    head=h[0].header
                    if (int(head["PLATEID"]),int(head["MJD"]),int(head["FIBERID"]))!=(plate,mjd,fiber):
                        raise ValueError("Spectrum plate/MJD/fibre identity mismatch")
                    spec=h["COADD"].data
                    w,f,iv=[np.asarray(spec[k],float) for k in ("loglam","flux","ivar")]
                    e=np.divide(1.,np.sqrt(np.maximum(iv,0)),out=np.full_like(iv,np.nan),where=(iv>0)&(spec["or_mask"]==0))
                    arrays=(10**w,f,e)
                tracefile["attempts"].append(dict(url=url,status="measured",sha256=hashlib.sha256(buf).hexdigest(),bytes=len(buf)))
                break
            except Exception as exc:
                tracefile["attempts"].append(dict(url=url,status="unavailable",error=str(exc)[:300]))
        trace["files"].append(tracefile)
        for j,ds in hits:
            m=dict(base,sep_arcsec=ds,status="spectrum_unavailable")
            if arrays is not None:
                m["diagnostics"]=safe(measure_arrays,*arrays,points[j]["wave"],"blob")
                m["status"]="measured" if m["diagnostics"].get("lines") else "unusable spectrum"
                if number(m["redshift_warning"])==0 and number(m["redshift"])>=-.005:
                    m["foreground_at_external_z"]=safe(classical_test,*arrays,points[j],number(m["redshift"]))
            out[j]["matches"].append(m)
    return out,trace


def center_tests(p,info,visits,ext,sd,aux):
    ext=dict(ext,matches=ext.get("matches",[])+sd.get("matches",[]))
    vac=safe(vac_info,aux,p)
    tests=[]
    for h in properties(p,ext,vac)["spectra"]:
        xy=h["ra"],h["dec"]
        masks=[bridge_dist(xy,o) for o in visits] if None not in xy else []
        inside=len(masks)==2 and all(g["inside"] for g in masks)
        ls=h["lines"]; main=ls.get("primary",{})
        line=number(main.get("snr"))>=4
        uv=number(ls.get("lya",{}).get("snr"))>=5 and number(ls.get("heii",{}).get("snr"))>=4
        zeq=h["secure"] and h["z"] is not None and abs(h["z"]-p["z_lya"])<.006
        fg=(h.get("foreground") or {}).get("supported_primary") or []
        doublet=number(ls.get("alt_oii",{}).get("snr"))>=4 and any(number(ls.get(k,{}).get("snr"))>=5 for k in
                              ("alt_halpha","alt_hbeta","alt_neiii","alt_oiii4959","alt_oiii5007"))
        tests.append(dict(source=h["source"],id=h["id"],spectype=h["spectype"],fiber_ra=h["ra"],fiber_dec=h["dec"],
            sep_arcsec=h["sep_arcsec"],systemic_z=h["z"],secure=h["secure"],nominated_line=main,
            same_bridge_lya=bool(inside and line and (zeq or uv)),same_bridge_foreground=bool(inside and line and (fg or doublet)),
            foreground_primary=fg,masks=masks,UV_pair=uv))
    return dict(ra=p["ra"],dec=p["dec"],origins=info["origins"],published_counterparts=info["published_counterparts"],
        query_radius_arcsec=1.25,air_wave_A=p["wave"],trial_z_Lya=p["z_lya"],requested_position_masks=[bridge_dist((p["ra"],p["dec"]),o) for o in visits],
        external_status=ext.get("status"),native_DESI=ext.get("DESI_static"),SDSS_status=sd.get("status"),HETDEX_DESI_VAC=vac,
        spectral_tests=tests,attribution="A redshift association needs an on-bridge fibre and a same-wavelength line detection; a central aperture can still include compact-galaxy emission.")


def run(d,cache_dir):
    d,cache_dir=map(Path,(d,cache_dir))
    out=d/"blob-followup"; out.mkdir(exist_ok=True)
    rd=d/"repeat-blobs"
    geom=json.loads((rd/"summary.json").read_text())
    original=sha(d/"vetted.jsonl.gz")
    if geom["input_sha256"]!=original:
        raise ValueError("Repeated-geometry systems belong to a different original catalog")
    geometry={}
    for _,sysrow in ledger(rd/"systems.jsonl.gz"):
        if not (sysrow["status"].startswith("source_vetting_") or sysrow["status"].startswith("unresolved;")):
            continue
        obs={x["first_pass_line"]:x for x in sysrow["observations"]}
        edges=[ed for ed in sysrow["geometry_edges"] if all(clean_bridge(obs[j].get("cube",{})) and
               clean_bridge(obs[j].get("secondary_centered_cube",{})) for j in ed["first_pass_lines"])]
        if not edges:
            continue
        sysrow["verified_edges"]=edges
        sysrow["past_LAN_unavailable"]=bool(geom.get("ancillary_catalog_errors",{}).get("LAN"))
        geometry[sysrow["group"]]=sysrow
    expected_count=geom["counts"].get("systems_with_both_knots_and_bridge_confirmed_twice",0)
    if len(geometry)!=expected_count:
        raise ValueError("Archived two-shot bridge counts disagree with per-system masks and statistics")
    registry={}
    for _,r in ledger(rd/"vetted.jsonl.gz"):
        g=r["repeat_blob_group"]
        if g in geometry:
            registry[r["first_pass_line"]]=(g,r["geometry_edge_lines"])
    known={g for g,_ in registry.values()}
    for g,obj in geometry.items():
        if g not in known:
            edge=obj["verified_edges"][0]["first_pass_lines"]
            registry[edge[0]]=g,edge
    work=[]
    if registry:
        for first,(line,row,path,_) in enumerate(adjudications(d,include_blob_followup=False),1):
            if first not in registry:
                continue
            g,edges=registry[first]
            # If spectral restoration failed in the frozen stage, an archival
            # redshift cannot certify the missing two-knot detector/AGN checks.
            missing=not row.get("vet",{}).get("confirmed_knot_catalogs")
            work.append((first,line,path,row,geometry[g],edges,missing))
    if len(work)!=len(geometry):
        raise ValueError("A spatially confirmed source cannot be traced to its original row")
    previous={}
    previous_hash=None
    if (out/"summary.json").exists() and (out/"receipts.json.gz").exists() and (out/"systems.json").exists():
        if json.loads((out/"summary.json").read_text())["input_sha256"]!=original:
            raise ValueError("Previous bridge follow-up belongs to a different frozen release")
        previous_hash=sha(out/"receipts.json.gz")
        with gzip.open(out/"receipts.json.gz","rt") as infile:
            archived=json.load(infile)["archive_checks"]
        old_sys=json.loads((out/"systems.json").read_text())
        for oldsystem,audit in zip(old_sys,archived,strict=True):
            previous[oldsystem["first_pass_line"]]=audit["archive_cones"]
    arc_points=[]; priors=[]; group_visits=[]
    for first,_,_,row,g,edges,_ in work:
        p,v=row["peak"],row["vet"]
        obs={o["first_pass_line"]:o for o in g["observations"]}
        visits=[obs[n] for n in edges];group_visits.append(visits)
        midpoint=next(iter(bridge_targets(p,[],visits)))
        mid=dict(p,ra=midpoint[0],dec=midpoint[1],wave=(p["wave"]+p["secondary_wave"])/2)
        old=previous.get(first) or [{s:v[s] for s in ("HSC","Legacy","NED","SIMBAD") if s in v},
             dict(v.get("secondary_archive_cones") or {},HSC=v.get("secondary_HSC") or {})]
        arc_points.extend((p,swapped_view(p),mid));priors.extend((old[0],old[1],old[2] if len(old)>2 else {}))
    arc=archive_cones(arc_points,priors,cache_dir) if work else []
    points=[];old=[];spans=[];center_points=[];center_meta=[]
    for i,(_,_,_,row,g,edges,_) in enumerate(work):
        r=3*i;start=len(points);center_start=len(center_points)
        p,v=row["peak"],row["vet"]
        for q,prior_ext in zip(arc_points[r:r+2],(v.get("SPARCL"),v.get("secondary_spectroscopy")),strict=True):
            points.append(dict(q,z_lya=air_to_vac(q["wave"])/REST["Ly_alpha"]-1));old.append(prior_ext or {})
        targets=bridge_targets(p,arc[r:r+3]+priors[r:r+3],group_visits[i])
        for (ra,dec),info in targets.items():
            q=dict(arc_points[r+2],ra=ra,dec=dec)
            q["z_lya"]=air_to_vac(q["wave"])/REST["Ly_alpha"]-1
            q.pop("sigma_A",None)
            points.append(q);old.append({});center_points.append(q);center_meta.append(info)
        spans.append((start,center_start,len(targets)))
    external,master=direct_spectra(points,old,cache_dir) if points else ([],{})
    sd,master["native_SDSS_dr17"]=sdss_positions(center_points)
    archive_cache=Cache(cache_dir)
    aux=Ancillary(archive_cache) if work else None
    summaries=[]; receipts=[]; overlays=[]
    try:
        for i,(first,line,path,row,g,edges,missing) in enumerate(work):
            original_snapshot=frozen_digest(row["peak"],row["vet"]).hex()
            item=dict(row,geometry_edge_lines=edges)
            start,cstart,ntargets=spans[i]
            centers=[center_tests(center_points[cstart+j],center_meta[cstart+j],group_visits[i],
                                  external[start+2+j],sd[cstart+j],aux) for j in range(ntargets)]
            updated,summary,receipt=assess(item,g,line,path,first,aux,external[start:start+2],arc[3*i:3*i+2],centers)
            receipt["archive_cones"]=arc[3*i:3*i+3]
            updated["previous_snapshot_sha256"]=original_snapshot
            if missing and updated["vet"]["disposition"] in ("candidate","unresolved"):
                reason="original source-level spectral and detector checks were incomplete"
                summary["disposition"]=updated["vet"]["disposition"]="unresolved"
                summary["reasons"]=updated["vet"]["reasons"]=sorted(set(summary["reasons"]+[reason]))
            summaries.append(summary);receipts.append(receipt);overlays.append(updated)
            print(d.name,"bridge",g["group"],summary["disposition"],summary["reasons"],flush=True)
    finally:
        if aux is not None:aux.close()
    write_json(out/"systems.json",summaries)
    flat=[condensed(s) for s in summaries]
    with (out/"systems.csv").open("w",newline="") as f:
        cols=flat[0].keys() if flat else ("release","group","first_pass_line","prior_disposition","disposition",
                                          "A_ra_deg","A_dec_deg","B_ra_deg","B_dec_deg","remaining_reasons")
        w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows(flat)
    with gzip.open(out/"vetted.jsonl.gz","wt") as f:
        for r in overlays:f.write(dump(r)+"\n")
    with gzip.open(out/"receipts.json.gz","wt") as f:
        f.write(dump(dict(DESI_static_and_spectra=master,archive_checks=receipts,ancillary=archive_cache.inputs))+"\n")
    details=dict(input_sha256=original,repeat_geometry_summary_sha256=sha(rd/"summary.json"),
        repeat_geometry_systems_sha256=sha(rd/"systems.jsonl.gz"),previous_targeted_receipts_sha256=previous_hash,confirmed_geometry=len(geometry),archives_rechecked=len(work),
        dispositions=dict(Counter(x["disposition"] for x in summaries)),products=["systems.csv","systems.json","vetted.jsonl.gz","receipts.json.gz"],
        scope="Every system with two distinct shots passing central and bridge thresholds when centered on both knots. Static DESI, targeted HETDEX-DESI spectroscopy, HETDEX LAN, NED, SIMBAD and optical counterparts are evaluated at both endpoints. Native DESI and SDSS DR17 were also searched at the midpoint and at catalog positions between the knots. Only an independently detected line in a fibre centred within both shot masks can contribute a bridge identification. An unclassified or unavailable external spectrum is not an astrophysical exclusion.")
    write_json(out/"summary.json",details)
    return details


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("release_dir",type=Path)
    ap.add_argument("--cache",type=Path,required=True)
    a=ap.parse_args()
    print(run(a.release_dir,a.cache))
