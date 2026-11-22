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
from astropy.io.votable import parse_single_table

from desi_lines import air_to_vac, desi_lines, fits_for_group
from desi_static import search as desi_search
from evidence_summary import adjudications, frozen_digest, ledger, stable_key
from focus import dump, fetch_spectrum
from repeat_blobs import clean_bridge, swapped_view
from vetting import (AGN_GAP, AGN_PASSED, Ancillary, Cache, atom, chord, fixed_external_pairs,
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


def archive_cones(points,prior):
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
            cones[s]=q
        h=prior[i].get("HSC",{})
        cones["HSC"]=h if h.get("status")=="measured" else safe(hsc_photometry,p)
        out.append(cones)
    return out


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
        m={k:atom(r[cols[k]]) for k in ("targetid","detectid","vi_z","vi_quality","dex_flag","source_type") if k in cols}
        hetsep=float(separation(r[cols["ra_hetdex"]],r[cols["dec_hetdex"]],p["ra"],p["dec"]))
        m.update(target_sep_arcsec=float(sep[j]),hetdex_sep_arcsec=hetsep)
        m["secure_visual_association"]=bool(hetsep<=1.25 and number(m.get("dex_flag"))==1 and
                                           number(m.get("vi_quality"))>=3 and number(m.get("vi_z"))>=0)
        if hetsep<=1.25 and number(m.get("dex_flag"))==1:
            m["diagnostics"]=safe(desi_lines,h,int(j),p["wave"],"blob")
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
            z=number(m.get("redshift"),None),secure=number(m.get("redshift_warning"))==0,
            spectype=str(m.get("spectype","")).strip().upper(),lines=(m.get("diagnostics") or {}).get("lines",{}),
            foreground=m.get("foreground_at_external_z"),supplementary_uv=m.get("supplementary_uv",{})))
    for m in vac.get("matches",[]):
        if number(m.get("dex_flag"))==1 and number(m.get("hetdex_sep_arcsec"),math.inf)<=1.25:
            hits.append(dict(source="HETDEX-DESI visual VAC",id=m.get("targetid"),sep_arcsec=m.get("target_sep_arcsec"),
                z=number(m.get("vi_z"),None),secure=m["secure_visual_association"],spectype=None,
                lines=(m.get("diagnostics") or {}).get("lines",{}),foreground=None,supplementary_uv={}))
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
            "z","rvz_redshift","refcode","rvz_bibcode","ls_id","type","flux_r","flux_ivar_r","sep_arcsec") if k in m}
            for m in q.get("matches",[])],**{k:q[k] for k in ("flux","flux_err","unit","layer") if k in q}}
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
                redshifts.append(dict(service=service,z=z,sep_arcsec=m["sep_arcsec"],
                                      reference=m.get("refcode",m.get("rvz_bibcode"))))
    return types,redshifts


def assess(row,group,line,path,first,aux,externals,receipts):
    p=row["peak"]; prior=row["vet"]; v=dict(prior)
    ps=[p,swapped_view(p)]
    for k in range(2):
        ps[k]=dict(ps[k],z_lya=air_to_vac(ps[k]["wave"])/REST["Ly_alpha"]-1)
    # sigma_A is the primary knot's fitted line width, never infer a partner width.
    ps[1].pop("sigma_A",None)
    source_knots=prior.get("confirmed_knot_catalogs") or []
    old=[{s:prior[s] for s in ("NED","SIMBAD","Legacy","HSC") if s in prior},
         dict(prior.get("secondary_archive_cones") or {},HSC=prior.get("secondary_HSC") or {})]
    cones=archive_cones(ps,old)
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
        agn,zlit=cone_flags(q,c)
        loc="primary" if k==0 else "secondary"
        if prop["multi_line_foreground"]:
            hard.append(loc+" knot has an independently measured foreground line and companion")
        if any(a["classification"]=="QSO" and a["same_line_and_lya_redshift"] for a in prop["external_active_or_star"]):
            hard.append(loc+" knot has a spectroscopic QSO counterpart at the Ly-alpha redshift with the same observed line")
        elif prop["external_active_or_star"] or agn:
            gaps.add(loc+" knot lies near a stellar or active counterpart; source-line association unresolved")
        if prop["competing_redshifts"] or zlit:
            gaps.add(loc+" knot has a competing counterpart redshift; association unresolved")
        for service in ("NED","SIMBAD","Legacy"):
            if c[service].get("status")=="queried":
                gaps.discard(service+(" search unavailable" if k==0 else " query unavailable at secondary knot"))
        if k==0 and (c["HSC"].get("status")=="measured" or c["Legacy"].get("status")=="queried"):
            gaps.discard("direct imaging checks unavailable")
        knots.append(dict(knot=loc,ra=q["ra"],dec=q["dec"],wave_air_A=q["wave"],trial_redshifts=expected(q),
             catalog_detections=[dict(m,flux_cgs=number(m.get("flux"),0)*1e-17 if m.get("flux") is not None else None,
                  flux_err_cgs=number(m.get("flux_err"),0)*1e-17 if m.get("flux_err") is not None else None) for m in detections],
             rest_EW_A=p.get("spectral_EW0_A") if k==0 else None,
             rest_EW_lower_A=p.get("spectral_EW0_lower_A") if k==0 else None,cone_results=concise_cones(c),
             published_redshifts=zlit,published_active_types=agn,HETDEX_DESI_VAC=vac,independent_spectroscopy=prop,
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
        archived_two_shot_bridge=bridge,knots=knots,LAN=lan,
        ordinary_AGN_screen="width and UV limits disfavor ordinary AGN; weak-line and obscured activity can remain"
        if agn_pass else "narrow width and/or UV constraints inadequate to exclude ordinary AGN")
    v.update(disposition=disposition,reasons=summary["reasons"],passed=summary["passed"],NED=cones[0]["NED"],
             SIMBAD=cones[0]["SIMBAD"],Legacy=cones[0]["Legacy"],HSC=cones[0]["HSC"],
             secondary_archive_cones=cones[1],secondary_spectroscopy=externals[1],targeted_spatial_source=summary)
    merged=dict(row,peak=p,vet=v,first_pass_line=first,key=stable_key(p),previous_ledger=path,
        previous_ledger_line=line,previous_disposition=prior["disposition"],previous_snapshot_sha256=frozen_digest(p,prior).hex())
    return merged,summary,dict(archive_cones=cones,spectra=receipts)


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
             "foreground":props["multi_line_foreground"],"high_z_evidence":props["independent_lya_evidence"],
             "external_line_flux_or_upper3_cgs":[{name:{v:h.get(v) for v in ("coverage","flux","err","upper3sigma")}
                 for name,h in m["lines"].items() if name in ("primary","heii","alt_oii","alt_hbeta","alt_oiii5007")} for m in props["spectra"]],
             "NED_SIMBAD_zs":k["published_redshifts"],
             "cone_status":{n:[c["status"],len(c["nearby"])] for n,c in k["cone_results"].items()},
             "external_spectra_status":k["external_cone_status"],"HETDEX_DESI_VAC_status":k["HETDEX_DESI_VAC"]["status"]}
        result.update({f"{role}_{n}":v for n,v in vals.items()})
    result.update(knot_separation_arcsec=s["archived_two_shot_bridge"][0]["knot_separation_arcsec"],
        shots=[m["shotid"] for m in s["archived_two_shot_bridge"]],
        bridge_empirical_snr=[[m["primary_bridge"]["sn"],m["secondary_bridge"]["sn"]] for m in s["archived_two_shot_bridge"]],
        bridge_flux_cgs=[[m["primary_bridge"]["value"]*1e-17,m["primary_bridge"]["scatter"]*1e-17] for m in s["archived_two_shot_bridge"]],
        LAN_status=s["LAN"]["status"],LAN_matches=s["LAN"].get("matches",[]),ordinary_AGN_screen=s["ordinary_AGN_screen"],
        prior_reasons=s["prior_reasons"],passed=s["passed"],remaining_reasons=s["reasons"])
    return {k:dump(v) if isinstance(v,(list,dict)) else v for k,v in result.items()}


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
    points=[];old=[]
    for _,_,_,row,_,_,_ in work:
        p=row["peak"];v=row["vet"]
        for q,oldext in ((p,v.get("SPARCL")),(swapped_view(p),v.get("secondary_spectroscopy"))):
            points.append(dict(q,z_lya=air_to_vac(q["wave"])/REST["Ly_alpha"]-1));old.append(oldext or {})
    external,master=direct_spectra(points,old,cache_dir) if points else ([],{})
    archive_cache=Cache(cache_dir)
    aux=Ancillary(archive_cache) if work else None
    summaries=[]; receipts=[]; overlays=[]
    try:
        for i,(first,line,path,row,g,edges,missing) in enumerate(work):
            original_snapshot=frozen_digest(row["peak"],row["vet"]).hex()
            item=dict(row,geometry_edge_lines=edges)
            updated,summary,receipt=assess(item,g,line,path,first,aux,external[2*i:2*i+2],{})
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
        repeat_geometry_systems_sha256=sha(rd/"systems.jsonl.gz"),confirmed_geometry=len(geometry),archives_rechecked=len(work),
        dispositions=dict(Counter(x["disposition"] for x in summaries)),products=["systems.csv","systems.json","vetted.jsonl.gz","receipts.json.gz"],
        scope="Every system with two distinct shots passing central and bridge thresholds when centered on both knots. Static DESI, targeted HETDEX-DESI spectroscopy, HETDEX LAN, NED, SIMBAD and optical counterparts are evaluated at the actual endpoint positions. An unclassified or unavailable external spectrum is not an astrophysical exclusion.")
    write_json(out/"summary.json",details)
    return details


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("release_dir",type=Path)
    ap.add_argument("--cache",type=Path,required=True)
    a=ap.parse_args()
    print(run(a.release_dir,a.cache))
