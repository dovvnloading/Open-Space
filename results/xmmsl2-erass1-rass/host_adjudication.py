"""Independent host localization and high-state rival audit; never assign novelty."""
import argparse
import gzip
import json
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord

from catalogs import Archive, NASA
from full_only_review import plain
from vet import is_galaxy, stellar_astrometry

EXTENDED={"EXP","DEV","SER","REX","COMP"}

def truth(value):
    return str(value).lower().strip() in {"1","1.0","true","t","yes"}

def finite(value):
    try:
        x=float(value)
        return x if np.isfinite(x) else np.nan
    except (ValueError,TypeError):
        return np.nan

def position(r,d):
    r,d=finite(r),finite(d)
    return SkyCoord(r*u.deg,d*u.deg) if np.isfinite(r) and abs(d)<=90 else None

def xname(value):
    match=re.search(r"J[0-9]{6}(?:\.[0-9]+)?[+-][0-9]{6}(?:\.[0-9]+)?",str(value),re.I)
    return match[0].upper() if match else str(value).strip()

def at(item,r="ra",d="dec"):
    return position(item.get(r),item.get(d))

def sep(a,b):
    return float(a.separation(b).arcsec) if a is not None and b is not None else np.nan

def records(path):
    try:
        return pd.read_csv(path,dtype=str,keep_default_na=False,low_memory=False).to_dict("records")
    except (pd.errors.EmptyDataError,FileNotFoundError):
        return []

class FieldCache(Archive):
    def save(self):
        # Inputs remain independent of the original catalogue manifest.
        pass

def fields(highs,work,out,kind):
    arc=FieldCache(work/kind,out)
    matches={sid:[] for sid,_,_,_ in highs}
    states={sid:"not queried" for sid in matches}
    try:
        if kind=="Gaia":
            url="https://gea.esac.esa.int/tap-server/tap/sync"
            table="gaiadr3.gaia_source"
            columns="source_id,ra,dec,parallax,parallax_error,pmra,pmra_error,pmdec,pmdec_error,ruwe,phot_variable_flag,phot_g_mean_mag"
            extra=(" AND (phot_variable_flag='VARIABLE' OR parallax>5*parallax_error OR "
                   "ABS(pmra)>5*pmra_error OR ABS(pmdec)>5*pmdec_error)")
            check=arc.tap(url,"SELECT TOP 1 source_id FROM "+table+" WHERE source_id=4472832130942575872")
            if not len(check) or str(check.iloc[0,0])!="4472832130942575872":
                raise ValueError("Gaia positive control not recovered")
        else:
            url,table=NASA,"erass1main"
            schema=arc.tap(url,"SELECT column_name,description FROM TAP_SCHEMA.columns WHERE table_name='erass1main'")
            lower={str(n).lower().replace("_",""):n for n in schema.iloc[:,0]}
            def field(aliases,required=False):
                for f in aliases.split():
                    if f in lower:
                        return lower[f]
                if required:
                    raise ValueError("unrecognized eROSITA archive fields "+str(list(lower))[:260])
            en=field("iauname erass1name name",True)
            er=field("ra",True)
            ed=field("dec",True)
            ml=field("detlike0 detectionlikelihood detlikelihood broadbandlikelihood likelihood")
            columns=",".join([f"{en} AS source_name",f"{er} AS ra",f"{ed} AS dec"]+([f"{ml} AS detlike"] if ml else []))
            extra=""
            test=arc.tap(url,"SELECT "+columns+" FROM "+table+" WHERE "
                         "1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',7.405249,-31.036103,0.001))")
            if not any("002937" in str(v) for v in test.iloc[:,0]):
                raise ValueError("Could not validate eRASS1 published source J002937.2-310209")
        failed=0
        for start in range(0,len(highs),12):
            batch=highs[start:start+12]
            if failed>=2:
                for sid,_,_,_ in batch:
                    states[sid]="unqueried after repeated archive failures"
                continue
            clauses=[]
            for sid,p,gate,year in batch:
                reach=gate+(90. if kind=="Gaia" else 0.)
                clauses.append(f"1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{p.ra.deg:.8f},{p.dec.deg:.8f},{reach/3600:.8f}))")
            try:
                data=arc.tap(url,"SELECT "+columns+" FROM "+table+" WHERE ("+" OR ".join(clauses)+")"+extra,maxrec=100000)
                if len(data)>=100000:
                    raise ValueError("archive row ceiling")
                rows=data.to_dict("records")
                sky=SkyCoord(pd.to_numeric(data["ra"]).to_numpy()*u.deg,
                             pd.to_numeric(data["dec"]).to_numpy()*u.deg)
                for sid,p,gate,year in batch:
                    reach=gate+(90. if kind=="Gaia" else 0.)
                    nearby=np.flatnonzero(p.separation(sky).arcsec<=reach)
                    for j in nearby:
                        row=rows[j]
                        loc=sky[j]
                        if kind=="Gaia":
                            # Gaia DR3 uses epoch 2016.0. Query an extra 90 arcsec
                            # before transporting the two tangent-plane components.
                            g=dict(plx=row.get("parallax"),e_plx=row.get("parallax_error"),pmra=row.get("pmra"),
                                   e_pmra=row.get("pmra_error"),pmdec=row.get("pmdec"),e_pmdec=row.get("pmdec_error"),ruwe=row.get("ruwe"))
                            variable=str(row.get("phot_variable_flag","")).upper()=="VARIABLE"
                            star=stellar_astrometry(g)
                            if not (star or variable) or loc is None:
                                continue
                            east,north=finite(g["pmra"]),finite(g["pmdec"])
                            if np.isfinite(east+north):
                                loc=loc.spherical_offsets_by(east*(year-2016.)*u.mas,north*(year-2016.)*u.mas)
                            if sep(p,loc)<=gate:
                                matches[sid].append(dict(id=str(row["source_id"]),ra=loc.ra.deg,dec=loc.dec.deg,
                                    sep_high_arcsec=sep(p,loc),stellar_astrometry=star,variable=variable,**g))
                        elif sep(p,loc)<=gate:
                            matches[sid].append(dict(name=str(row["source_name"]),ra=loc.ra.deg,dec=loc.dec.deg,
                                                     sep_high_arcsec=sep(p,loc),detlike=row.get("detlike")))
                    states[sid]="ok"
                failed=0
            except Exception as exc:
                failed+=1
                for sid,_,_,_ in batch:
                    states[sid]=str(exc)[:250]
    except Exception as exc:
        states={k:str(exc)[:300] for k in states}
    return matches,states,arc.records

def obs_year(row,mission):
    if mission=="rass":
        return 1991.
    mjd=finite(row.get(mission+("_mjd" if mission=="ero" else "_mjd_min")))
    return 2000.+(mjd-51544.5)/365.25 if np.isfinite(mjd) else (2020. if mission=="ero" else 2008.)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--work",type=Path,required=True)
    a=p.parse_args()
    out=a.out
    sources=records(out/"review_register.csv.gz")
    if not sources:
        sources=records(out/"input.csv")
    with gzip.open(out/"vetting/counterpart_evidence.jsonl.gz","rt") as f:
        evidence={r["id"]:r for line in f if (r:=json.loads(line))}
    local=records(out/"published_hosts/local_galaxies.csv")
    pub=records(out/"published_hosts/LS10.csv.gz")
    all_xmm=pd.read_csv(out/"catalog_sources.csv.gz",usecols=["detection_uid","source_group","xmm_obs","ra","dec","xmm_native_quality"],dtype={"xmm_obs":str})
    centers={}
    highs=[]
    for r in sources:
        sid=r["detection_uid"]
        high=r.get("robust_high") or r["high"]
        hp=at(r,high+"_ra",high+"_dec")
        hs=finite(r.get(high+"_poserr"))
        if hp is None or not np.isfinite(hs):
            raise ValueError("Unavailable high-state position/error: "+sid)
        # Search one RASS beam for confusing soft X-ray sources. Keep the
        # high-state 95% localization separate from that wider context.
        gate=2.45*hs
        field=max(gate,90. if high=="rass" else 35.)
        centers[sid]=(high,hp,hs,gate,field)
        highs.append((sid,hp,field,obs_year(r,high)))
    with ThreadPoolExecutor(max_workers=2) as pool:
        fs={k:pool.submit(fields,highs,a.work/"host-audit",out,k) for k in ("Gaia","eRASS1")}
        catalog={k:f.result() for k,f in fs.items()}
    hypotheses=[]
    summaries=[]
    dossier=[]
    for r in sources:
        sid=r["detection_uid"]
        high,hp,hs,hgate,field=centers[sid]
        detail=evidence[sid]
        raw=detail["alternatives_within_35arcsec"]
        refined=at(r)
        sigma=finite(r["poserr_arcsec"])
        gate=finite(detail["gate_arcsec"])
        gal=[g for g in local if g["detection_uid"]==sid]
        nw=[g for g in pub if g.get("ero_name") and g["ero_name"]==r.get("ero_name")]
        # Use the published bestCTP flag; a competitor's p_i alone is not
        # evidence of a primary optical counterpart.
        nwbest=[g for g in nw if truth(g.get("pub_bestCTP"))]
        neds=[j for j in raw.get("NED",[]) if is_galaxy(str(j.get("type","")))]
        sim=[j for j in raw.get("SIMBAD",[]) if is_galaxy(str(j.get("otype","")))]
        localstars=catalog["Gaia"][0][sid]
        otherx=[s for s in catalog["eRASS1"][0][sid]
                if xname(s["name"])!=xname(r.get("ero_name",""))]
        gs=catalog["Gaia"][1][sid]
        es=catalog["eRASS1"][1][sid]
        if r.get("ero_ul_status")=="outside_dr1":
            es="outside public DR1 footprint at the slew position"
        if es=="ok":
            galactic=hp.galactic
            side=np.cos(galactic.b.radian)*np.sin(galactic.l.radian)
            if side>0:
                es="outside the Galactic-longitude half-sky of public DR1"
            elif np.degrees(np.arcsin(abs(side)))*3600<field:
                es="public DR1 footprint boundary intersects high-state search field"
        other_same_slew=[]
        if high=="xmm":
            same=all_xmm[(all_xmm.xmm_obs==r["xmm_obs"]) & (all_xmm.source_group!=r["source_group"]) & all_xmm.xmm_native_quality]
            if len(same):
                close=hp.separation(SkyCoord(same.ra.to_numpy()*u.deg,same.dec.to_numpy()*u.deg)).arcsec
                other_same_slew=[{"id":str(s),"sep_arcsec":float(d)} for s,d in zip(same.detection_uid,close) if d<=field]
        seeds=[]
        for obj in gal:
            seeds.append(("Legacy morphology",at(obj,"host_ra","host_dec"),obj["ls_id"]))
        for obj in nwbest:
            if obj.get("pub_LS10_TYPE","").strip().upper() in EXTENDED:
                seeds.append(("NWAY primary extended",at(obj,"pub_LS10_RA","pub_LS10_DEC"),str(obj.get("pub_LS10_FULLID",""))))
        for obj in neds+sim:
            seeds.append(("identified galaxy",at(obj),str(obj.get("preferred_name") or obj.get("main_id"))))
        used=[]
        proposed=[]
        for source,loc,gid in seeds:
            if loc is None or sep(refined,loc)>gate or any(sep(loc,q)<1. for q in used):
                continue
            used.append(loc)
            gmatch=sorted((j for j in gal if sep(loc,at(j,"host_ra","host_dec"))<=2.),key=lambda j:finite(j["pcc"]))
            pmatch=sorted((j for j in nwbest if sep(loc,at(j,"pub_LS10_RA","pub_LS10_DEC"))<=2.),
                          key=lambda j:-finite(j.get("pub_NWAY_p_i")))
            nmatch=[j for j in neds if sep(loc,at(j))<=2.]
            smatch=[j for j in sim if sep(loc,at(j))<=2.]
            wmatch=[j for j in raw.get("WISE",[]) if sep(loc,at(j))<=2.]
            g=gmatch[0] if gmatch else {}
            n=pmatch[0] if pmatch else {}
            pcc=finite(g.get("pcc"))
            pany,pi=finite(n.get("pub_NWAY_p_any")),finite(n.get("pub_NWAY_p_i"))
            optname=next((j["main_id"] for j in smatch),next((j["preferred_name"] for j in nmatch),gid))
            hpsep=sep(hp,loc)
            rivals=[z for z in localstars if z["sep_high_arcsec"]<=hgate]
            column=dict(detection_uid=sid,name=r["name"],host=optname,origin=source,host_ra=loc.ra.deg,host_dec=loc.dec.deg,
                 same_NED_Legacy_galaxy=bool(nmatch and gmatch),same_SIMBAD_Legacy_galaxy=bool(smatch and gmatch),
                 bestCTP_position_match=bool(pmatch),published_primary_count=len(nwbest),NWAY_p_any=pany,NWAY_p_i=pi,
                 matched_legacy_id=g.get("ls_id",""),legacy_type=g.get("type",""),local_pcc=pcc,
                 shifted_as_good=g.get("shifted_at_least_as_good"),n_shifted=g.get("shifted_total"),
                 n_competing_galaxies=g.get("number_of_galaxy_alternatives"),
                 WISE_within_2as=len(wmatch),refined_sep_arcsec=sep(loc,refined),
                 high=high,high_ra=hp.ra.deg,high_dec=hp.dec.deg,high_sigma_arcsec=hs,high_r95_approx_arcsec=hgate,
                 host_high_sep_arcsec=hpsep,host_inside_high95=hpsep<=hgate,NED_name=" | ".join(j["preferred_name"] for j in nmatch),
                 SIMBAD_name=" | ".join(j["main_id"] for j in smatch),
                 Gaia_astrometric_or_variable_in_high95=len(rivals) if gs=="ok" else None,
                 Gaia_at_galaxy=sum(sep(loc,at(j))<=2. for j in rivals) if gs=="ok" else None,
                 other_eROSITA_main_in_high95=sum(j["sep_high_arcsec"]<=hgate for j in otherx) if es=="ok" else None,
                 other_eROSITA_main_in_field=len(otherx) if es=="ok" else None,field_radius_arcsec=field,same_OBSID_slew_rivals=len(other_same_slew),
                 original_stellar_evidence=r.get("stellar_evidence",""),original_associations=r.get("known_associations",""),
                 Gaia_high_field=gs,eRASS1_high_field=es)
            # Tags describe evidence for a galaxy localization. A nearby
            # classifiable star or other X-ray source remains visible regardless.
            independent=bool(gmatch and (nmatch or smatch))
            column["galaxy_support"]=("published primary + independently identified bright extended galaxy"
               if independent and pmatch and pi>=.95 and pany>=.8 and pcc<=.01 else
               "independently identified bright extended galaxy" if independent and pcc<=.01 else
               "published primary extended counterpart" if pmatch and pi>=.95 and pany>=.8 else
               "positional alternative")
            proposed.append(column)
            hypotheses.append(column)
        proposed.sort(key=lambda j:(not j["host_inside_high95"],not j["galaxy_support"].startswith("published primary +"),
                                   j["galaxy_support"]=="positional alternative",not np.isfinite(j["local_pcc"]),
                                   j["local_pcc"] if np.isfinite(j["local_pcc"]) else np.inf))
        top=proposed[0] if proposed else {}
        summaries.append(dict(detection_uid=sid,name=r["name"],high=high,high_r95_arcsec=hgate,
            original_disposition=r.get("disposition",""),original_stellar_evidence=r.get("stellar_evidence",""),
            original_associations=r.get("known_associations",""),original_novelty_checks=r.get("service_outcomes",""),
            Gaia_high_field=gs,eRASS1_high_field=es,number_galaxy_hypotheses=len(proposed),**{k:top.get(k) for k in ("host","host_ra","host_dec","galaxy_support",
            "host_inside_high95","host_high_sep_arcsec","same_NED_Legacy_galaxy","same_SIMBAD_Legacy_galaxy","bestCTP_position_match",
            "NWAY_p_any","NWAY_p_i","local_pcc","shifted_as_good","n_shifted","other_eROSITA_main_in_high95",
            "Gaia_astrometric_or_variable_in_high95","Gaia_at_galaxy")}))
        dossier.append(dict(detection_uid=sid,name=r["name"],high=high,high_field_radius_arcsec=field,
           high_gaia=catalog["Gaia"][0][sid],other_eRASS1_sources=otherx,same_OBSID_other_slew=other_same_slew,
           published_relations=nw,original_stellar_and_event_evidence={k:r.get(k) for k in ("stellar_evidence","known_associations","event_review")},
           high_gaia_status=gs,high_ero_status=es,galaxy_hypotheses=proposed))
    pd.DataFrame(summaries).to_csv(out/"host_adjudication.csv",index=False)
    pd.DataFrame(hypotheses).to_csv(out/"host_adjudication_alternatives.csv",index=False)
    with gzip.open(out/"host_adjudication_fields.jsonl.gz","wt") as f:
        for row in dossier:
            f.write(json.dumps(plain(row),allow_nan=False)+"\n")
    audit={"input_detections":len(sources),"hypotheses":len(hypotheses),
       "galaxy_evidence_counts":pd.Series([r.get("galaxy_support") or "no eligible galaxy hypothesis" for r in summaries]).value_counts().to_dict(),
       "strongest_cross_identifications":[r for r in summaries if str(r.get("galaxy_support","")).startswith("published primary +")],
       "coverage":{k:pd.Series(v[1]).value_counts().to_dict() for k,v in catalog.items()},
       "query_receipts":{k:v[2] for k,v in catalog.items()},
       "scope":"Each detection keeps all galaxy hypotheses inside its refined search circle. Galaxy records join at radius 2 arcsec. Galaxy support requires the stated independent catalogues, local z-band number density and, where available, the published primary NWAY identification with large absolute and relative probabilities. The different catalogues and match statistics can share data and are not multiplied into a posterior. High-state R95=2.45 times the catalogue matching scale is approximate, particularly for ROSAT. Gaia is propagated from J2016.0 with 90 arcsec of extra search margin. Gaia astrometry and variability identify possible alternatives, not proof that any star caused an event. Other eROSITA sources are later-time competitors for historical positions. A blank or failed coverage is not a veto or proof of absence. This audit does not adjudicate prior reporting or change any source classification."}
    (out/"host_adjudication_summary.json").write_text(json.dumps(plain(audit),indent=2,allow_nan=False)+"\n")
    print("Audited host and high-state alternatives in",len(sources),"slew detections",flush=True)


if __name__=="__main__":
    main()
