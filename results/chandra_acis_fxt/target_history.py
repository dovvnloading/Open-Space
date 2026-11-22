"""All-years pointed ACIS imaging and CSC 2.2 follow-up at supplied fixed positions.

python results/chandra_acis_fxt/target_history.py --leads selected_leads.json --out history

Only target characterization is performed. These observations never extend the blind
survey's exposure. Raw CXC responses, source file receipts, usable time segments,
and individual exclusions are retained. Raw EVT2 files are available by receipt.
"""
import argparse
from collections import defaultdict
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
from tempfile import TemporaryDirectory
import urllib.parse
import xml.etree.ElementTree as ET

from astropy.time import Time
from archive import FIELDS, TAP, download, dump, fetch_events
from vet import separation
from visits import POINTING_RADIUS_ARCSEC, _lead_id, _measure


CSC = "https://cda.cfa.harvard.edu/csccli/getProperties"
HERE = Path(__file__).resolve().parent


def _query(where, endpoint, parameters, required):
    """Archive the actual HTTP entity, including service error documents."""
    where.parent.mkdir(parents=True, exist_ok=True)
    rec = {"endpoint": endpoint, "parameters": parameters}
    # Use GET for the CSC service, exactly as in its documented wget examples.
    encoded = urllib.parse.urlencode(parameters)
    url, body = ((endpoint + "?" + encoded, None) if endpoint == CSC
                 else (endpoint, encoded.encode("ascii")))
    xml = where.with_suffix(".xml")
    try:
        rec["response"] = download(url, xml, data=body)
        raw = xml.read_bytes()
        packed = where.with_suffix(".xml.gz")
        packed.write_bytes(gzip.compress(raw, mtime=0))
        rec["response"].update(raw_entity_gzip=packed.name,
                               gzip_sha256=hashlib.sha256(packed.read_bytes()).hexdigest())
        xml.unlink()
        tree = ET.fromstring(raw)
        states = [x for x in tree.findall(".//{*}INFO") if x.get("name","").upper()=="QUERY_STATUS"]
        if endpoint == TAP and not states:
            raise ValueError("CXC TAP response omits query status")
        if any(x.get("value","").upper() != "OK" for x in states):
            raise ValueError("query incomplete: " + str([(x.get("value"),x.text) for x in states])[:500])
        tables = tree.findall(".//{*}TABLE")
        if len(tables) != 1 or tables[0].find(".//{*}TABLEDATA") is None:
            raise ValueError("expected one VOTable TABLEDATA result")
        table = tables[0]
        cols = [str(e.get("name") or e.get("ID") or "").lower().split(".")[-1]
                for e in table.findall("{*}FIELD")]
        if not set(required).issubset(cols) or len(cols) != len(set(cols)):
            raise ValueError("unexpected result columns " + str(cols))
        rows = [dict(zip(cols,((e.text or "").strip() for e in tr),strict=True))
                for tr in table.findall(".//{*}TR")]
        rec.update(status="ok",row_count=len(rows),rows=rows)
    except Exception as exc:
        rec.update(status="failed", error=f"{type(exc).__name__}: {exc}")
    dump(where.with_suffix(".json"), rec)
    return rec


def _bounding_box(ra, dec, radius, prefix):
    """Padded rectangular prefilter for a spherical cap, including RA wrap."""
    low = max(-90.,dec-radius-1.e-7)
    high = min(90.,dec+radius+1.e-7)
    parts = [f"({prefix}.dec BETWEEN {low:.9f} AND {high:.9f})"]
    if low > -90 and high < 90:
        width=math.degrees(math.asin(min(1.,math.sin(math.radians(radius))/math.cos(math.radians(dec)))))+1.e-7
        left,right=(ra-width)%360,(ra+width)%360
        if left<right:
            parts.append(f"({prefix}.ra BETWEEN {left:.9f} AND {right:.9f})")
        else:
            parts.append(f"({prefix}.ra >= {left:.9f} OR {prefix}.ra <= {right:.9f})")
    return " AND ".join(parts)


def _pointings(folder, lead, cutoff):
    ra,dec=lead["ra_deg"],lead["dec_deg"]
    pred=_bounding_box(ra,dec,POINTING_RADIUS_ARCSEC/3600,"o")
    query=("SELECT "+", ".join("o."+k for k in FIELDS)+" FROM cxc.observation o "
           "WHERE o.instrument LIKE 'ACIS%' AND "+pred)
    rec=_query(folder/"pointings",TAP,{"REQUEST":"doQuery","LANG":"ADQL","FORMAT":"votable",
               "MAXREC":50000,"QUERY":query},FIELDS)
    accepted=defaultdict(list)
    excluded=[]
    for r in rec.get("rows",[]):
        oid=int(r["obsid"])
        distance=separation(ra,dec,float(r["ra"]),float(r["dec"]))
        if distance>POINTING_RADIUS_ARCSEC:
            continue
        reason=None
        if oid==int(lead["obsid"]):
            reason="own trigger observation"
        elif r["status"].upper()!="ARCHIVED":
            reason="not marked archived"
        elif not r["public_avail_date"] or r["public_avail_date"][:10]>cutoff:
            reason="unknown or later public release date"
        elif r["grating"].upper() in {"HETG","LETG"} or r["exposure_mode"].upper() in {"CC","CONTINUOUS","CONTINUOUS CLOCKING"}:
            reason="grating or continuous-clock observation"
        entry={**r,"pointing_separation_arcsec":distance}
        if reason:
            excluded.append({**entry,"reason":reason})
        else:
            accepted[oid].append(entry)
    return {"status":rec["status"],"error":rec.get("error"),"raw_query":"pointings.json",
            "included":dict(accepted),"excluded":excluded}


def _csc(folder,lead):
    ra,dec=lead["ra_deg"],lead["dec_deg"]
    radius=max(15.,3*float(lead["position_error_arcsec"]))
    cond=(_bounding_box(ra,dec,radius/3600,"m")+
          f" AND dbo.cone_distance(m.ra,m.dec,{ra:.10f},{dec:.10f}) <= {radius/60:.12f}")
    master=("SELECT m.name,m.ra,m.dec,m.err_ellipse_r0,m.err_ellipse_r1,m.err_ellipse_ang,"
            "m.significance,m.likelihood_class,m.conf_flag,m.sat_src_flag,m.streak_src_flag,"
            "m.acis_num,m.hrc_num,m.flux_aper_b FROM master_source m WHERE "+cond)
    # CXC's released CSC 2.2 Source Observation Summary joins (cli/adql.html).
    visits=("SELECT DISTINCT m.name,m.ra,m.dec,a.match_type,o.obsid,o.obi,o.region_id,o.gti_obs,"
            "o.gti_end,o.livetime,o.detector,o.flux_significance_b,o.src_cnts_aper_b,o.src_rate_aper_b,"
            "o.flux_aper_b,o.var_index_b,o.conf_code,o.streak_src_flag FROM master_source m, "
            "master_stack_assoc a, stack_source s, stack_observation_assoc b, observation_source o "
            "WHERE m.name=a.name AND s.detect_stack_id=a.detect_stack_id AND s.region_id=a.region_id "
            "AND s.detect_stack_id=b.detect_stack_id AND s.region_id=b.region_id "
            "AND o.obsid=b.obsid AND o.obi=b.obi AND o.region_id=b.region_id AND "+cond)
    result={"release":"2.2","cone_arcsec":radius}
    for name,query,required in (("master",master,("name","ra","dec")),
                                ("detections",visits,("name","ra","dec","obsid","gti_obs"))):
        rec=_query(folder/("csc_"+name),CSC,{"version":"rel2.2","query":query,
                   "coordFormat":"decimal","floatFormat":"native","outputFormat":"votable"},required)
        out=[]
        for r in rec.get("rows",[]):
            row={**r,"separation_arcsec":separation(ra,dec,float(r["ra"]),float(r["dec"]))}
            if "obsid" in row:
                row["same_trigger_obsid"]=int(row["obsid"])==int(lead["obsid"])
            out.append(row)
        result[name]={"status":rec["status"],"error":rec.get("error"),"record":f"csc_{name}.json","rows":out}
    return result


def _visit(folder,oid,leads):
    folder.mkdir(parents=True,exist_ok=True)
    result={"obsid":oid,"lead_ids":[_lead_id(r) for r in leads],"products":[]}
    with TemporaryDirectory() as directory:
        work=Path(directory)
        try:
            products=fetch_events(oid,work)
            for receipt in products:
                path=Path(receipt["path"])
                row={"input":{k:v for k,v in receipt.items() if k!="path"}}
                result["products"].append(row)
                try:
                    row["measurements"]=_measure(path,oid,leads)
                    for measure in row["measurements"].values():
                        measure["visit_iso_utc"]=Time(measure["visit_mjd_utc"],format="mjd",scale="utc").isot.tolist()
                    row["status"]="measured_or_geometry_excluded"
                except Exception as exc:
                    row.update(status="unassessable",error=f"{type(exc).__name__}: {exc}")
            result["status"]="acquired"
        except Exception as exc:
            result.update(status="acquisition_failed",error=f"{type(exc).__name__}: {exc}")
        if (work/"acquisition.json").is_file():
            result["archive_attempts"]=json.loads((work/"acquisition.json").read_text())
        for response in work.glob("primary_*.html"):
            raw=response.read_bytes()
            (folder/(response.name+".gz")).write_bytes(gzip.compress(raw,mtime=0))
    dump(folder/"visit.json",result)
    return result


def run(leads, out, cutoff):
    if (os.getenv("GITHUB_ACTIONS")!="true" or
        os.getenv("GITHUB_REF")!="refs/heads/research/chandra-acis-fast-transients" or
        os.getenv("GITHUB_EVENT_NAME") not in {"workflow_dispatch","push"}):
        raise RuntimeError("Targeted retrieval requires the campaign workflow")
    out.mkdir(parents=True,exist_ok=True)
    if not leads or len(set(_lead_id(r) for r in leads))!=len(leads):
        raise ValueError("expected unique, nonempty fixed-position lead list")
    summary={"release_cutoff_utc_inclusive":cutoff,"pointing_radius_arcsec":POINTING_RADIUS_ARCSEC,
             "git_commit":os.getenv("GITHUB_SHA"),"selected_leads":leads,"targets":{},"visits":[],
             "code_sha256":{s:hashlib.sha256((HERE/s).read_bytes()).hexdigest() for s in
                            ("target_history.py","visits.py","search.py","archive.py","vet.py")},
             "limitations":"Observation-catalog target coordinates are a pointing proxy. The wide cone "
             "and independent detector geometry are required. Missing catalog coordinates, offset "
             "pointings, HRC, gratings and CC mode limit coverage. Aperture counts use only supported "
             "dither/GTI bins. The 95% limit applies to mean detected aperture rate in those bins, "
             "without a PSF, vignetting or bad-pixel correction. Raw single-visit excess probabilities "
             "include chance excesses and diffuse emission. CSC matches are positional alternatives. "
             "A CSC per-observation photometry row can be a weak or forced measurement."}
    plan=defaultdict(list)
    for i,lead in enumerate(leads):
        label=f"target-{i:02d}"
        folder=out/label
        sky=[lead["ra_deg"],lead["dec_deg"],lead["position_error_arcsec"],lead["event_mjd"]]
        if not all(math.isfinite(v) for v in sky) or not 0<=sky[0]<360 or not -90<=sky[1]<=90 or sky[2]<=0:
            raise ValueError("invalid fixed target coordinate or epoch")
        p=_pointings(folder,lead,cutoff)
        dump(folder/"selection.json",p)
        try:
            csc=_csc(folder,lead)
        except Exception as exc:
            csc={"status":"failed","error":f"{type(exc).__name__}: {exc}"}
        dump(folder/"csc.json",csc)
        key=_lead_id(lead)
        summary["targets"][key]={"folder":label,"pointing_query_status":p["status"],
                                 "nearby_obsids":sorted(p["included"]),"supported":[],"gaps":[]}
        for oid in p["included"]:
            plan[oid].append(lead)
    dump(out/"history.json",summary)
    for oid,target_rows in sorted(plan.items()):
        folder=out/"visits"/str(oid)
        v=_visit(folder,oid,target_rows)
        record=(folder/"visit.json").relative_to(out).as_posix()
        summary["visits"].append({"obsid":oid,"status":v["status"],"record":record})
        for r in target_rows:
            key=_lead_id(r)
            target=summary["targets"][key]
            if not v["products"]:
                target["gaps"].append({"obsid":oid,"reason":v["status"],"record":record})
            for prod in v["products"]:
                meas=prod.get("measurements",{}).get(key,{})
                stamp={"obsid":oid,"record":record,"filename":prod["input"]["filename"],
                       "input_sha256":prod["input"]["sha256"],"status":meas.get("status",prod["status"])}
                if meas.get("status")=="measured_aperture":
                    stamp.update({k:meas[k] for k in ("visit_iso_utc","temporal_relation","accepted_live_s",
                                 "gross_counts","annulus_counts","net_aperture_rate_ct_s","local_excess_p",
                                 "conditional_aperture_rate_upper95_ct_s")})
                    target["supported"].append(stamp)
                else:
                    target["gaps"].append(stamp)
        dump(out/"history.json",summary)
    for t in summary["targets"].values():
        trials=len(t["supported"])
        for product in t["supported"]:
            product["per_target_bonferroni_p"]=min(1.,trials*product["local_excess_p"])
    dump(out/"history.json",summary)


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--leads",required=True,type=Path)
    cli.add_argument("--out",required=True,type=Path)
    cli.add_argument("--release-through",default=datetime.now(timezone.utc).date().isoformat())
    args=cli.parse_args()
    datetime.strptime(args.release_through,"%Y-%m-%d")
    run(json.loads(args.leads.read_text()),args.out,args.release_through)
