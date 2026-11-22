"""Audit CSC 2.2 catalog variability at frozen real and scrambled trigger positions.

python results/chandra_acis_fxt/csc_variability.py --evidence PATH --out OUT --group 0
python results/chandra_acis_fxt/csc_variability.py --evidence PATH --out OUT --obsids 23324,32380
An explicit list can also be supplied as a JSON array of integers in --obsids-file.
Without a selector all deposited ledgers are processed. --replay uses only saved catalog replies.
"""
import argparse
from collections import Counter, defaultdict
import gzip
import hashlib
import json
import math
import os
from pathlib import Path

from archive import dump
from campaign import read_obs, save_obs
from target_history import CSC, _query
from vet import separation


# Column names and table joins follow the CXC 2.2 Source Observation Summary.
# https://cxc.cfa.harvard.edu/csc2.2/cli/adql.html
MASTER = ("SELECT m.name,m.ra,m.dec,m.err_ellipse_r0,m.err_ellipse_r1,m.err_ellipse_ang,"
          "m.significance,m.likelihood_class,m.conf_flag,m.sat_src_flag,m.streak_src_flag,"
          "m.acis_num,m.hrc_num,m.flux_aper_b FROM master_source m WHERE ")
OBSERVATION = ("SELECT DISTINCT m.name,m.ra,m.dec,a.match_type,o.obsid,o.obi,o.region_id,o.gti_obs,"
               "o.gti_end,o.livetime,o.detector,o.flux_significance_b,o.src_cnts_aper_b,o.src_rate_aper_b,"
               "o.flux_aper_b,o.var_index_b,o.conf_code,o.streak_src_flag FROM master_source m, "
               "master_stack_assoc a, stack_source s, stack_observation_assoc b, observation_source o "
               "WHERE m.name=a.name AND s.detect_stack_id=a.detect_stack_id AND s.region_id=a.region_id "
               "AND s.detect_stack_id=b.detect_stack_id AND s.region_id=b.region_id "
               "AND o.obsid=b.obsid AND o.obi=b.obi AND o.region_id=b.region_id AND ")
RELEASE = "rel2.2"
RELEASE_COVERAGE_END = "2023-12-31"


def valid(event):
    vals = [event.get(k) for k in ("ra_deg","dec_deg","position_error_arcsec")]
    return (all(isinstance(x,(int,float)) and math.isfinite(x) for x in vals)
            and 0<=vals[0]<360 and -90<=vals[1]<=90 and vals[2]>0)


def triggers(ledger):
    rows=[]
    ensembles=[]
    for prod in ledger["products"]:
        if "search" not in prod:
            continue
        search=prod["search"]
        groups=[("real",search["triggers"])]
        groups.extend((f"control-{c['replica']}",c["triggers"]) for c in search["controls"])
        inp=prod["input"]
        ensembles.append({"input_sha256":inp["sha256"],"filename":inp["filename"],
                          "completed_product":prod.get("complete"),"streams":[label for label,_ in groups]})
        for label,events in groups:
            for i,t in enumerate(events):
                if t.get("instrument_pass"):
                    rows.append({"stream":label,"trigger_key":t.get("trigger_key",f"{label}/{i}"),
                                 "input_sha256":inp["sha256"],"input_url":inp["url"],
                                 "input_filename":inp["filename"],"input_bytes":inp["bytes"],
                                 "frozen_vet_status":t.get("vet",{}).get("status"),
                                 **{k:t.get(k) for k in ("ra_deg","dec_deg","event_mjd","position_error_arcsec",
                                                        "on_counts","peak_counts","seed_score")}})
    return rows,ensembles


def search_radius(event):
    return max(15.,3.*event["position_error_arcsec"])


def sky_predicate(rows):
    """One enclosing rectangle per ObsID, padded for every trigger's review cone."""
    ra0=rows[0]["ra_deg"]
    low,high=90.,-90.
    left,right=180.,-180.
    full_ra=False
    for row in rows:
        ra,dec,r=row["ra_deg"],row["dec_deg"],search_radius(row)/3600
        low,high=min(low,dec-r),max(high,dec+r)
        if dec-r<=-90 or dec+r>=90:
            full_ra=True
            continue
        width=math.degrees(math.asin(min(1.,math.sin(math.radians(r))/math.cos(math.radians(dec)))))
        center=(ra-ra0+180)%360-180
        left,right=min(left,center-width),max(right,center+width)
    # Small outward paddings keep decimal serialization from clipping cone edges.
    lower=max(-90.,low-1.e-7)
    upper=min(90.,high+1.e-7)
    slab=f"(m.dec BETWEEN {lower:.9f} AND {upper:.9f})"
    if full_ra or right-left>=359.99999:
        return slab
    a=(ra0+left-1.e-7)%360
    b=(ra0+right+1.e-7)%360
    if a<=b:
        return slab+f" AND (m.ra BETWEEN {a:.9f} AND {b:.9f})"
    return slab+f" AND (m.ra >= {a:.9f} OR m.ra <= {b:.9f})"


def request(cache, query, required, replay):
    params={"version":RELEASE,"query":query,"coordFormat":"decimal","floatFormat":"native","outputFormat":"votable"}
    token=hashlib.sha256(json.dumps(params,sort_keys=True,separators=(",",":")).encode()).hexdigest()
    base=cache/token[:2]/token
    receipt=base.with_suffix(".json")
    if receipt.exists():
        rec=json.loads(receipt.read_text())
        if rec["parameters"]!=params:
            raise ValueError("catalog request cache identity mismatch")
        if rec["status"]=="ok":
            raw=gzip.decompress(base.with_suffix(".xml.gz").read_bytes())
            if hashlib.sha256(raw).hexdigest()!=rec["response"]["sha256"]:
                raise ValueError("catalog response bytes differ from archived receipt")
            return rec,receipt.relative_to(cache).as_posix()
    if replay:
        return {"status":"failed","error":"verified reply unavailable in replay"},receipt.relative_to(cache).as_posix()
    rec=_query(base,CSC,params,required)
    return rec,receipt.relative_to(cache).as_posix()


def value(raw):
    if raw not in ("",None):
        try:
            val=float(raw)
            return val if math.isfinite(val) else None
        except (ValueError,TypeError):
            pass
    return None


def characterize(event, masters, observations, original_status):
    if not valid(event):
        return {"status":"invalid_trigger_localization","matches":[]}
    reviewed=[]
    for name,master in masters.items():
        sep=separation(event["ra_deg"],event["dec_deg"],float(master["ra"]),float(master["dec"]))
        if sep>search_radius(event):
            continue
        major=value(master.get("err_ellipse_r0"))
        # Circular outer envelope of the published ellipse, for review only.
        limit=math.hypot(event["position_error_arcsec"],major) if major is not None else None
        associations=observations.get(name,[])
        close=limit is not None and sep<=limit
        variable=any((value(det.get("var_index_b")) or 0)>=8 for det in associations)
        unique_variable=any(det.get("match_type")=="u" and (value(det.get("var_index_b")) or 0)>=8
                            for det in associations)
        reviewed.append({"master":master,"separation_arcsec":sep,"positional_envelope_arcsec":limit,
                         "within_positional_envelope":close,"original_obsid_rows":associations,
                         "review_same_obsid_variability":close and variable,
                         "unique_association_var_index_at_least_8":close and unique_variable})
    reviewed.sort(key=lambda r:r["separation_arcsec"])
    if original_status!="ok":
        status="incomplete_catalog_queries"
    elif not reviewed:
        status="no_master_within_search_radius"
    elif not any(r["original_obsid_rows"] for r in reviewed):
        status="nearby_master_no_same_obsid_row"
    else:
        status="nearby_master_with_same_obsid_row"
    return {"status":status,"search_radius_arcsec":search_radius(event),"matches":reviewed}


def query_pair(oid, valid_rows, cache, replay):
    cond=sky_predicate(valid_rows)
    specs=[("masters",MASTER+cond,("name","ra","dec","err_ellipse_r0")),
           ("original_obsid",OBSERVATION+f"o.obsid={int(oid)} AND "+cond,
            ("name","ra","dec","obsid","var_index_b","flux_significance_b"))]
    receipts={}
    results={}
    for name,query,required in specs:
        try:
            result,path=request(cache,query,required,replay)
            receipts[name]={"status":result["status"],"error":result.get("error"),"receipt":path,
                            "raw_sha256":result.get("response",{}).get("sha256"),"row_count":result.get("row_count")}
            results[name]=result.get("rows",[])
        except Exception as exc:
            receipts[name]={"status":"failed","error":f"{type(exc).__name__}: {exc}"}
            results[name]=[]
    masters={r["name"]:r for r in results["masters"]}
    det=defaultdict(list)
    for r in results["original_obsid"]:
        if int(r["obsid"])!=int(oid):
            raise ValueError("returned CSC observation differs from request")
        det[r["name"]].append(r)
        # Preserve joined observation evidence even if the separate master query failed.
        masters.setdefault(r["name"],{k:r[k] for k in ("name","ra","dec")})
    ok=all(r["status"]=="ok" for r in receipts.values())
    return receipts,masters,det,("ok" if ok else "incomplete")


def audit(evidence, out, obsids, replay=False):
    if (os.environ.get("GITHUB_ACTIONS")!="true" or
        os.environ.get("GITHUB_REF")!="refs/heads/research/chandra-acis-fast-transients" or
        os.environ.get("GITHUB_EVENT_NAME") not in {"workflow_dispatch","push"}):
        raise RuntimeError("CSC audit requires the campaign workflow")
    out.mkdir(parents=True,exist_ok=True)
    cache=out/"catalog_replies"
    (out/"observations").mkdir(exist_ok=True)
    inv=evidence/"inventory/inventory.json"
    metadata={int(r["obsid"]):r for r in json.loads(inv.read_text())} if inv.exists() else {}
    summary={"csc_release":RELEASE,"release_date":"2026-10","public_observation_coverage_through":RELEASE_COVERAGE_END,
             "catalog_coordinate_system":"ICRS, original event and detection epochs kept separately",
             "github_commit":os.environ.get("GITHUB_SHA"),"implementation_sha256":{
                 n:hashlib.sha256(Path(__file__).with_name(n).read_bytes()).hexdigest()
                 for n in ("csc_variability.py","target_history.py","archive.py","vet.py")},
             "selected_obsids":sorted(set(obsids)),"observations":[],"limitations":
             "All instrument-passing real and control triggers use the same bounded cone and same-observation "
             "catalog join. The circular positional envelope is a review aid, not a false-association "
             "probability. Each association and its CSC quality flags are retained. Weak or forced CSC "
             "per-observation photometry is possible. Positional overlap and a catalog variability index "
             "are evidence about prior measurement, not the physical origin or burst epoch. Empty "
             "master or observation rows can reflect catalog coverage, sensitivity or confusion."}
    for oid in summary["selected_obsids"]:
        ledger=evidence/"observations"/f"{oid}.json.gz"
        report={"obsid":oid,"status":"not_deposited"}
        if ledger.is_file():
            try:
                obj=read_obs(ledger)
                if int(obj["obsid"])!=oid:
                    raise ValueError("ledger ObsID mismatch")
                events,ensembles=triggers(obj)
                loc=[r for r in events if valid(r)]
                receipts,masters,dets,status=query_pair(oid,loc,cache,replay) if loc else ({},{},{},"ok")
                states=Counter()
                priority=Counter()
                for e in events:
                    e["csc"]=characterize(e,masters,dets,status)
                    states[e["csc"]["status"]]+=1
                    if any(m["review_same_obsid_variability"] for m in e["csc"]["matches"]):
                        priority[e["stream"]]+=1
                release=metadata.get(oid,{}).get("release_date")
                details={"obsid":oid,"ledger_sha256":hashlib.sha256(ledger.read_bytes()).hexdigest(),
                         "ledger_status":obj["status"],"cxc_public_release_date":release,
                         "public_after_csc_coverage":release[:10]>RELEASE_COVERAGE_END if release else None,
                         "requested_catalog":RELEASE,"query_receipts":receipts,"ensembles":ensembles,
                         "events":events}
                dest=out/"observations"/f"{oid}.json.gz"
                save_obs(dest,details)
                report.update(status="audited" if status=="ok" else "catalog_incomplete",
                              record="observations/"+dest.name,sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),
                              instrument_pass_events=len(events),result_counts=dict(states),
                              same_obsid_variable_review_by_stream=dict(priority))
            except Exception as exc:
                report.update(status="failed",error=f"{type(exc).__name__}: {exc}")
        summary["observations"].append(report)
        dump(out/"csc_audit.json",summary)
    if not summary["selected_obsids"]:
        dump(out/"csc_audit.json",summary)
    return summary


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--evidence",type=Path,required=True)
    cli.add_argument("--out",type=Path,required=True)
    ch=cli.add_mutually_exclusive_group()
    ch.add_argument("--group",type=int)
    ch.add_argument("--obsids")
    ch.add_argument("--obsids-file",type=Path)
    cli.add_argument("--replay",action="store_true")
    a=cli.parse_args()
    if a.group is not None:
        groups=json.loads((a.evidence/"inventory/shards.json").read_text())
        if not 0<=a.group<len(groups):
            cli.error("invalid frozen group")
        ids=groups[a.group]
    elif a.obsids is not None:
        ids=[int(i.strip()) for i in a.obsids.split(",")]
    elif a.obsids_file is not None:
        ids=json.loads(a.obsids_file.read_text())
    else:
        ids=[int(p.name.split(".")[0]) for p in (a.evidence/"observations").glob("*.json.gz")]
    audit(a.evidence,a.out,list(map(int,ids)),a.replay)
