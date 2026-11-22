"""Offline notice and exposure replay of the immutable 27089 field follow-up."""
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import tarfile
import tempfile

from astropy.time import Time
from scipy.stats import chi2

import alerts
import archive
import campaign
import ingest_primary as ingest
import visits

ORIGINAL="30dfe8cf4027ac8c8fed691cc51bab672cb48930"
SUBDIR="45008126305-1"
GCN_ARCHIVE="6aafe58397639e773f00765e087c05f8071af00e9d94ebf2e2dc8b6bf145837c"
need,sha=ingest.require,ingest.sha


def encoded(obj):
    return (json.dumps(obj,sort_keys=True,indent=2,allow_nan=False)+"\n").encode()


def digest(b):
    return hashlib.sha256(b).hexdigest()


def put(p,b):
    p.parent.mkdir(parents=True,exist_ok=True)
    need(not p.is_symlink(),"Unexpected result symlink")
    if p.exists():
        need(sha(p)==digest(b),"Different bytes already occupy "+str(p))
        return
    with tempfile.NamedTemporaryFile(dir=p.parent,delete=False) as f:
        temp=Path(f.name)
        f.write(b)
    try:
        os.link(temp,p)
    except FileExistsError:
        need(sha(p)==digest(b),"Different concurrent result")
    finally:
        temp.unlink(missing_ok=True)


def wanted(review):
    out={}
    for row in (r for lead in review["leads"] for r in lead["matches"]):
        path,tar,h,ident=row["file"],row["archive_member"],row["sha256"],int(row["id"])
        ref=(tar,h,ident)
        need(ingest.HEX.fullmatch(h) and path==f"members/{ident}-{h[:16]}.json" and tar and
             not Path(tar).is_absolute() and ".." not in Path(tar).parts and out.get(path,ref)==ref,
             "Conflicting or unsafe notice identity")
        out[path]=ref
    return out


def copy_members(source,out,review):
    for name,(_,h,ident) in wanted(review).items():
        body=(source/name).read_bytes()
        need(digest(body)==h and int(json.loads(body)["circularId"])==ident,"Absent or altered circular")
        put(out/name,body)


class Snapshot:
    def __init__(self,root,field,src):
        self.root=root
        repo=campaign.HERE.parents[1]
        paths=[str(p.relative_to(repo)) for p in (field,src)]
        raw=subprocess.check_output(["git","-C",str(repo),"ls-tree","-r","-z",ORIGINAL,"--",*paths],
                                    env={**os.environ,"GIT_NO_LAZY_FETCH":"1"})
        self.blobs={}
        self.checked={}
        for line in raw.split(b"\0"):
            if line:
                meta,name=line.split(b"\t",1)
                mode,kind,h=meta.split()
                if (mode,kind)==(b"100644",b"blob"):
                    self.blobs[repo/name.decode()]=h.decode()

    def check(self,p):
        need(p in self.blobs and p.is_file() and not p.is_symlink(),"Missing original regular file: "+str(p))
        size=p.stat().st_size
        git=hashlib.sha1(f"blob {size}\0".encode())
        body=hashlib.sha256()
        with p.open("rb") as f:
            while block:=f.read(1024**2):
                git.update(block)
                body.update(block)
        need(git.hexdigest()==self.blobs[p],"Changed file in original field deposit: "+str(p))
        self.checked[p.relative_to(self.root).as_posix()]={"git_blob":git.hexdigest(),"sha256":body.hexdigest(),"bytes":size}
        return p


def restore(src,out,review,check,expected_tar):
    receipts=review["receipts"]
    for label in ("index","archive"):
        r=receipts[label]
        stem=label+"-"+digest(r["url"].encode())[:24]
        need(r["ok"] and ingest.HEX.fullmatch(r["sha256"]) and
             r["body_file"]==stem+"-"+r["sha256"][:20]+".bin" and campaign.read(check(src/(stem+".json")))==r,
             "Saved GCN HTTP memo differs")
        file=check(src/r["body_file"])
        need(sha(file)==r["sha256"] and file.stat().st_size==r["size_bytes"] and
             any(a.get("http_status")==200 and all(a.get(k)==r[k] for k in ("body_file","sha256","size_bytes"))
                 for a in r["attempts"]),"Changed GCN HTTP response")
    index,tar=receipts["index"],receipts["archive"]
    need(index["url"]==alerts.INDEX and alerts.archive_link((src/index["body_file"]).read_bytes())==tar["url"]
         and tar["sha256"]==expected_tar,"Archive snapshot differs")
    entries=wanted(review)
    match={a:(name,h,ident) for name,(a,h,ident) in entries.items()}
    need(len(match)==len(entries),"Duplicate requested tar member")
    found={}
    with tarfile.open(src/tar["body_file"],mode="r|gz") as z:
        for member in z:
            if member.name in match:
                name,h,ident=match[member.name]
                need(name not in found and member.isfile() and 0<member.size<=2000000,"Invalid shortlisted member")
                b=z.extractfile(member).read()
                need(digest(b)==h and len(b)==member.size and int(json.loads(b)["circularId"])==ident,"Changed circular")
                found[name]=b
    need(found.keys()==entries.keys(),"Missing shortlisted circular")
    for name,b in found.items():
        put(out/name,b)
    return found


def circulars(review,lead,bodies):
    [old]=review["leads"]
    need(all(old[k]==lead[k] for k in ("obsid","input_sha256","trigger_key","ra_deg","dec_deg","event_mjd","trial_duration_s")),
         "Lead identity differs")
    peak=Time(lead["event_mjd"],format="mjd",scale="utc").to_datetime(timezone=timezone.utc)
    rows=[]
    for m in old["matches"]:
        d=json.loads(bodies[m["file"]])
        title,body=d["subject"],d["body"]
        event=str(d.get("eventId") or "")
        published=datetime.fromtimestamp(float(d["createdOn"])/1000,timezone.utc)
        days=alerts.dates_in_name(event+" "+title)
        pos=list(alerts.positions(body))
        times=alerts.onset_times(body)
        replay=alerts.review(lead,peak,days,published,bool(alerts.EP.search(event+" "+title) or alerts.EP.search(body)),pos,times)
        need(replay and all(m[k]==v for k,v in replay.items()) and
             (m["event"],m["subject"],m["published_utc"])==(event,title,published.isoformat()),"Circular parser replay differs")
        rows.append({"id":m["id"],"file":m["file"],"archive_member":m["archive_member"],"sha256":m["sha256"],"subject":title,
            "shortlist_relation":m["relation"],"position_status":"parsed" if pos else "unlocalized_by_parser",
            "positions":[{**p,"offset_arcsec":alerts.separation(lead["ra_deg"],lead["dec_deg"],p["ra_deg"],p["dec_deg"])} for p in pos],
            "published_utc":published.isoformat(),"publication_minus_peak_s":(published-peak).total_seconds(),
            "event_name_days":[{"utc_day":v.isoformat(),"calendar_days_from_peak":(v-peak.date()).days} for v in sorted(days)],
            "explicit_onsets":[{"utc":v.isoformat(),"seconds_from_peak":(v-peak).total_seconds()} for v in times]})
    return {"peak_utc":peak.isoformat(),"target":old["status"],"notices":rows,
            "limits":"Positions, named event dates, publication dates and explicit onsets remain separate. Unlocalized entries "
            "lack coordinates recognized by the frozen parser. A circular may discuss unrelated sources."}


def checked_input(folder,oid,r,attempts,check):
    name=r["filename"]
    need(r["obsid"]==oid and ingest.HEX.fullmatch(r["sha256"]) and ingest.HEX.fullmatch(r["fits_sha256"]) and
         r["bytes"]>0 and archive.product_names(f'<a href="{name}">EVT2</a>',oid)==[name],"Invalid event receipt")
    archive.check_archive_url(r["url"])
    for i,a in enumerate(attempts):
        ix=a.get("index",{})
        if a["index_url"]!=r["index_url"] or ix.get("sha256")!=r["index_sha256"] or name not in a.get("selected_files",[]):
            continue
        matches=[p for p in a.get("downloaded",[]) if p["filename"]==name and all(p[k]==r[k] for k in ("url","sha256","bytes"))]
        need(len(matches)==1,"Missing EVT2 download receipt")
        listing=gzip.decompress(check(folder/f"primary_{i}.html.gz").read_bytes())
        need(digest(listing)==ix["sha256"] and len(listing)==ix["bytes"] and
             name in archive.product_names(listing.decode(errors="replace"),oid),"Original EVT2 index changed")
        return {k:r[k] for k in ("obsid","filename","url","sha256","fits_sha256","bytes","revision","index_sha256")}
    raise ValueError("No archived primary index for EVT2 receipt")


def combine(rows):
    by={}
    for group in ("before","after","overlapping_epoch"):
        selected=[r for r in rows if r["temporal_relation"]==group]
        live=math.fsum(r["accepted_live_s"] for r in selected)
        n=sum(r["gross_counts"] for r in selected)
        off=sum(r["annulus_counts"] for r in selected)
        scaled=math.fsum(r["annulus_counts"]/r["annulus_area_ratio"] for r in selected)
        upper=float(.5*chi2.ppf(.95,2*(n+1))) if live else None
        by[group]={"supported_products":len(selected),"distinct_obsids":len({r["obsid"] for r in selected}),
                   "summed_live_s":live,"gross_counts":n,"annulus_counts":off,"scaled_annulus_counts":scaled,
                   "mean_net_aperture_rate_ct_s":(n-scaled)/live if live else None,
                   "gross_poisson_upper95_counts":upper,"gross_poisson_upper95_ct_s":upper/live if live else None}
    return by


def visit_records(field,lead,execution,check):
    folder=field/"history"
    summary=campaign.read(check(folder/"history.json"))
    key=visits._lead_id(lead)
    need(summary["selected_leads"]==[lead] and set(summary["targets"])=={key},"Foreign visit target")
    t=summary["targets"][key]
    need(t["folder"]=="target-00" and t["pointing_query_status"]=="ok","Unresolved pointing query")
    selection=campaign.read(check(folder/"target-00/selection.json"))
    need(sorted(map(int,selection["included"]))==sorted(t["nearby_obsids"]) and
         {r["obsid"] for r in summary["visits"]}==set(t["nearby_obsids"]) and
         len(summary["visits"])==len(t["nearby_obsids"]) and lead["obsid"] not in t["nearby_obsids"],"Visit set changed")
    allrows,gaps,pins,actual=[],[],set(),{}
    ntrials=len(t["supported"])
    for v in summary["visits"]:
        oid=v["obsid"]
        rel=f"visits/{oid}/visit.json"
        need(v["record"]==rel,"Unrecognized visit record")
        path=check(folder/rel)
        obs=campaign.read(path)
        need((obs["obsid"],obs["lead_ids"],obs["status"])==(oid,[key],v["status"]),"Visit identities differ")
        if not obs["products"]:
            gaps.append({"obsid":oid,"record":rel,"reason":obs["status"]})
        for entry in obs["products"]:
            r=entry["input"]
            src=checked_input(path.parent,oid,r,obs.get("archive_attempts",[]),check)
            need(r["sha256"] not in pins,"Duplicate visit product")
            pins.add(r["sha256"])
            m=entry.get("measurements",{}).get(key,{})
            row={"obsid":oid,"record":rel,"filename":r["filename"],"input_sha256":r["sha256"],"status":m.get("status",entry["status"])}
            if row["status"]!="measured_aperture":
                gaps.append(row)
                continue
            ratio=m["annulus_area_ratio"]
            start_mjd,end_mjd=m["visit_mjd_utc"]
            epoch="before" if end_mjd<=lead["event_mjd"] else "after" if start_mjd>lead["event_mjd"] else "overlapping_epoch"
            need(ratio==12 and start_mjd<=end_mjd and m["event_mjd_utc"]==lead["event_mjd"] and
                 m["temporal_relation"]==epoch and m["accepted_live_s"]>0,
                 "Invalid measured aperture or exposure")
            spans=sorted((a,b,s["duty"]) for s in m["segments"] for a,b in s["accepted_bins"])
            need(all(a<b and 0<d<=1 for a,b,d in spans) and all(stop<=start for (_,stop,_),(start,_,_) in zip(spans,spans[1:]))
                 and math.isclose(m["accepted_live_s"],math.fsum((b-a)*d for a,b,d in spans),rel_tol=1e-10) and
                 all(m[k]==sum(s[k] for s in m["segments"]) and m[k]>=0 for k in ("gross_counts","annulus_counts")),
                 "Saved visit counts or GTI live exposure disagree")
            n,bg,live=m["gross_counts"],m["annulus_counts"],m["accepted_live_s"]
            need(math.isclose((n-bg/ratio)/live,m["net_aperture_rate_ct_s"],abs_tol=1e-16) and
                 math.isclose(float(.5*chi2.ppf(.95,2*(n+1)))/live,m["conditional_aperture_rate_upper95_ct_s"]),
                 "Saved visit rate or limit differs")
            for k in ("visit_iso_utc","temporal_relation","accepted_live_s","gross_counts","annulus_counts",
                      "net_aperture_rate_ct_s","local_excess_p","conditional_aperture_rate_upper95_ct_s"):
                row[k]=m[k]
            row["per_target_bonferroni_p"]=min(1.,ntrials*row["local_excess_p"])
            actual[(oid,r["sha256"])]=row
            allrows.append({**row,"annulus_area_ratio":ratio,"input":src,"visit_record_sha256":sha(path)})
    need(len(actual)==ntrials and all(actual.get((r["obsid"],r["input_sha256"]))==r for r in t["supported"])
         and gaps==t["gaps"] and execution["visit_coverage"]=={"pointing_query_status":"ok","nearby_obsids":len(t["nearby_obsids"]),
         "supported_products":len(allrows),"coverage_gaps":len(gaps)},"Gaps, trial counts or measured products disagree")
    return {"summary_by_epoch":combine(allrows),"supported":allrows,"per_target_trials":ntrials,
            "pointing_exclusions":selection["excluded"],"nearby_obsids":t["nearby_obsids"],"coverage_gaps":gaps,
            "method":"Sum gross and scaled annulus counts across supported product livetime. Treat all gross counts as source "
            "for the one-sided exact Poisson upper bound. These are exposure-weighted mean detected aperture rates over "
            "separate observations; short peaks can be substantially brighter. Every single-visit test retains its original "
            "multiplicity across all supported products. Original Git blobs, exact archive directory listings and EVT2 "
            "download receipts are verified here. Compressed and decompressed EVT2 checksums identify the remote source bytes."}


def render(notices,v):
    def f(n): return "—" if n is None else f"{n:.6g}"
    lines=["# Archived field follow-up for ObsID 27089","",
           "Independent visits give the following exposure-weighted mean aperture rates at the Chandra flare position.","",
           "| Epoch | Products | Live time (s) | Gross counts | Scaled annulus counts | Mean net rate (ct/s) | Gross 95% upper rate (ct/s) |",
           "|---|---:|---:|---:|---:|---:|---:|"]
    for group,b in v["summary_by_epoch"].items():
        lines.append(f"| {group} | {b['supported_products']} | {f(b['summed_live_s'])} | {b['gross_counts']} | "
                     f"{f(b['scaled_annulus_counts'])} | {f(b['mean_net_aperture_rate_ct_s'])} | {f(b['gross_poisson_upper95_ct_s'])} |")
    rows=notices["notices"]
    miss=sum(r["position_status"]=="unlocalized_by_parser" for r in rows)
    lines+=["",f"The original comparison included {v['per_target_trials']} supported product tests and retained "
            f"{len(v['coverage_gaps'])} coverage gaps. The rates average the supported intervals in other observations. "
            "Their upper bounds apply to mean detected aperture counts over those intervals.", "",
            f"The archived GCN shortlist contains {len(rows)} circulars. The frozen parser found no labelled position in "
            f"{miss} of them. Publication dates, named event days, explicit onsets and angular offsets are kept separately "
            "for each circular. The source text remains available for manual interpretation.","",
            "| Circular | Parsed offsets (arcsec) | Publication minus peak (s) | Named day offsets | Explicit onset offsets (s) |",
            "|---|---|---:|---|---|"]
    for r in rows:
        coord=", ".join(f(p["offset_arcsec"]) for p in r["positions"]) or "Unlocalized"
        dates=", ".join(str(p["calendar_days_from_peak"]) for p in r["event_name_days"]) or "—"
        onsets=", ".join(f(p["seconds_from_peak"]) for p in r["explicit_onsets"]) or "—"
        lines.append(f"| [{r['id']}](../../{r['file']}) | {coord} | {f(r['publication_minus_peak_s'])} | {dates} | {onsets} |")
    lines+=["","[Circular offsets, dates and exact source paths](notices.json) · [Visit inputs, exclusions and multiplicity](visits.json) · "
            "[Original and derived checksums](manifest.json)",""]
    return "\n".join(lines).encode()


def backfill(root,run):
    import field_followup as field
    need(run==SUBDIR,"Only the pinned original 27089 field run is supported")
    out=root/"field_followup"/str(field.OBSID)/run
    src=root/"priority"/str(ingest.RUN)/"alerts/source"
    snap=Snapshot(root,out,src)
    check=snap.check
    report=campaign.read(check(out/"gcn_review.json"))
    execution=campaign.read(check(out/"execution.json"))
    [lead]=campaign.read(check(out/"leads.json"))
    [inputrec]=campaign.read(check(out/"event_input.json"))
    need(execution["source_artifact_id"]==field.ARTIFACT and execution["source_zip_digest"]==field.ZIP_SHA and
         execution["primary_run"]==ingest.RUN and execution["source_workflow_commit"]==ingest.CODE and
         execution["gcn_source_directory"]==src.relative_to(root).as_posix() and
         (lead["obsid"],lead["input_sha256"],lead["trigger_key"])==(field.OBSID,field.PRODUCT,field.TRIGGER) and
         inputrec["sha256"]==field.PRODUCT and report["leads_sha256"]==sha(out/"leads.json") and
         report["status"]=="parsed_archive" and len(report["leads"][0]["matches"])==13,"Different original field follow-up")
    bodies=restore(src,out,report,check,GCN_ARCHIVE)
    n=circulars(report,lead,bodies)
    v=visit_records(out,lead,execution,check)
    folder=out/"source_replay"/sha(Path(__file__))[:16]
    payload={"notices.json":encoded(n),"visits.json":encoded(v),"README.md":render(n,v)}
    for k,b in payload.items():
        put(folder/k,b)
    manifest={"schema":1,"original_commit":ORIGINAL,"original_followup":SUBDIR,"primary_run":ingest.RUN,
              "original_files":snap.checked,"recovered_members":{k:digest(b) for k,b in sorted(bodies.items())},
              "derived_files":{k:digest(b) for k,b in payload.items()},
              "code":{k:sha(campaign.HERE/k) for k in ("field_replay.py","alerts.py","visits.py","requirements.txt")},
              "location_of_member_paths":"Relative to the original field directory, two levels above this manifest."}
    put(folder/"manifest.json",encoded(manifest))
    print((folder/"README.md").read_text(),flush=True)
