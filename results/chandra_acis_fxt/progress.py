"""Summarize the checksum-verified dated Chandra evidence present in this checkout."""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import io
import os
from pathlib import Path
import subprocess

import campaign
import ingest_primary as ingest
import report
from search import DURATIONS, _time_support

LENGTHS = tuple(str(n) for n in DURATIONS)
SUPPORT = tuple("eligible_start_s_" + n for n in LENGTHS)
OBS_FIELDS = ("obsid","artifact_group","ledger_origin","ledger_sha256","status","scanned_all_products",
              "complete_followup","acquired_products","real_rows_recorded","detector_pass_rows_recorded",
              "real_rows_included","detector_pass_rows_included","distinct_detector_passes_included",*SUPPORT,"errors")
INPUT_FIELDS = ("obsid","filename","url","bytes","sha256","fits_sha256","revision","cycle",
                "included_in_support","ledger_origin","ledger_sha256")


def local_git(*args, binary=True):
    # A partial clone must report a missing object rather than doing a lazy network fetch.
    return subprocess.check_output(["git","-C",str(campaign.HERE.parents[1]),*args], text=not binary,
                                   env={**os.environ,"GIT_NO_LAZY_FETCH":"1","GIT_ALLOW_PROTOCOL":""})


def checkpoint_tree():
    raw = local_git("ls-tree","-r","-z",ingest.PLAN+":"+ingest.PREFIX+"/observations")
    out = {}
    for entry in raw.split(b"\0"):
        if entry:
            metadata, name = entry.split(b"\t",1)
            mode, kind, blob = metadata.split()
            ingest.require(kind == b"blob" and name.endswith(b".json.gz"), "Unexpected checkpoint observation tree")
            out[int(name.split(b".")[0])] = blob.decode()
    return out


def native_hashes(path, with_git):
    checksum = hashlib.sha256()
    oid = hashlib.sha1(usedforsecurity=False) if with_git else None
    if oid:
        oid.update(f"blob {path.stat().st_size}\0".encode())
    with path.open("rb") as f:
        while chunk := f.read(1024 * 1024):
            checksum.update(chunk)
            if oid:
                oid.update(chunk)
    return checksum.hexdigest(), oid.hexdigest() if oid else None


def registered(root, audit, source):
    known, groups = {}, []
    hist = root/ingest.HISTORY
    ingest.require(audit["source_run"]==ingest.RUN and audit["source_code"]==ingest.CODE
                   and audit["assignment_commit"]==ingest.PLAN and audit["source_status"]==source["run"],
                   "Refresh the offline audit for this artifact index")
    for g in map(str,range(192)):
        recpath=ingest.receipt_file(root,g)
        group=audit["groups"][g]
        verified=group["artifact_status"]=="verified"
        ingest.require(recpath.exists()==verified, "Refresh audit after importing group "+g)
        if not verified:
            continue
        rec=campaign.read(recpath)
        entry=source["artifacts"]["chandra-evidence-"+g]
        ingest.require(rec["artifact_id"]==group["artifact_id"]==entry["id"] and rec["zip_digest"]==entry["digest"],
                       "Producer identity differs from the ingest audit")
        manifest=hist/"members"/(g+".json.gz")
        ingest.require(ingest.sha(manifest)==rec["member_manifest_sha256"], "Original artifact membership changed")
        for m in campaign.read_obs(manifest):
            p=Path(m["source"])
            if p.parent==Path("observations") and m["source"]==m["destination"]:
                oid=int(p.name.split(".")[0])
                ingest.require(oid not in known and audit["dated_obsids"][str(oid)]["group"]==int(g),
                               "Conflicting primary observation ownership")
                known[oid]=(m["sha256"],m["bytes"])
        groups.append({"group":int(g),"artifact_id":entry["id"],"zip_digest":entry["digest"],
                       "created_at":entry["created_at"],"member_manifest_sha256":rec["member_manifest_sha256"]})
    return known, groups


def describe(row):
    acquired={(x["filename"],x["sha256"]) for x in row.get("acquired",[])}
    completed={(p["input"]["filename"],p["input"]["sha256"]) for p in row["products"] if p.get("complete")}
    finished=bool(acquired) and acquired==completed and row["status"] in ("complete","followup_incomplete")
    intervals=defaultdict(list)
    recorded=passed=0
    positions=[]
    for product in row["products"]:
        s=product.get("search") or {}
        events=s.get("triggers",[])
        recorded+=len(events)
        passed+=sum(e.get("instrument_pass") is True for e in events)
        if finished:
            for d in LENGTHS:
                intervals[d].extend(s["start_support_by_duration"][d])
            for e in events:
                if e["instrument_pass"]:
                    positions.append({**{k:e.get(k) for k in ("seed_score","start_met_s","stop_met_s","ra_deg","dec_deg",
                                                            "position_error_arcsec")},"obsid":row["obsid"]})
    support={field:sum((b-a)*duty for a,b,duty in _time_support(intervals[d]))
             for d,field in zip(LENGTHS,SUPPORT)}
    distinct=None if any(None in p.values() for p in positions) else len(report.distinct(positions))
    return {"scanned_all_products":finished,"complete_followup":finished and row["status"]=="complete",
            "acquired_products":len(acquired),"real_rows_recorded":recorded,"detector_pass_rows_recorded":passed,
            "real_rows_included":recorded if finished else 0,"detector_pass_rows_included":passed if finished else 0,
            "distinct_detector_passes_included":distinct,**support}


def csv_gz(path,fields,rows):
    temp=path.with_name(path.name+".new")
    with temp.open("wb") as f, gzip.GzipFile(filename="",fileobj=f,mode="wb",mtime=0) as gz, \
            io.TextIOWrapper(gz,encoding="utf-8",newline="") as txt:
        writer=csv.DictWriter(txt,fieldnames=fields,extrasaction="ignore",lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temp.replace(path)


def render(summary, early_link, supplement_link):
    n=summary
    lines=["# Chandra primary search progress","",
           f"**All event products have been scanned in {n['all_product_scans']:,} of {n['census_obsids']:,} dated ObsIDs. "
           f"The accumulated eligible start time for 60-second flares is {n['start_support_s']['60']/1e6:,.6f} Ms.**","",
           f"The deposited files contain {n['known_ledgers']:,} authenticated observation ledgers. Catalog checks and "
           f"required photometry finished in {n['complete_followup_obsids']:,} ObsIDs. The tally spans the immutable "
           f"checkpoint and {n['verified_primary_artifacts']} of 192 primary worker artifacts. "
           f"The latest included source material is dated {n['latest_source_utc']}.","",
           f"The completed scans produced {n['real_trigger_rows']:,} real trigger rows. "
           f"{n['detector_pass_rows']:,} passed the instrumental tests before counterpart review and cross-product deduplication. "
           "The observation table keeps trigger rows from partial scans in separate columns.","",
           "| Trial duration (s) | Eligible flare-start live time (s) |","|---:|---:|"]
    lines.extend(f"| {d} | {n['start_support_s'][d]:,.1f} |" for d in LENGTHS)
    lines.extend(["","Support is the union of available start intervals within each ObsID with its live-time fraction, summed across "
                      "ObsIDs. Counts and support above require scans of every acquired product, including the 24 scrambled streams "
                      "and the 12 scheduled injection draws per product.","",
                  f"{n['census_obsids']-n['all_product_scans']:,} dated ObsIDs remain outside this exposure. Which fields appear in "
                  "a checkpoint depends on worker scheduling, run time and archive or catalog availability. "
                  "These totals describe the deposited fields. The 14 supplemental ObsIDs with missing release dates are reported separately.","",
                  "[Observation status and trigger counts](progress_obsids.csv.gz) · [Input URLs and checksums](progress_inputs.csv.gz) · "
                  "[Machine summary](progress.json) · [Artifact coverage and gaps](coverage.json)","",
                  f"[Prespecified 128-ObsID cohort]({early_link}) · [Separate 14-ObsID supplement]({supplement_link}) · "
                  f"[Original Actions run](https://github.com/{ingest.REPO}/actions/runs/{ingest.RUN})"])
    return "\n".join(lines)+"\n"


def run(root):
    hist=root/ingest.HISTORY
    audit=campaign.read(hist/"coverage.json")
    source=campaign.read(hist/"source.json")
    inventory=root/"inventory/inventory.json"
    ingest.require(ingest.sha(inventory)==ingest.INVENTORY,"Different release census")
    population=[e["obsid"] for e in campaign.read(inventory)]
    ingest.require(len(set(population))==len(population) and set(map(str,population))==set(audit["dated_obsids"]),
                   "Census and audit differ")
    imported,groups=registered(root,audit,source)
    original=checkpoint_tree()
    statuses,origins=Counter(),Counter()
    observations,inputs=[],[]
    for oid in sorted(population):
        path=root/"observations"/f"{oid}.json.gz"
        group=audit["dated_obsids"][str(oid)]["group"]
        row={"obsid":oid,"status":"unattempted","products":[],"errors":[]}
        origin="no_ledger"
        checksum=None
        if path.exists():
            checksum, git_oid=native_hashes(path,oid in original and oid not in imported)
            if oid in imported:
                ingest.require((checksum,path.stat().st_size)==imported[oid],"Imported ledger changed for "+str(oid))
                origin="primary_worker"
            else:
                ingest.require(oid in original and git_oid==original[oid],"Ledger lacks matching import or checkpoint identity: "+str(oid))
                origin="checkpoint" if group is None else "checkpoint_with_primary_retry_assigned"
            row=ingest.ledger(path,oid)
        else:
            ingest.require(oid not in imported and oid not in original,"A previously deposited ledger is missing: "+str(oid))
        entry={"obsid":oid,"artifact_group":group,"ledger_origin":origin,"ledger_sha256":checksum,"status":row["status"],
               "errors":" | ".join(row["errors"]),**describe(row)}
        observations.append(entry)
        statuses[row["status"]]+=1
        origins[origin]+=1
        for rec in row.get("acquired",[]):
            inputs.append({**rec,"obsid":oid,"included_in_support":entry["scanned_all_products"],"ledger_origin":origin,
                           "ledger_sha256":checksum})
    csv_gz(hist/"progress_obsids.csv.gz",OBS_FIELDS,observations)
    csv_gz(hist/"progress_inputs.csv.gz",INPUT_FIELDS,inputs)
    def total(field):
        return sum(r[field] for r in observations)
    timestamps=[r["created_at"] for r in groups]
    checkpoint_utc=local_git("show","-s","--format=%cI",ingest.PLAN,binary=False).strip()
    from datetime import datetime, timezone
    latest=max(datetime.fromisoformat(t.replace("Z","+00:00")) for t in [checkpoint_utc,*timestamps])
    latest_utc=latest.astimezone(timezone.utc).isoformat(timespec="seconds").replace("+00:00","Z")
    summary={"source_run":ingest.RUN,"source_code":ingest.CODE,"checkpoint_commit":ingest.PLAN,
             "latest_source_utc":latest_utc,"census_obsids":len(population),"known_ledgers":sum(n for k,n in origins.items() if k!="no_ledger"),
             "all_product_scans":total("scanned_all_products"),"complete_followup_obsids":total("complete_followup"),
             "verified_primary_artifacts":len(groups),"status_counts":dict(statuses),"ledger_origins":dict(origins),
             "acquired_product_receipts":len(inputs),"included_product_receipts":sum(i["included_in_support"] for i in inputs),
             "real_trigger_rows":total("real_rows_included"),"detector_pass_rows":total("detector_pass_rows_included"),
             "distinct_detector_passes":None if any(r["distinct_detector_passes_included"] is None for r in observations)
                  else total("distinct_detector_passes_included"),
             "real_trigger_rows_in_partial_scans":total("real_rows_recorded")-total("real_rows_included"),
             "detector_pass_rows_in_partial_scans":total("detector_pass_rows_recorded")-total("detector_pass_rows_included"),
             "start_support_s":{d:total(field) for d,field in zip(LENGTHS,SUPPORT)},
             "files_sha256":{p:ingest.sha(hist/p) for p in ("progress_inputs.csv.gz","progress_obsids.csv.gz","coverage.json","source.json")},
             "source_artifacts":groups, "inventory_sha256":ingest.INVENTORY, "label":"Measured"}
    summary["catalog_headline"]=(f"Interim Chandra ACIS search covered {summary['all_product_scans']:,} of {len(population):,} dated ObsIDs "
                                 f"and {summary['start_support_s']['60']/1e6:.6f} Ms of eligible 60-second flare-start time.")
    ingest.save(hist/"progress.json",summary)
    early=f"../../provisional/moderate-count-v1/{ingest.RUN}-{ingest.ATTEMPT}/selection.json"
    ingest.require((hist/early).is_file(),"Missing prespecified 128-ObsID cohort selection")
    supplemental=f"../../supplemental_search/{ingest.RUN}-{ingest.ATTEMPT}/README.md"
    if not (hist/supplemental).is_file():
        artifact=source["artifacts"].get(ingest.SUP)
        supplemental=(f"https://github.com/{ingest.REPO}/actions/runs/{ingest.RUN}/artifacts/{artifact['id']}"
                      if artifact else "../../release_supplement/supplement.json")
    page=hist/"PROGRESS.md"
    staging=page.with_name("PROGRESS.md.new")
    staging.write_text(render(summary,early,supplemental))
    staging.replace(page)
    return summary


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence",type=Path,default=campaign.EVIDENCE)
    args=parser.parse_args()
    ingest.require(os.getenv("GITHUB_ACTIONS")=="true" and os.getenv("GITHUB_REPOSITORY")==ingest.REPO
                   and os.getenv("GITHUB_REF")=="refs/heads/"+ingest.BRANCH,"Use the Chandra research workflow")
    print(run(args.evidence.resolve())["catalog_headline"])
