"""Materialize the registered sample and archive the primary scheduler state before evaluation."""
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

import campaign
import heldout
import ingest_primary as ingest


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def read_jobs(raw, mapping):
    page=json.loads(raw)
    for job in page["jobs"]:
        match=re.fullmatch(r"observation \(([0-9]+)\)",job["name"])
        if not match:
            continue
        ident=int(match[1])
        if not 96<=ident<192:
            continue
        if (job["run_id"]!=ingest.RUN or job["head_sha"]!=ingest.CODE or
                job.get("run_attempt")!=ingest.ATTEMPT or ident in mapping):
            raise ValueError("Held-out scheduler provenance or group identity changed")
        mapping[ident]={k:job.get(k) for k in ("id","name","status","conclusion","runner_id","started_at","completed_at")}
        mapping[ident]["step_count"]=len(job.get("steps") or [])
    return len(page["jobs"])


def queued(job):
    # GitHub sets started_at to the enqueue time even on never-assigned matrix jobs.
    return (job["status"]=="queued" and not job["runner_id"] and not job["conclusion"]
            and not job["completed_at"] and job["step_count"]==0)


def run():
    if (os.getenv("GITHUB_ACTIONS")!="true" or os.getenv("GITHUB_REPOSITORY")!=ingest.REPO or
            os.getenv("GITHUB_REF")!="refs/heads/"+ingest.BRANCH):
        raise RuntimeError("Use the Chandra research workflow")
    root=campaign.EVIDENCE
    dest=root/"held_out"
    sample=heldout.freeze(root,dest)
    folder=dest/"registration"
    folder.mkdir(parents=True,exist_ok=True)
    receipt=folder/"scheduler.json"
    selection_sha=heldout.sha(dest/"selection.json")
    jobs={}
    if receipt.exists():
        state=campaign.read(receipt)
        if state["selection_sha256"]!=selection_sha:
            raise ValueError("The previously attested assignment changed")
        for item in state["pages"]:
            raw=gzip.decompress((folder/item["file"]).read_bytes())
            if digest(raw)!=item["response_sha256"]:
                raise ValueError("The saved scheduling response changed")
            read_jobs(raw,jobs)
        if {str(k):v for k,v in jobs.items()}!=state["groups"]:
            raise ValueError("The saved scheduler rows disagree with the preserved pages")
        return
    pages=[]
    for page in range(1,10):
        route=f"repos/{ingest.REPO}/actions/runs/{ingest.RUN}/attempts/{ingest.ATTEMPT}/jobs?per_page=100&page={page}"
        raw=subprocess.check_output(["gh","api",route])
        count=read_jobs(raw,jobs)
        path=folder/f"jobs-{page}.json.gz"
        path.write_bytes(gzip.compress(raw,mtime=0))
        pages.append({"api_path":route,"file":path.name,"response_sha256":digest(raw),"response_bytes":len(raw)})
        if count<100:
            break
    if set(jobs)!=set(range(96,192)):
        raise ValueError("Missing or repeated registered groups")
    all_queued=all(queued(v) for v in jobs.values())
    state={"checked_utc":datetime.now(timezone.utc).isoformat(),"source_run":ingest.RUN,
           "source_attempt":ingest.ATTEMPT,"source_code":ingest.CODE,"selection_sha256":selection_sha,
           "protocol_sha256":heldout.PROTOCOL,"registered_utc":"2026-11-22T06:35:00Z",
           "attestation_run":os.environ["GITHUB_RUN_ID"],"attestation_attempt":os.environ["GITHUB_RUN_ATTEMPT"],
           "attestation_code":os.environ["GITHUB_SHA"],"groups":dict(sorted(jobs.items())),"pages":pages,
           "all_96_jobs_still_queued_without_assigned_runner_or_steps":all_queued,
           "assigned_rows":len(sample),"included_after_fixed_exclusions":sum(not r["exclusions"] for r in sample),
           "scope":"This records the later materialization of the exact sample and the returned GitHub scheduling "
                   "state. The detailed statistical protocol and source assignment were already fixed separately."}
    campaign.dump(receipt,state)
    print("Registered sample:",state["included_after_fixed_exclusions"],"ObsIDs; all groups queued:",all_queued)


if __name__=="__main__":
    run()
