"""Bounded, resumable prefill of the original independent-visit product checkpoints."""
import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
from unittest.mock import patch

import campaign
import followup
import ingest_primary as ingest
import visits

DEPS=("visits.py","search.py","archive.py","vet.py","requirements.txt")
need,sha,save=ingest.require,ingest.sha,ingest.save


def implementation():
    sources={k:sha(campaign.HERE/k) for k in DEPS}
    return hashlib.sha256(json.dumps(sources,sort_keys=True).encode()).hexdigest()


def key(code,receipt,selected):
    # Precisely the native visits.run_visits key, including original lead order.
    identity=[code,receipt["sha256"],[[visits._lead_id(r),r["ra_deg"],r["dec_deg"],r["position_error_arcsec"],
                                      r["event_mjd"]] for r in selected]]
    return hashlib.sha256(json.dumps(identity,separators=(",",":")).encode()).hexdigest()


def partition(jobs,n):
    buckets=[[] for _ in range(min(n,len(jobs)))]
    weight=[0]*len(buckets)
    costs=[max(1024**2,j["receipt"]["bytes"])*(1+len(j["selected"])) for j in jobs]
    for k in sorted(range(len(jobs)),key=lambda k:(-costs[k],k)):
        b=min(range(len(buckets)),key=weight.__getitem__)
        buckets[b].append(k)
        weight[b]+=costs[k]
    return [sorted(v) for v in buckets]


def plan(root,groups):
    leads=campaign.read(root/"leads.json")
    lookup={visits._lead_id(v):i for i,v in enumerate(leads)}
    need(len(lookup)==len(leads),"Duplicate corrected lead")
    inventory,matches,excluded=visits._receipt_plan(root,leads) if leads else ([],[],[])
    code=implementation()
    jobs=[]
    names=set()
    for oid,rec,rows in matches:
        k=key(code,rec,rows)
        name=f"{oid}-{k[:20]}.json"
        need(name not in names and int(rec["obsid"])==oid and ingest.HEX.fullmatch(rec["sha256"]),"Invalid visit membership")
        names.add(name)
        jobs.append({"obsid":oid,"receipt":rec,"selected":[lookup[visits._lead_id(r)] for r in rows],
                     "key":k,"record":name})
    return {"schema":1,"implementation":code,"adapter_sha256":sha(Path(__file__)),"leads":leads,
            "leads_sha256":sha(root/"leads.json"),"inventory_sha256":sha(root/"inventory/inventory.json"),
            "inventory_obsids":len(inventory),"unassessable":excluded,"jobs":jobs,"groups":partition(jobs,groups)}


def digest(work):
    return hashlib.sha256(json.dumps(work,sort_keys=True,separators=(",",":"),allow_nan=False).encode()).hexdigest()


def state(out,job,leads):
    path=out/job["record"]
    if not path.exists():
        return "missing"
    try:
        rec=campaign.read(path)
    except (ValueError,OSError):
        return "invalid_checkpoint"
    expected=[visits._lead_id(leads[i]) for i in job["selected"]]
    if (rec.get("key")!=job["key"] or rec.get("lead_ids")!=expected or rec.get("obsid")!=job["obsid"] or
        any(rec.get("input",{}).get(k)!=job["receipt"].get(k) for k in ("filename","bytes","sha256","fits_sha256"))):
        return "invalid_checkpoint"
    if rec.get("status")=="complete" and set(rec.get("measurements",{}))==set(expected):
        return "complete"
    return "error" if rec.get("status")=="error" else "invalid_checkpoint"


def preserve(path,directory):
    if path.is_file():
        dst=directory/(sha(path)+".json")
        if not dst.exists():
            ingest.install(path,dst)


def perform(root,out,area,job,leads,retained):
    previous=state(out,job,leads)
    if previous=="complete":
        return
    path=out/job["record"]
    preserve(path,area/"previous"/job["record"])
    if previous=="invalid_checkpoint":
        path.unlink()  # Exact earlier bytes remain in previous/.
    selected=[leads[n] for n in job["selected"]]
    oid,receipt=job["obsid"],job["receipt"]
    wrote,attempts,acquisition=[],[],[]
    old_dump,old_replay=visits.dump,visits.replay
    def write(p,item):
        if p.name!="independent_visits.json":
            need(p==path,"Native visit key differs")
            old_dump(p,item)
            wrote.append(p)
    def take(manifest,dest):
        need(campaign.read(manifest)==[receipt],"Wrong pinned product request")
        start=time.monotonic()
        try:
            if retained:
                followup.event_bytes(receipt,oid,dest,retained,attempts)
                return
            try:
                old_replay(manifest,dest)
                attempts.append({"kind":"pinned_archive","url":receipt["url"],"sha256":receipt["sha256"],"bytes":receipt["bytes"]})
            except Exception as exc:
                attempts.append({"kind":"pinned_archive","url":receipt["url"],"error":str(exc)})
                followup.event_bytes(receipt,oid,dest,(),attempts)
        finally:
            acquisition.append(time.monotonic()-start)
    def one(evidence,rows):
        need(evidence==root and rows==selected,"Native worker changed assigned visit")
        return [],[(oid,receipt,selected)],[]
    started=time.monotonic()
    try:
        with patch.object(visits,"_receipt_plan",one),patch.object(visits,"dump",write),patch.object(visits,"replay",take):
            visits.run_visits(root,selected,out)
        need(wrote==[path] and state(out,job,leads) in {"complete","error"},"Missing native checkpoint")
    finally:
        stamp=area/"acquisition"/job["record"]
        preserve(stamp,area/"previous_acquisition"/job["record"])
        save(stamp,{"key":job["key"],"status":state(out,job,leads),"attempts":attempts,"wall_s":time.monotonic()-started,
                    "acquisition_s":sum(acquisition),"run":os.environ["GITHUB_RUN_ID"],"attempt":os.environ["GITHUB_RUN_ATTEMPT"]})


def worker(root,file,group,selection,stop,retained):
    work=campaign.read_obs(file)
    need(file.parent.name==digest(work) and work["implementation"]==implementation() and
         work["adapter_sha256"]==sha(Path(__file__)) and work["leads_sha256"]==sha(root/"leads.json") and
         work["inventory_sha256"]==sha(root/"inventory/inventory.json"),"Changed independent-visit inputs or program")
    need(0<=group<len(work["groups"]) and len(selection)==len(set(selection)) and
         set(selection)<=set(work["groups"][group]),"Foreign visit unit")
    with tempfile.TemporaryDirectory(dir=os.getenv("RUNNER_TEMP")) as tmp:
        os.environ["TMPDIR"]=tmp
        tempfile.tempdir=tmp
        for k in selection:
            if time.time()>=stop:
                break
            perform(root,root/"independent_visits",file.parent,work["jobs"][k],work["leads"],retained)


def report(work,out,area,finished,timed,start):
    totals,volume=Counter(),Counter()
    rows={}
    for j in work["jobs"]:
        s=state(out,j,work["leads"])
        p=out/j["record"]
        rows[j["record"]]={"status":s,"sha256":sha(p) if p.is_file() else None}
        totals[s]+=1
        volume[s]+=j["receipt"]["bytes"]
    r={"plan_sha256":digest(work),"nearby_products":len(rows),"lead_product_pairs":sum(len(j["selected"]) for j in work["jobs"]),
       "record_status":dict(totals),"input_bytes_by_status":dict(volume),"worker_exit_codes":finished,
       "workers_at_deadline":timed,"wall_s":time.monotonic()-start,"run":os.environ["GITHUB_RUN_ID"],
       "attempt":os.environ["GITHUB_RUN_ATTEMPT"],"records":rows,
       "scope":"Execution progress. The native reducer computes per-lead comparisons across the full fixed visit set."}
    file=area/f"progress-{r['run']}-{r['attempt']}.json.gz"
    if file.exists():
        old=area/"previous_progress"/(sha(file)+".json.gz")
        if not old.exists():
            ingest.install(file,old)
    save(file,r,zipped=True)
    print(f"{len(rows)} native products. Statuses {dict(totals)}; compressed bytes {dict(volume)}.",flush=True)
    return totals["complete"]==len(rows) and not any(finished.values())


def run(root,n,groups,minutes,retained):
    need(1<=n<=4 and n<=groups<=192 and 0<minutes<=260,"Invalid parallelism or soft budget")
    start=time.monotonic()
    work=plan(root,groups)
    out=root/"independent_visits"
    area=out/"execution"/digest(work)
    file=area/"plan.json.gz"
    if file.exists():
        need(campaign.read_obs(file)==work,"Conflicting visit plan")
    else:
        save(file,work,zipped=True)
    fresh,retries=[],[]
    for g,items in enumerate(work["groups"]):
        cached=[(k,state(out,work["jobs"][k],work["leads"])) for k in items]
        for label,queue in (("new",fresh),("retry",retries)):
            selected=[k for k,s in cached if s=="missing"] if label=="new" else [
                k for k,s in cached if s not in {"complete","missing"}]
            if selected:
                queue.append((f"{label}-{g}",g,selected))
    queue=fresh+retries
    deadline=start+minutes*60
    active,done,timed={},{},[]
    print(f"{len(work['jobs'])} products, {sum(j['receipt']['bytes'] for j in work['jobs'])} unique compressed bytes, "
          f"{len(queue)} pending work batches, {n} worker slots.",flush=True)
    try:
        while (queue or active) and time.monotonic()<deadline:
            while queue and len(active)<n:
                label,g,indices=queue.pop(0)
                args=[sys.executable,"-u",str(Path(__file__).resolve()),"_worker","--evidence",str(root),
                      "--plan",str(file),"--group",str(g),"--indices",",".join(map(str,indices)),
                      "--stop-at",str(time.time()+max(0,deadline-time.monotonic()-15))]
                for p in retained:
                    args+=["--retained",str(p)]
                active[label]=subprocess.Popen(args)
            for label,p in list(active.items()):
                if p.poll() is not None:
                    done[label]=p.returncode
                    del active[label]
            if active:
                time.sleep(min(1,max(0,deadline-time.monotonic())))
    finally:
        for label,p in active.items():
            if p.poll() is None:
                timed.append(label)
                p.terminate()
        for label,p in active.items():
            try:
                p.wait(timeout=3)
            except subprocess.TimeoutExpired:
                p.kill()
                p.wait()
            done[label]=p.returncode
    if not report(work,out,area,done,timed,start):
        raise SystemExit("Partial visit products preserved. Retry remaining records.")


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("mode",choices=("run","_worker"))
    cli.add_argument("--evidence",type=Path,default=campaign.EVIDENCE)
    cli.add_argument("--workers",type=int,default=2)
    cli.add_argument("--groups",type=int,default=48)
    cli.add_argument("--budget-minutes",type=float,default=230)
    cli.add_argument("--retained",type=Path,action="append",default=[])
    cli.add_argument("--plan",type=Path)
    cli.add_argument("--group",type=int)
    cli.add_argument("--indices")
    cli.add_argument("--stop-at",type=float)
    a=cli.parse_args()
    need(os.getenv("GITHUB_ACTIONS")=="true" and os.getenv("GITHUB_REPOSITORY")==ingest.REPO and
         os.getenv("GITHUB_REF")=="refs/heads/"+ingest.BRANCH and
         os.getenv("GITHUB_WORKFLOW_REF","").split("@",1)[0]==ingest.REPO+"/"+ingest.FLOW,"Use the Chandra workflow")
    if a.mode=="run":
        run(a.evidence.resolve(),a.workers,a.groups,a.budget_minutes,a.retained)
    else:
        worker(a.evidence.resolve(),a.plan.resolve(),a.group,[int(v) for v in a.indices.split(",")],a.stop_at,a.retained)
