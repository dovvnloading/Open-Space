"""Evaluate or verify the preregistered held-out analysis on original primary worker bytes."""
import argparse
from collections import defaultdict
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile

import campaign
import followup
import heldout
import heldout_attest
import ingest_primary as ingest

need=ingest.require
sha=ingest.sha
EVALUATOR="171943c76e8bad8057e07c1667abf00d2aa1be2cf4a49832466285f3980ff69c"
OUTPUTS=("coverage.json","controls.json","survivors.json","trigger_audit.csv.gz","summary.json","README.md")
RECEIPT="original_result.json"
INPUTS="original_ledgers.json.gz"


def digest(data):
    return hashlib.sha256(data).hexdigest()


def encoded(value):
    return (json.dumps(value,sort_keys=True,separators=(",",":"),allow_nan=False)+"\n").encode()


def locate(root, member, earlier, checkpoints):
    """Only the upstream producer SHA and size decide byte identity, across any repair cycles."""
    logical=member["source"]
    name=Path(logical).name
    paths=[("producer",followup.safe(root,member["destination"]),None)]
    if logical in checkpoints:
        paths.append(("checkpoint",checkpoints[logical],None))
    paths.extend(("before",folder/"before"/"observations"/name,folder) for folder in earlier)
    h=member["sha256"]
    for kind,path,revision in paths:
        if not path.is_file() or path.is_symlink() or path.stat().st_size!=member["bytes"] or sha(path)!=h:
            continue
        if kind=="before":
            oid=name.split(".")[0]
            plan_file=revision/"plan.json"
            plan=campaign.read(plan_file)
            need(plan["obsids"][oid]["primary_sha256"]==h,"Repair plan does not authenticate its before-ledger")
            if (revision/"result.json").exists():
                result=campaign.read(revision/"result.json")
                need(result["plan_sha256"]==sha(plan_file) and result["changed_or_attempted"][oid]["before_sha256"]==h,
                     "Completed repair receipt disagrees with its before-ledger")
        if kind=="checkpoint":
            need(digest(ingest.blob(ingest.PLAN,ingest.PREFIX+"/"+logical))==h,
                 "Checkpoint before-copy differs from the immutable original")
        return path
    raise ValueError(f"Original producer ledger exists in the archive but cannot be verified locally: {logical} {h}")


def registration(folder, sample):
    receipt=folder/"registration/scheduler.json"
    r=campaign.read(receipt)
    need((r["source_run"],r["source_attempt"],r["source_code"],r["protocol_sha256"]) ==
         (ingest.RUN,ingest.ATTEMPT,ingest.CODE,heldout.PROTOCOL) and
         r["selection_sha256"]==sha(sample)==sha(folder/"selection.json"),"Source sample differs from preregistration")
    jobs={}
    for page in r["pages"]:
        path=followup.safe(folder/"registration",page["file"])
        raw=gzip.decompress(path.read_bytes())
        need(digest(raw)==page["response_sha256"] and len(raw)==page["response_bytes"],"Changed scheduler page")
        heldout_attest.read_jobs(raw,jobs)
    need(set(jobs)==set(range(96,192)) and r["groups"]=={str(k):v for k,v in jobs.items()} and
         all(heldout_attest.queued(v) for v in jobs.values())==r["all_96_jobs_still_queued_without_assigned_runner_or_steps"],
         "Original scheduler rows or attestation changed")
    return {"sha256":sha(receipt),"registered_utc":r["registered_utc"],
            "all_jobs_queued_at_attestation":r["all_96_jobs_still_queued_without_assigned_runner_or_steps"]}


def original_campaign(root,view):
    need(heldout.COMMIT==ingest.PLAN and sha(Path(heldout.__file__))==EVALUATOR and
         sha(campaign.HERE/"HELD_OUT.md")==heldout.PROTOCOL,"Changed registered selection or evaluator")
    ingest.immutable(ingest.PLAN)
    ingest.immutable(ingest.CODE)
    name=ingest.PREFIX+"/inventory/campaign.json"
    raw=ingest.blob(ingest.PLAN,name)
    freeze=json.loads(raw)
    original={**followup.BEFORE,"vet.py":followup.OLD}
    need(freeze["implementation"]==original and freeze["inventory_sha256"]==followup.INVENTORY and
         freeze["n_controls"]==24 and freeze["n_injections_per_product"]==12 and
         digest(ingest.blob(ingest.CODE,name))==digest(raw) and
         {**campaign.identity(),"vet.py":followup.OLD}==original and
         campaign.identity()["vet.py"] in (followup.OLD,followup.NEW),"Original survey implementation changed")
    folder=view/"inventory"
    folder.mkdir(parents=True)
    (folder/"campaign.json").write_bytes(raw)
    file=root/"inventory/inventory.json"
    need(sha(file)==followup.INVENTORY,"Changed dated archive census")
    (folder/file.name).symlink_to(file)
    return digest(raw),original


def registry(root):
    area=root/ingest.HISTORY
    source=campaign.read(area/"source.json")
    coverage=campaign.read(area/"coverage.json")
    run=source["run"]
    need((run["id"],run["run_attempt"],run["head_sha"],run["head_branch"],run["status"])==
         (ingest.RUN,ingest.ATTEMPT,ingest.CODE,ingest.BRANCH,"completed") and
         run["path"].split("@")[0]==ingest.FLOW and source["plan_commit"]==ingest.PLAN and
         coverage["source_status"]==run and coverage["assignment_commit"]==ingest.PLAN and
         coverage["source_run"]==ingest.RUN and coverage["source_code"]==ingest.CODE,
         "Original source run and its terminal ingest audit are required")
    return source,coverage


def members(root,g,ids,source,coverage):
    """Authenticate this original ZIP and its member map without opening corrected ledgers."""
    k=str(g)
    name=f"chandra-evidence-{g}"
    entry=source["artifacts"].get(name)
    state=coverage["groups"][k]
    need(all(coverage["dated_obsids"][str(i)]["group"]==g for i in ids),"Wrong original group ownership")
    if entry is None:
        need(state["artifact_status"]=="missing_artifact" and state["artifact_id"] is None,
             "Inconsistent terminal missing-artifact receipt")
        return {},{},{"status":"missing_original_artifact"}
    rec_file=ingest.receipt_file(root,k)
    rec=campaign.read(rec_file)
    table=root/ingest.HISTORY/"members"/(k+".json.gz")
    need(state["artifact_status"]=="verified" and entry["id"]==state["artifact_id"]==rec["artifact_id"] and
         rec["group"]==k and rec["artifact_name"]==entry["name"]==name and
         rec["code_commit"]==ingest.CODE and rec["plan_commit"]==ingest.PLAN and
         entry["size_in_bytes"]==rec["zip_bytes"] and entry["digest"]==rec["zip_digest"] and
         re.fullmatch("sha256:[0-9a-f]{64}",entry["digest"]) and
         rec["assignment_sources"]==source["assignment_sources"] and
         sha(table)==rec["member_manifest_sha256"] and set(map(int,rec["worker_ledger_status"]))==set(ids),
         "Changed original artifact identity, assignment or member manifest")
    data=campaign.read_obs(table)
    need(len(data)==rec["source_member_count"],"Changed original artifact member count")
    result={}
    backups={}
    for member in data:
        logical=member["source"]
        if not logical.startswith("observations/"):
            continue
        m=re.fullmatch(r"observations/(\d+)\.json\.gz",logical)
        need(m and int(m[1]) in ids and logical not in result,"Foreign or duplicate producer ledger")
        allowed={logical,(ingest.HISTORY/"interrupted_worker"/k/logical).as_posix()}
        need(member["destination"] in allowed and re.fullmatch("[0-9a-f]{64}",member["sha256"]) and
             type(member["bytes"]) is int and member["bytes"]>0,"Wrong original producer member path or identity")
        result[logical]=member
        checkpoint=rec["checkpoint_before"].get(logical)
        if checkpoint and checkpoint["sha256"]==member["sha256"]:
            backups[logical]=root/ingest.HISTORY/"checkpoint_before"/k/logical
    for oid in ids:
        observed=rec["worker_ledger_status"][str(oid)]["worker_status"]
        need((observed!="unattempted")== (f"observations/{oid}.json.gz" in result)
             and coverage["dated_obsids"][str(oid)]["worker_status"]==observed,"Missing original worker ledger accounting")
    return result,backups,{"status":"verified","artifact_id":entry["id"],"artifact_name":name,
                           "zip_digest":entry["digest"],"zip_bytes":entry["size_in_bytes"],
                           "member_manifest_sha256":sha(table),"ingest_receipt_sha256":sha(rec_file)}


def mount_group(view,root,rows,source_members,backups,earlier):
    answer=[]
    for row in rows:
        oid=row["obsid"]
        logical=f"observations/{oid}.json.gz"
        m=source_members.get(logical)
        record={**row,"member":None,"included_bytes":False}
        if m:
            record["member"]={k:m[k] for k in ("source","destination","bytes","sha256")}
        # Exclusions and original absences are symmetric across all 25 streams.
        if m and not row["exclusions"]:
            file=locate(root,m,earlier,backups)
            (view/logical).symlink_to(file)
            record["included_bytes"]=True
        answer.append(record)
    return answer


def inputs(root,view,work):
    campaign_hash,science=original_campaign(root,view)
    selection=heldout.freeze(view,work)
    folder=root/"held_out"
    scheduler=registration(folder,work/"selection.json")
    registered=work/"assignment_sources"
    for f in registered.iterdir():
        need(sha(f)==sha(folder/"assignment_sources"/f.name),"Preregistered assignment source changed")
    source,coverage=registry(root)
    assignment=campaign.read(work/"selection.json")["sources"]
    shards=source["assignment_sources"]["inventory/shards.json"]
    need(assignment["evidence/inventory/shards.json"]["git_blob"]==shards["git_blob"] and
         assignment["evidence/inventory/shards.json"]["sha256"]==shards["sha256"],"Registry uses different fixed shards")
    earlier=sorted((root/followup.SCOPED).glob("*"))
    (view/"observations").mkdir()
    by_group=defaultdict(list)
    for r in selection:
        by_group[r["group"]].append(r)
    # Include every assigned group even when it has zero observations.
    groups,items={},[]
    for g in range(96,192):
        batch=by_group[g]
        live,backups,receipt=members(root,g,{r["obsid"] for r in batch},source,coverage)
        groups[str(g)]=receipt
        items.extend(mount_group(view,root,batch,live,backups,earlier))
    snapshot={"schema":1,"source_run":ingest.RUN,"source_attempt":ingest.ATTEMPT,"source_code":ingest.CODE,
              "assignment_commit":ingest.PLAN,"selection_sha256":sha(work/"selection.json"),
              "source_campaign_sha256":campaign_hash,"source_implementation":science,"protocol_sha256":heldout.PROTOCOL,
              "registration":scheduler,"groups":groups,"ledgers":items}
    return selection,snapshot


def coverage_check(rows,coverage):
    need(len(coverage)==len(rows),"Changed held-out coverage length")
    for original,r in zip(rows,coverage,strict=True):
        need(all(original[k]==r[k] for k in ("group","obsid","exclusions")), "Held-out cohort or order changed")
        h=original["member"]["sha256"] if original["included_bytes"] else None
        need(r["ledger_sha256"]==h,"Outcome read another version of its ledger")
        if original["exclusions"]:
            need(r["disposition"]=="previously_exposed","Previously exposed ObsID included")
        elif h is None:
            need(r["disposition"]=="unattempted","A later ledger filled an originally missing observation")
        else:
            need(r["disposition"] in {"incomplete","evaluable"},"Unaudited present producer ledger")


def run(root,mode):
    dest=root/"held_out"
    sealed=dest/RECEIPT
    if mode=="verify":
        need(sealed.is_file(),"Original held-out evaluation has not yet been sealed")
    with tempfile.TemporaryDirectory(dir=os.getenv("RUNNER_TEMP")) as tmp:
        base=Path(tmp)
        view=base/"view"
        work=base/"outputs"
        selection,snapshot=inputs(root,view,work)
        projected=encoded(snapshot)
        if sealed.is_file():
            s=campaign.read(sealed)
            need(s["schema"]==1 and s["program_sha256"]==sha(Path(__file__)) and s["evaluator_sha256"]==EVALUATOR
                 and s["input_json_sha256"]==digest(projected) and s["input_file_sha256"]==sha(dest/INPUTS)
                 and gzip.decompress((dest/INPUTS).read_bytes())==projected and
                 s["result_artifacts"]=={f:sha(dest/f) for f in OUTPUTS},"Original held-out inputs or result changed")
            summary=campaign.read(dest/"summary.json")
            need(summary["primary_implementation"]==snapshot["source_implementation"]
                 and summary["selection_sha256"]==snapshot["selection_sha256"],"Held-out campaign identity changed")
            coverage_check(snapshot["ledgers"],campaign.read(dest/"coverage.json"))
            return
        heldout.evaluate(selection,view,work)
        coverage_check(snapshot["ledgers"],campaign.read(work/"coverage.json"))
        doc=work/"README.md"
        with doc.open("a") as out:
            out.write("\n[Original producer ledgers and artifact checksums](original_ledgers.json.gz) · "
                      "[Sealed evaluation](original_result.json)\n"
                      "The source manifest identifies every original worker ledger. The replay command verifies its "
                      "archived bytes, including observations retained in pre-correction history.\n")
        (work/INPUTS).write_bytes(gzip.compress(projected,compresslevel=6,mtime=0))
        record={"schema":1,"program_sha256":sha(Path(__file__)),"evaluator_sha256":EVALUATOR,
                "input_json_sha256":digest(projected),"input_file_sha256":sha(work/INPUTS),"source_run":ingest.RUN,
                "source_code":ingest.CODE,"scope":"Original primary producer records from groups 96–191.",
                "result_artifacts":{f:sha(work/f) for f in OUTPUTS},
                "evaluation_run":os.environ["GITHUB_RUN_ID"],"evaluation_attempt":os.environ["GITHUB_RUN_ATTEMPT"],
                "evaluation_code":os.environ["GITHUB_SHA"],
                "checked_commit":ingest.git("rev-parse","HEAD").decode().strip()}
        (work/RECEIPT).write_bytes(encoded(record))
        # Preserve any earlier unsealed output, including an interrupted invocation.
        for name in (*OUTPUTS,INPUTS,RECEIPT):
            target=dest/name
            if target.is_file() and sha(target)!=sha(work/name):
                prior=dest/"prior_unsealed"/sha(target)/name
                if not prior.exists():
                    followup.copy(target,prior)
            followup.copy(work/name,target)
    print("Verified the original held-out producer manifest and sealed evaluation.")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode",choices=("evaluate","verify"))
    p.add_argument("--evidence",type=Path,default=campaign.EVIDENCE)
    args=p.parse_args()
    need(os.getenv("GITHUB_ACTIONS")=="true" and os.getenv("GITHUB_REPOSITORY")==ingest.REPO and
         os.getenv("GITHUB_REF")=="refs/heads/"+ingest.BRANCH,"Use the Chandra research workflow")
    run(args.evidence.resolve(),args.mode)
