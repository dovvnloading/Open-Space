"""Resumable final reduction of the frozen primary run, executed by the Chandra workflow."""
import argparse
import csv
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from unittest.mock import patch

import campaign
import final_assembly
import followup
import ingest_primary as ingest

HERE=campaign.HERE
ROOT=campaign.EVIDENCE
MARK="stage_manifest.json.gz"
SEAL="stage_seal.json"
SUP=ROOT/"supplemental_search"/f"{ingest.RUN}-{ingest.ATTEMPT}"
TABLES=("summary.json","obsid_coverage.csv","event_receipts.csv","control_trials.csv","leads.json","completeness.csv")
SCRIPTS=("finalize.py","final_assembly.py","followup.py","catalog_pack.py","heldout.py","publication.py","candidate_inputs.py")
need=ingest.require
sha=ingest.sha
save=ingest.save


def use(*args):
    subprocess.run([sys.executable,"-u",*map(str,args)],check=True)


def tree(path, omit=()):
    result={}
    for p in sorted(path.rglob("*")):
        rel=p.relative_to(path).as_posix()
        if rel in omit:
            continue
        need(not p.is_symlink(),"A stage cannot reference a symlink: "+rel)
        if p.is_file():
            result[rel]=[p.stat().st_size,sha(p)]
    return result


def stamped(path,mode,meta):
    files=tree(path,{MARK,SEAL})
    save(path/MARK,files,zipped=True)
    rec={"mode":mode,"files_sha256":sha(path/MARK),"file_count":len(files),
         "stage_code_sha256":sha(Path(__file__)),"producer_run":os.environ["GITHUB_RUN_ID"],
         "producer_attempt":os.environ["GITHUB_RUN_ATTEMPT"],**meta}
    save(path/SEAL,rec)
    return rec


def verified(path,mode):
    rec=campaign.read(path/SEAL)
    need(rec["mode"]==mode and rec["files_sha256"]==sha(path/MARK)
         and rec["stage_code_sha256"]==sha(Path(__file__)),"Mismatched finalization stage")
    files=campaign.read_obs(path/MARK)
    need(len(files)==rec["file_count"] and tree(path,{MARK,SEAL})==files,"Stage contents differ from sealed bytes")
    return rec


def terminal():
    hist=ROOT/ingest.HISTORY
    source=campaign.read(hist/"source.json")
    census=campaign.read(hist/"coverage.json")
    need(source["run"]["id"]==ingest.RUN and source["run"]["head_sha"]==ingest.CODE
         and source["run"]["status"]=="completed" and census["source_status"]==source["run"]
         and census["assignment_commit"]==ingest.PLAN and sha(ROOT/"inventory/inventory.json")==followup.INVENTORY,
         "Refresh the source artifact index and terminal ingest audit")
    for name,g in ingest.EXPECTED.items():
        state=census["groups"][g]["artifact_status"]
        if name in source["artifacts"]:
            need(state=="verified" and census["groups"][g]["artifact_id"]==source["artifacts"][name]["id"]
                 and ingest.receipt_file(ROOT,g).is_file(),"Still-available source artifact needs importing: "+name)
        else:
            need(state=="missing_artifact","Undocumented missing original source artifact: "+name)
    need((SUP/"summary.json").is_file() and (SUP/"release_bounds.json").is_file(),
         "Recover the original supplemental archive or rerun the fourteen cases separately before finalizing")
    return {"source_run":ingest.RUN,"source_code":ingest.CODE,"assignment_commit":ingest.PLAN,
            "source_registry_sha256":sha(hist/"source.json"),"source_coverage_sha256":sha(hist/"coverage.json"),
            "inventory_sha256":sha(ROOT/"inventory/inventory.json"),
            "supplemental_summary_sha256":sha(SUP/"summary.json"),"gaps":census["transport_gaps"],
            "incomplete_original_groups":census["incomplete_worker_groups"]}


def independent(path):
    base=path.resolve()
    need(not base.is_relative_to(ROOT) and not ROOT.is_relative_to(base),"Use a staging directory outside deposited evidence")
    return base


def plan_job(bundle,groups):
    context=terminal()
    if (bundle/SEAL).exists():
        check_plan(bundle)
    else:
        bundle.parent.mkdir(parents=True,exist_ok=True)
        with tempfile.TemporaryDirectory(dir=bundle.parent) as tmp:
            stage=Path(tmp)/"plan"
            followup.prepare(ROOT,stage,groups)
            save(stage/"source_context.json",context)
            stamped(stage,"plan",{"plan_sha256":sha(stage/"plan.json")})
            if bundle.exists():
                backup=Path(tempfile.mkdtemp(dir=bundle.parent,prefix="unsealed-plan-"))
                bundle.rename(backup/"previous")
            stage.rename(bundle)
    plan=campaign.read(bundle/"plan.json")
    groups=list(range(len(plan["groups"])))
    if os.getenv("GITHUB_OUTPUT"):
        with open(os.environ["GITHUB_OUTPUT"],"a") as f:
            print("matrix="+json.dumps(groups),file=f)
            print("has_work="+str(bool(groups)).lower(),file=f)


def check_plan(bundle):
    seal=verified(bundle,"plan")
    need(seal["plan_sha256"]==sha(bundle/"plan.json") and campaign.read(bundle/"source_context.json")==terminal(),
         "Plan provenance and terminal archive differ")
    plan=campaign.read(bundle/"plan.json")
    need(plan["repair_sha256"]==sha(HERE/"followup.py") and plan["catalog_reader_sha256"]==sha(HERE/"catalog_pack.py"),
         "Changed retry implementation")
    return plan


def owner(name, run, attempt):
    p=Path(name)
    first=p.parts[0]
    if first=="observations" and len(p.parts)==2 and p.name.endswith(".json.gz"):
        return p.name[:-8]
    if first=="photometry" and len(p.parts)>=3:
        return p.parts[1]
    if first=="catalogs" and len(p.parts)>=4 and p.parts[1]=="followup-"+run+"-"+attempt:
        return p.parts[2]
    if first=="repair_records":
        if len(p.parts)==2 and p.name.endswith(".json"):
            return p.name[:-5]
        if len(p.parts)>=3:
            return p.parts[1]
    raise ValueError("Unexpected worker output: "+name)


def worker_ownership(files,group,plan,run,attempt):
    allowed={str(i) for i in plan["groups"][group]}
    for name in files:
        need(owner(name,run,attempt) in allowed,"Foreign worker file: "+name)


def seal_worker(bundle,group,stage):
    plan=check_plan(bundle)
    need(0<=group<len(plan["groups"]),"Unknown residual group")
    stage.mkdir(parents=True,exist_ok=True)
    ri,at=os.environ["GITHUB_RUN_ID"],os.environ["GITHUB_RUN_ATTEMPT"]
    for area in ("repair_records","observations"):
        for path in (stage/area).glob("*"):
            match=re.fullmatch(r"(\d+)\.json\.(?:new|tmp)",path.name)
            if path.is_file() and match:
                oid=int(match[1])
                need(oid in plan["groups"][group],"Unowned interrupted write")
                dest=stage/"repair_records"/str(oid)/"interrupted_writes"/area/(sha(path)+"-"+path.name)
                dest.parent.mkdir(parents=True,exist_ok=True)
                path.replace(dest)
    worker_ownership(tree(stage,{MARK,SEAL}),group,plan,ri,at)
    receipts={}
    for path in (stage/"repair_records").glob("*.json"):
        r=campaign.read(path)
        oid=r["obsid"]
        need(str(oid)+".json"==path.name and r["run_id"]==ri and r["run_attempt"]==at and
             r["plan_sha256"]==sha(bundle/"plan.json") and r["primary_sha256"]==plan["obsids"][str(oid)]["primary_sha256"],
             "Worker receipt is inconsistent with the sealed plan")
        if r.get("revised_sha256"):
            need(sha(followup.old_path(stage,oid))==r["revised_sha256"],"Changed revised worker ledger")
        receipts[str(oid)]=sha(path)
    return stamped(stage,"workers",{"group":group,"plan_sha256":sha(bundle/"plan.json"),"receipts":receipts,
                                     "missing":sorted(set(plan["groups"][group])-set(map(int,receipts)))})


def cycle_path():
    freeze=campaign.read(ROOT/"inventory/campaign.json")
    need(freeze["implementation"]==campaign.identity() and freeze["implementation"]["vet.py"]==followup.NEW,
         "A completed corrected reduction is required")
    rev=followup.safe(ROOT,freeze["followup_revision"])
    final=campaign.read(rev/"result.json")
    need(rev.is_relative_to(ROOT/followup.SCOPED) and final["plan_sha256"]==sha(rev/"plan.json"),"Missing corrected reduction")
    return rev,ROOT/"finalization"/rev.name


def stage_copy(src,dst):
    dst.parent.mkdir(parents=True,exist_ok=True)
    try:
        os.link(src,dst)
    except OSError:
        shutil.copyfile(src,dst)


def reduce(bundle,incoming,partial):
    plan=check_plan(bundle)
    ph=sha(bundle/"plan.json")
    stages={}
    selected={}
    for path in sorted(incoming.rglob(SEAL)):
        if campaign.read(path).get("mode")!="workers":
            continue
        d=path.parent
        rec=verified(d,"workers")
        g=rec["group"]
        need(0<=g<len(plan["groups"]) and rec["plan_sha256"]==ph,"Alien worker stage")
        files=campaign.read_obs(d/MARK)
        worker_ownership(files,g,plan,rec["producer_run"],rec["producer_attempt"])
        seen={Path(n).stem for n in files if Path(n).parent==Path("repair_records")}
        need(seen==set(rec["receipts"]) and all(rec["receipts"][n]==sha(d/"repair_records"/(n+".json")) for n in seen),
             "Missing or extra worker completion receipt")
        need(sorted(set(plan["groups"][g])-set(map(int,seen)))==rec["missing"],"Inconsistent missing-worker accounting")
        stage_id=sha(path)
        if stage_id in stages:
            continue
        stages[stage_id]=(rec,d,files)
        for n in seen:
            record=campaign.read(d/"repair_records"/(n+".json"))
            need(record["obsid"]==int(n) and record["plan_sha256"]==ph and record["run_id"]==rec["producer_run"]
                 and record["run_attempt"]==rec["producer_attempt"] and record["primary_sha256"]==plan["obsids"][n]["primary_sha256"],
                 "Supplied worker record is outside the sealed identity")
            rank=(record["outcome"]=="complete",bool(record.get("revised_sha256")),int(rec["producer_run"]),int(rec["producer_attempt"]))
            current=selected.get(n)
            if current and current[0][0] and rank[0]:
                need(current[2]==record.get("revised_sha256"),"Two completed worker outputs disagree for "+n+
                     "; choose the appropriate sealed stage or run a new corrective plan")
            if not current or rank>current[0]:
                selected[n]=(rank,stage_id,record.get("revised_sha256"))
    unresolved={str(g):[oid for oid in ids if str(oid) not in selected] for g,ids in enumerate(plan["groups"])}
    unresolved={g:ids for g,ids in unresolved.items() if ids}
    need(partial or not unresolved,"Workers missing. Retry the missing ObsIDs or use --checkpoint-partial for a bounded repair.")
    tag=plan["run_id"]+"-"+plan["run_attempt"]
    logs=ROOT/"finalization"/tag
    logs.mkdir(parents=True,exist_ok=True)
    summary={"plan_sha256":ph,"source_context":campaign.read(bundle/"source_context.json"),
             "plan_seal":campaign.read(bundle/SEAL),"received":{s:r for s,(r,_,_) in stages.items()},
             "chosen_stage_by_obsid":{n:item[1] for n,item in selected.items()},"missing_worker_obsids":unresolved}
    selection=logs/"reduction_inputs.json"
    if selection.exists():
        need(campaign.read(selection)==summary,"An earlier reduction used different worker stages; start a continuation plan")
    save(selection,summary)
    followup.copy(bundle/MARK,logs/"stage_manifests"/"plan.json.gz")
    with tempfile.TemporaryDirectory(dir=os.getenv("RUNNER_TEMP")) as temp:
        staging=Path(temp)
        for stage_id,(rec,d,files) in stages.items():
            followup.copy(d/MARK,logs/"stage_manifests"/(stage_id+".json.gz"))
            for name in files:
                n=owner(name,rec["producer_run"],rec["producer_attempt"])
                if n in selected and selected[n][1]==stage_id:
                    stage_copy(d/name,staging/stage_id/name)
                else:
                    # Preserve receipts, negatives and interrupted non-selected attempts without merging namespaces.
                    target=logs/"alternative_workers"/stage_id/name
                    if target.exists():
                        need(sha(target)==files[name][1],"Conflicting historical worker output")
                    else:
                        followup.copy(d/name,target)
        result=followup.finish(ROOT,bundle,staging,True)
    save(logs/"reduction_status.json",{"plan_sha256":ph,"result_sha256":sha(cycle_path()[0]/"result.json"),
                                       "residual_obsids":result["remaining_obsids"],"missing_worker_obsids":unresolved})
    return result


def provenance():
    rev,out=cycle_path()
    return {"corrected_revision":rev.name,"revision_receipt_sha256":sha(rev/"result.json"),
            "plan_sha256":sha(rev/"plan.json"),"campaign_sha256":sha(ROOT/"inventory/campaign.json"),
            "implementation_sha256":{name:sha(HERE/name) for name in SCRIPTS},
            "inventory_sha256":sha(ROOT/"inventory/inventory.json")}


def repeatable(receipt,key):
    if not receipt.is_file():
        return False
    old=campaign.read(receipt)
    return old["key"]==key and old["status"]=="complete" and all(
        (ROOT/name).is_file() and sha(ROOT/name)==h for name,h in old["outputs"].items())


def output_receipt(path,key,outputs,status="complete",errors=()):
    save(path,{"key":key,"status":status,"errors":list(errors),
               "outputs":{p.relative_to(ROOT).as_posix():sha(p) for p in outputs if p.is_file()}})


def measured_report():
    import report
    rev,out=cycle_path()
    receipt=out/"measured_report.json"
    key=provenance()
    plan=campaign.read(rev/"plan.json")
    revised=campaign.read(rev/"result.json")["changed_or_attempted"]
    for n,prior in plan["obsids"].items():
        p=followup.old_path(ROOT,int(n))
        expected=revised.get(n,{}).get("after_sha256",prior["primary_sha256"])
        need((sha(p) if p.exists() else None)==expected,"Census ledger changed since the repair: "+n)
    if repeatable(receipt,key):
        return
    with tempfile.TemporaryDirectory(dir=os.getenv("RUNNER_TEMP")) as tmp:
        stage=Path(tmp)
        for tree_name in ("inventory","observations"):
            (stage/tree_name).symlink_to(ROOT/tree_name,target_is_directory=True)
        with patch.object(report,"EVIDENCE",stage),patch.object(report,"HERE",stage),patch.object(report,"update_catalog",lambda *_:None), \
                patch.dict(os.environ,{"GITHUB_RUN_ID":plan["run_id"]}):
            report.run()
        need(all((stage/n).is_file() for n in TABLES),"Missing primary reduction products")
        outputs=[]
        for p in stage.iterdir():
            if p.is_file() and p.name!="README.md":
                followup.archive_replace(p,ROOT/p.name,out/"before_report"/p.name)
                outputs.append(ROOT/p.name)
        output_receipt(receipt,key,outputs)
    print("Regenerated full corrected primary measurement tables; the catalog awaits written adjudication.")


def csc_context(out):
    import csc_variability
    area=ROOT/"csc_variability"
    with (ROOT/"obsid_coverage.csv").open(newline="") as f:
        ids=[int(r["obsid"]) for r in csv.DictReader(f) if r["scanned_all_products"].lower()=="true"]
    records={}
    pending=[]
    failures=[]
    def check(oid):
        ledger=followup.old_path(ROOT,oid)
        cached=area/"observations"/f"{oid}.json.gz"
        if cached.is_file():
            r=campaign.read_obs(cached)
            if r["ledger_sha256"]==sha(ledger):
                success=all(q.get("status")=="ok" for q in r["query_receipts"].values())
                records[str(oid)]={"ledger_sha256":r["ledger_sha256"],"record_sha256":sha(cached),
                                   "status":"complete" if success else "partial"}
                return success
        return False
    for oid in ids:
        if not check(oid):
            pending.append(oid)
    for start in range(0,len(pending),100):
        batch=pending[start:start+100]
        for oid in batch:
            cached=area/"observations"/f"{oid}.json.gz"
            if cached.is_file():
                history=out/"before_csc"/str(oid)/(sha(cached)+".json.gz")
                if not history.exists():
                    followup.copy(cached,history)
        csc_variability.audit(ROOT,area,batch)
        receipt=area/"csc_audit.json"
        followup.copy(receipt,out/"csc_batches"/(sha(receipt)+".json"))
        for oid in batch:
            if not check(oid):
                failures.append(str(oid))
                records.setdefault(str(oid),{"status":"missing_current_csc_record"})
    index=out/"csc_records.json.gz"
    save(index,records,zipped=True)
    return [index],failures


def context(task):
    rev,out=cycle_path()
    measured=out/"measured_report.json"
    need(measured.is_file() and repeatable(measured,provenance()),"First regenerate the corrected primary measurement")
    script={"csc":"csc_variability.py","optical":"lead_optical.py"}.get(task,task+".py")
    key={"measurement_sha256":sha(measured),"leads_sha256":sha(ROOT/"leads.json"),
         "program_sha256":sha(HERE/script)}
    if task=="review":
        for f in (ROOT/"independent_visits/independent_visits.json",ROOT/"alerts/review.json"):
            key[f.name]=sha(f) if f.is_file() else None
    path=out/(task+".json")
    if repeatable(path,key):
        return
    errors=[]
    files=[]
    try:
        if task=="optical":
            area=ROOT/"lead_optical"
            use(HERE/"lead_optical.py","--evidence",ROOT,"--out",area)
            files=[area/"index.csv"]
            need(files[0].is_file(),"Optical index missing")
            for lead in campaign.read(ROOT/"leads.json"):
                oid,h,tag=lead["obsid"],lead["input_sha256"],lead["trigger_key"]
                p=area/"observations"/str(oid)/h/(tag+".json.gz")
                if p.is_file():
                    info=campaign.read_obs(p)
                    need(info["schema"]==1 and (info["lead"]["obsid"],info["lead"]["input_sha256"],info["lead"]["trigger_key"])==(oid,h,tag),
                         "Optical per-trigger identity differs from corrected lead")
                    files.append(p)
                    if info.get("complete") is not True:
                        errors.append(f"{oid} {tag}: partial optical catalog coverage")
                else:
                    errors.append(f"{oid} {tag}: optical record unavailable")
        elif task=="csc":
            files,errors=csc_context(out)
        elif task=="heldout":
            use(HERE/"heldout.py","evaluate")
            area=ROOT/"held_out"
            files=[p for p in area.rglob("*") if p.is_file()]
        elif task=="visits":
            import visits
            area=ROOT/"independent_visits"
            data=visits.run_visits(ROOT,campaign.read(ROOT/"leads.json"),area)
            files=[area/"independent_visits.json",*[area/r["record"] for r in data["nearby_products"]]]
            errors.extend(f'{r["obsid"]}: {r["status"]}' for r in data["nearby_products"] if r["status"]!="complete")
        elif task=="alerts":
            use(HERE/"alerts.py","--evidence",ROOT)
            area=ROOT/"alerts"
            files=[p for p in area.rglob("*") if p.is_file()]
            r=campaign.read(area/"review.json")
            if r["status"] not in {"parsed_archive","no_leads"}:
                errors.append(r["status"])
        elif task=="review":
            use(HERE/"review.py","--evidence",ROOT)
            files=[ROOT/"LEAD_REVIEW.md"]
    except Exception as exc:
        errors.append(f"{type(exc).__name__}: {exc}")
        raise
    finally:
        output_receipt(path,key,files,"partial" if errors else "complete",errors)
    if errors:
        raise ValueError(f"{task} retained partial evidence. Retry this task or document its limits before publication.")


def reports_ready(out):
    measurement=out/"measured_report.json"
    need(repeatable(measurement,provenance()),"The full corrected report is missing or stale")
    for task in ("heldout","csc","optical","visits","alerts","review"):
        path=out/(task+".json")
        need(path.is_file(),"Run final context stage "+task)
        rec=campaign.read(path)
        need(rec["key"]["measurement_sha256"]==sha(measurement) and rec["key"]["leads_sha256"]==sha(ROOT/"leads.json")
             and rec["key"]["program_sha256"]==sha(HERE/{"csc":"csc_variability.py","optical":"lead_optical.py"}.get(task,task+".py"))
             and all((ROOT/p).is_file() and sha(ROOT/p)==h for p,h in rec["outputs"].items()),
             "Context no longer matches corrected primary reports: "+task)
        need(rec["status"]=="complete" or task in {"alerts","visits","csc","optical"} and rec["status"]=="partial",
             "Final statistical or review output is incomplete: "+task)
        if task=="review":
            for p in (ROOT/"independent_visits/independent_visits.json",ROOT/"alerts/review.json"):
                need(rec["key"][p.name]==(sha(p) if p.is_file() else None),"Regenerate human review after new context")
    return {task:campaign.read(out/(task+".json"))["status"] for task in ("heldout","csc","optical","visits","alerts","review")}


def assemble(extra):
    _,out=cycle_path()
    reports_ready(out)
    return final_assembly.assemble(SUP,out,extra)


def publish(extra,retained,bounded_reason):
    import publication
    import candidate_inputs as candidates
    rev,out=cycle_path()
    phases=reports_ready(out)
    residual=campaign.read(rev/"result.json")
    remaining=residual["remaining_obsids"]
    if remaining:
        need(bounded_reason and bounded_reason.strip(),"Retry the residual plan or document a bounded scope with --bounded-reason")
        need(all(residual["changed_or_attempted"][str(i)]["worker_outcome"]!="unattempted" for i in remaining),
             "Some residuals never received a corrective worker; resume the plan before publication")
    decision_file=final_assembly.assemble(SUP,out,extra)
    groups,leads,signature,files=publication.load(SUP)
    decisions,sources=publication.adjudicate(campaign.read(decision_file),leads,signature)
    originals={p:p.read_bytes() if p.is_file() else None for p in
               (HERE.parents[1]/"README.md",HERE/"README.md",ROOT/"publication.json")}
    try:
        pub=publication.publish(groups,leads,decisions,files,sources,signature,decision_file,SUP)
        store=candidates.preserve(decision_file,SUP,retained)
        held=campaign.read(ROOT/"held_out/summary.json")
        lines=["","## Prospective secondary analysis","",
               f"The fixed held-out fields produced {held['evaluable_real_survivors']} real triggers meeting the prespecified criteria. "
               f"The secondary analysis used {held['control_replicas']} scrambled replicas and reported the two-test multiplicity adjustment. "
               "[Protocol](HELD_OUT.md) · [Coverage, survivor scores and control maxima](evidence/held_out/README.md)","",
               f"The exact archive-byte deposit contains {len(store['products'])} manually retained candidate products. "
               "[Transport byte manifest and reassembly commands](evidence/candidate_inputs/manifest.json)."]
        if remaining:
            lines.extend(["",f"{len(remaining)} dated ObsIDs remained incomplete after corrective attempts. {bounded_reason.strip()}"])
        if any(v=="partial" for v in phases.values()):
            limits=", ".join(k for k,v in phases.items() if v=="partial")
            lines.extend(["",f"Supporting checks with incomplete acquisition or coverage include {limits}. Individual review notes preserve these limits."])
        lines.extend(["","[Residual work, measured scope and context receipts]("+out.relative_to(HERE).as_posix()+"/release.json).",""])
        with (HERE/"README.md").open("a") as f:
            f.write("\n".join(lines))
        save(out/"release.json",{"corrected_inputs":provenance(),"publication_sha256":sha(ROOT/"publication.json"),
             "label":pub["label"],"headline":pub["headline"],"campaign_before_and_after":rev.relative_to(ROOT).as_posix(),
             "remaining_obsids":remaining,"bounded_scope_reason":bounded_reason,"context_phases":phases,
             "context_receipts":{k:sha(out/(k+".json")) for k in phases},"manual_review_queue_sha256":sha(out/"review_queue.json"),
             "decisions_sha256":sha(decision_file),"candidate_store_manifest_sha256":sha(candidates.BASE/"manifest.json"),
             "candidate_products":sorted(store["products"]),"survey_readme_sha256":sha(HERE/"README.md")})
        print(pub["headline"])
        if os.getenv("GITHUB_OUTPUT"):
            with open(os.environ["GITHUB_OUTPUT"],"a") as f:
                print("candidate_products="+json.dumps(sorted(store["products"])),file=f)
                print("has_candidates="+str(bool(store["products"])).lower(),file=f)
    except BaseException:
        for p,data in originals.items():
            if data is None:
                p.unlink(missing_ok=True)
            else:
                p.write_bytes(data)
        raise


def candidate_replay(checksum):
    import candidate_inputs as ci
    _,out=cycle_path()
    release=campaign.read(out/"release.json")
    need(checksum in release["candidate_products"] and sha(ci.BASE/"manifest.json")==release["candidate_store_manifest_sha256"],
         "Unknown or changed final candidate input")
    entry=ci.manifest(ci.BASE)["products"][checksum]
    key={"product_sha256":checksum,"input_manifest_sha256":release["candidate_store_manifest_sha256"],
         "implementation_sha256":sha(HERE/"candidate_inputs.py"),"ledger_sha256":sha(pub_local(entry["ledger"]))}
    path=out/"candidate_replays"/(checksum+".json")
    if path.is_file() and campaign.read(path)==key:
        return
    with tempfile.TemporaryDirectory(dir=os.getenv("RUNNER_TEMP")) as temp:
        ci.replay(entry,Path(temp))
    save(path,key)


def pub_local(name):
    import publication
    return publication.local(name)


def verify_commit(receipt):
    import candidate_inputs as ci
    import publication
    repo=HERE.parents[1]
    need(not receipt.resolve().is_relative_to(repo),"Write Git verification receipts to runner temporary storage")
    for command in (["git","diff","--quiet","HEAD","--"],):
        subprocess.run(command,cwd=repo,check=True)
    _,out=cycle_path()
    release=campaign.read(out/"release.json")
    need(release["corrected_inputs"]==provenance() and sha(ROOT/"publication.json")==release["publication_sha256"]
         and sha(HERE/"README.md")==release["survey_readme_sha256"]
         and sha(ci.BASE/"manifest.json")==release["candidate_store_manifest_sha256"],"Final scientific publication changed")
    doc=campaign.read(out/"decisions.json")
    pub=campaign.read(ROOT/"publication.json")
    need(sha(ROOT/"fxt_candidates.csv")==pub["candidate_table_sha256"],"Published candidate table changed")
    need(sha(out/"decisions.json")==release["decisions_sha256"]==pub["decisions_sha256"],"Final written decisions changed")
    sidecars=doc["review_source_sha256"]
    standard={p.relative_to(HERE).as_posix() for p in [*(HERE/n for n in final_assembly.DEFAULTS),
              *(HERE/"manual_decisions").glob("*.json")] if p.is_file()}
    need(standard<=sidecars.keys() and all(sha(publication.local(p))==v for p,v in sidecars.items()),
         "Manual source material changed since assembling the final decisions")
    _,leads,signature,_=publication.load(SUP)
    _,refs=publication.adjudicate(doc,leads,signature)
    need(signature==pub["reviewed_inputs_sha256"] and refs==pub["cited_local_evidence_sha256"],
         "Final measurement or a reviewer’s cited local evidence changed")
    candidate_store=ci.manifest(ci.BASE)
    need(sorted(candidate_store["products"])==release["candidate_products"],"Candidate byte inventory differs from the release")
    for h in release["candidate_products"]:
        replay=campaign.read(out/"candidate_replays"/(h+".json"))
        entry=candidate_store["products"][h]
        need(replay["product_sha256"]==h and replay["input_manifest_sha256"]==release["candidate_store_manifest_sha256"]
             and replay["implementation_sha256"]==sha(HERE/"candidate_inputs.py")
             and replay["ledger_sha256"]==entry["ledger_sha256"]==sha(pub_local(entry["ledger"])),
             "An exact candidate EVT2 replay is still needed")
    commands=[[sys.executable,str(HERE/"candidate_inputs.py"),"check-index"],
              [sys.executable,"-m","unittest","discover","-s",str(HERE),"-p","test_*.py"],
              [sys.executable,"-m","open_space.verify"],[sys.executable,"-m","unittest","discover","-s","tests"]]
    tracked=subprocess.check_output(["git","ls-files","--others","--exclude-standard","--","results/chandra_acis_fxt"],cwd=repo)
    need(not tracked.strip(),"Untracked Chandra results must be committed before the final verifier")
    head=subprocess.check_output(["git","rev-parse","HEAD"],cwd=repo,text=True).strip()
    receipt.parent.mkdir(parents=True,exist_ok=True)
    results=[]
    for i,cmd in enumerate(commands):
        output=receipt.parent/f"final-gate-{i}.log"
        with output.open("wb") as target:
            proc=subprocess.run(cmd,cwd=repo,stdout=target,stderr=subprocess.STDOUT)
        results.append({"command":cmd,"returncode":proc.returncode,"log_sha256":sha(output),"log":output.name})
        save(receipt,{"verified_git_commit":head,"publication_sha256":release["publication_sha256"],"checks":results,
                      "complete":i==len(commands)-1 and all(x["returncode"]==0 for x in results)})
        if proc.returncode:
            with output.open("rb") as f:
                f.seek(max(0,output.stat().st_size-2500))
                print(f.read().decode(errors="replace"),file=sys.stderr)
            raise SystemExit(f"Final gate {i} failed; see {output.name}")


def main():
    need(os.getenv("GITHUB_ACTIONS")=="true" and os.getenv("GITHUB_REPOSITORY")==ingest.REPO and
         os.getenv("GITHUB_REF")=="refs/heads/"+ingest.BRANCH and os.getenv("GITHUB_EVENT_NAME") in {"push","workflow_dispatch"}
         and os.getenv("GITHUB_WORKFLOW_REF","").split("@")[0]==ingest.REPO+"/"+ingest.FLOW,"Use the existing Chandra workflow")
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest="command",required=True)
    for op in ("plan","worker","seal-worker","reduce","report","context","assemble","publish","candidate-replay","verify-commit"):
        c=sub.add_parser(op)
        if op in {"plan","worker","seal-worker","reduce"}:
            c.add_argument("--bundle",required=True,type=Path)
        if op=="plan":
            c.add_argument("--groups",type=int,default=48)
        if op in {"worker","seal-worker"}:
            c.add_argument("--group",type=int,required=True)
            c.add_argument("--stage",type=Path,required=True)
            if op=="worker":
                c.add_argument("--obsids",help="Retry only these comma-separated assigned ObsIDs in a new empty stage")
        if op in {"worker","publish"}:
            c.add_argument("--retained",type=Path,action="append",default=[])
        if op=="reduce":
            c.add_argument("--workers",type=Path,required=True)
            c.add_argument("--checkpoint-partial",action="store_true")
        if op=="context":
            c.add_argument("task",choices=("heldout","csc","optical","visits","alerts","review"))
        if op in {"assemble","publish"}:
            c.add_argument("--sidecar",action="append",default=[])
        if op=="publish":
            c.add_argument("--bounded-reason")
        if op=="candidate-replay":
            c.add_argument("sha256")
        if op=="verify-commit":
            c.add_argument("--receipt",required=True,type=Path)
    a=parser.parse_args()
    if getattr(a,"bundle",None):
        a.bundle=independent(a.bundle)
    if a.command=="plan":
        plan_job(a.bundle,a.groups)
    elif a.command in {"worker","seal-worker"}:
        stage=independent(a.stage)
        need(not stage.is_relative_to(a.bundle) and not a.bundle.is_relative_to(stage),"Separate plan and worker stage")
        if a.command=="worker":
            plan=check_plan(a.bundle)
            need(0<=a.group<len(plan["groups"]),"Invalid worker group")
            if (stage/SEAL).is_file():
                rec=verified(stage,"workers")
                need(rec["plan_sha256"]==sha(a.bundle/"plan.json") and rec["group"]==a.group,"Wrong resumed group")
                if not rec["missing"]:
                    return
            work=[("group",a.group)]
            if a.obsids:
                ids=list(map(int,a.obsids.split(",")))
                need(ids and len(ids)==len(set(ids)) and set(ids)<=set(plan["groups"][a.group]),"Retry outside assigned group")
                work=[("one",oid) for oid in ids]
            for mode,index in work:
                argv=[HERE/"followup.py",mode,ROOT,a.bundle,index,stage]
                for path in a.retained:
                    argv.extend(["--retained",path])
                use(*argv)
        seal_worker(a.bundle,a.group,stage)
    elif a.command=="reduce":
        reduce(a.bundle,independent(a.workers),a.checkpoint_partial)
    elif a.command=="report":
        measured_report()
    elif a.command=="context":
        context(a.task)
    elif a.command=="assemble":
        assemble(a.sidecar)
    elif a.command=="publish":
        publish(a.sidecar,a.retained,a.bounded_reason)
    elif a.command=="candidate-replay":
        candidate_replay(a.sha256)
    else:
        verify_commit(a.receipt)


if __name__=="__main__":
    main()
