"""Conditional recovery for the single primary group-1 product, after its source artifact arrives."""
import argparse
from copy import deepcopy
import gzip
import json
import math
import os
import re
import sys
from pathlib import Path
from urllib.parse import urlencode
from unittest.mock import patch

import campaign
import catalog_pack
import followup as fix
import ingest_primary as ingest
import scan_units as units
import search

OID = 27099
GROUP = "1"
EVENT = ("acisf27099N002_evt2.fits.gz", 127158960,
         "40405115dd4e6df3c7dc5d081ded4674c8a1631929cddeb882554b876a76fd24",
         "158372f2b47dd3eda7bdd67b40207ca1a70717b49f4ecd4c6a5863cda4db3113")
SOURCE_CACHE = f"catalogs/{ingest.RUN}-{ingest.ATTEMPT}-{GROUP}"
need = ingest.require
sha = fix.sha


def save(path, value, packed=False):
    raw = units.encode(value)
    units.atomic(path, gzip.compress(raw, compresslevel=6, mtime=0) if packed else raw)


def source(root):
    receipt = ingest.receipt_file(root,GROUP)
    if receipt.is_file():
        r = campaign.read(receipt)
        need(r["group"] == GROUP and r["artifact_name"] == "chandra-evidence-1" and
             r["plan_commit"] == ingest.PLAN and r["code_commit"] == ingest.CODE and
             set(r["worker_ledger_status"]) == {str(OID)}, "Wrong source run, artifact or assignment")
        path = root / ingest.HISTORY / "members/1.json.gz"
        need(sha(path) == r["member_manifest_sha256"], "Member manifest changed")
        items = campaign.read_obs(path)
        members = {m["source"]:m for m in items}
        need(len(members) == len(items),"Duplicate source member")
        metadata = {"group_receipt_sha256":sha(receipt),"member_manifest_sha256":sha(path),
                    "source_zip":r["zip_digest"],"source_artifact":"verified"}
    else:
        # The normal campaign permits a bounded terminal product when its runner
        # died before upload. Explicitly distinguish this from an active job.
        snap=root/ingest.HISTORY/"source.json"
        report=root/ingest.HISTORY/"coverage.json"
        src,cover=campaign.read(snap),campaign.read(report)
        need(src["run"]["id"] == ingest.RUN and src["run"]["head_sha"] == ingest.CODE and
             src["run"]["status"] == "completed" and cover["source_status"] == src["run"] and
             cover["assignment_commit"] == ingest.PLAN and "chandra-evidence-1" not in src["artifacts"]
             and cover["groups"][GROUP]["artifact_status"] == "missing_artifact"
             and cover["dated_obsids"][str(OID)]["group"] == int(GROUP), "Primary group still pending or unimported")
        members={}
        metadata={"group_receipt_sha256":None,"source_zip":None,"member_manifest_sha256":None,
                  "source_artifact":"missing_at_terminal_run","source_index_sha256":sha(snap),
                  "terminal_coverage_sha256":sha(report)}
    file = fix.old_path(root, OID)
    member = members.get(str(file.relative_to(root)))
    h = sha(file)
    if not member or h != member["sha256"]:
        ingest.immutable(ingest.PLAN)
        baseline=ingest.blob(ingest.PLAN,ingest.PREFIX+"/observations/27099.json.gz")
        need(h == units.digest(baseline), "Expected the immutable checkpoint or imported primary ledger")
    if member:
        archived = fix.safe(root, member["destination"])
        need(sha(archived) == member["sha256"], "Wrong source worker ledger")
    ledger = campaign.read_obs(file)
    need(ledger["obsid"] == OID and len(ledger.get("acquired", [])) == len(ledger["products"]) == 1 and
         ledger["acquired"][0] == ledger["products"][0]["input"], "Expected exactly one acquired ACIS product")
    event = ledger["acquired"][0]
    need(tuple(event[k] for k in ("filename","bytes","sha256","fits_sha256")) == EVENT,
         "Changed ObsID 27099 product revision or bytes")
    environment = members.get("environment-1.json")
    environment_group = GROUP
    if environment is None:
        # Group zero used the identical job definition. It can witness the
        # run's software versions even when the group-1 runner was lost.
        environment_group = "0"
        sibling = campaign.read(ingest.receipt_file(root,"0"))
        index = root/ingest.HISTORY/"members/0.json.gz"
        need(sibling["code_commit"] == ingest.CODE and sibling["plan_commit"] == ingest.PLAN and
             sha(index) == sibling["member_manifest_sha256"],"Missing authenticated source environment")
        [environment] = [f for f in campaign.read_obs(index) if f["source"] == "environment-0.json"]
    ef = fix.safe(root,environment["destination"])
    need(sha(ef) == environment["sha256"], "Changed original environment")
    env = campaign.read(ef)
    need(str(env["run_id"]) == str(ingest.RUN) and env["git_commit"] == ingest.CODE, "Foreign primary environment")
    meta = {**metadata,"obsid":OID, "group":GROUP, "ledger_sha256":sha(file),
            "environment_group":environment_group,"environment_sha256":sha(ef),
            "environment_path":environment["destination"],"event":event,
            "old_catalog_memos":sum(name.startswith(SOURCE_CACHE+"/") and name.endswith(".json") for name in members)}
    return ledger, meta, members


def raw_event(event):
    return {k:v for k,v in event.items() if k not in fix.FOLLOW}


def assignments(product, groups):
    """Assign every as-yet unclassified real/control row; do not select by significance or detector flags."""
    need(1 <= groups <= 32, "Use 1–32 catalog groups")
    rows = {}
    tiles = {}
    for stream, i, event in fix.triggers(product):
        key = f"{stream}/{i}"
        need(key not in rows, "Repeated trigger identity")
        eligible = not bool(event.get("vet", {}).get("queries"))
        row = {"raw_sha256":fix.value_sha(raw_event(event)), "group":None}
        rows[key] = row
        if eligible:
            ra,dec = event.get("ra_deg"),event.get("dec_deg")
            tile = (None,None) if ra is None or dec is None else (round(ra*120),round(dec*120))
            tiles.setdefault(tile, []).append(key)
    load = [0] * groups
    block = max(1, math.ceil(sum(map(len,tiles.values())) / groups))
    for _, keys in sorted(tiles.items(),key=lambda kv:(-len(kv[1]),kv[1][0])):
        for start in range(0,len(keys),block):
            subset = keys[start:start+block]
            g = min(range(groups), key=load.__getitem__)
            for k in subset:
                rows[k]["group"] = g
            load[g] += len(subset)
    return rows


def prepare(root, out, groups, manifest=None, evt=None, capsules=None):
    ledger, src, _ = source(root)
    p = ledger["products"][0]
    state = fix.search_state(p, OID)
    if ledger["status"] == "complete" and not fix.inspect(OID, ledger, root)[0]:
        print("Originally complete. No conditional search or catalog requests.")
        return
    need(state in {"complete","absent"}, "An inconsistent or partial raw checkpoint requires separate review")
    if state == "absent":
        need(all((manifest,evt,capsules)), "Use scan_units plan/work first, then supply manifest, EVT2 and all capsules")
        p = dict(p)
        p["search"] = units.execute(manifest,evt,out/"assembly",staged=capsules,ledger=fix.old_path(root,OID))
        units.atomic(out/"assembly/unit-plan.json",manifest.read_bytes())
    else:
        need(not any((manifest,evt,capsules)), "A full primary raw checkpoint already exists")
    need(fix.search_state(p,OID) == "complete", "Incomplete raw scan or different input/controls/injections")
    fix.source_pair(out)
    fix.checked(root,out)
    primary_env = campaign.read(fix.safe(root,src["environment_path"]))
    env = units.environment()
    packages = {k.lower():v for k,v in primary_env["packages"].items()}
    need(env["python"] == primary_env["python"].split()[0] and env["architecture"] in primary_env["platform"] and
         all(packages.get(k.lower()) == v for k,v in env["packages"].items()) and
         env["omp_threads"] == env["blas_threads"] == "1", "Recovery environment differs from the original worker")
    save(out/"raw.json.gz",p["search"],True)
    rows = assignments(p,groups)
    plan = {"source":src, "primary_scan":state, "raw_sha256":sha(out/"raw.json.gz"),
            "raw_science_sha256":fix.raw_scan(p), "implementation":units.environment(),
            "followup_sha256":sha(Path(fix.__file__)), "reader_sha256":sha(Path(catalog_pack.__file__)),
            "recovery_sha256":sha(Path(__file__)), "new_vet_sha256":fix.NEW, "groups":groups,"rows":rows,
            "assembly_sha256":sha(out/"assembly/assembly.json") if state == "absent" else None}
    save(out/"plan.json",plan)
    print(f"{sum(v['group'] is not None for v in rows.values())} unannotated rows in {groups} groups")


def open_plan(root, prep):
    ledger, src, members = source(root)
    plan = campaign.read(prep/"plan.json")
    fix.checked(root,prep)
    need(plan["source"] == src and plan["implementation"] == units.environment() and
         plan["followup_sha256"] == sha(Path(fix.__file__)) and plan["reader_sha256"] == sha(Path(catalog_pack.__file__))
         and plan["recovery_sha256"] == sha(Path(__file__)) and plan["new_vet_sha256"] == fix.NEW and
         sha(prep/"raw.json.gz") == plan["raw_sha256"], "Recovery source, code or raw scan changed")
    p = {**ledger["products"][0], "search":campaign.read_obs(prep/"raw.json.gz")}
    need(fix.search_state(p,OID) == "complete" and fix.raw_scan(p) == plan["raw_science_sha256"] and
         assignments(p,plan["groups"]) == plan["rows"], "Raw survey coverage, event order or assignments changed")
    if plan["primary_scan"] == "complete":
        need(p["search"] == ledger["products"][0]["search"], "Existing raw search changed")
    else:
        a = campaign.read(prep/"assembly/assembly.json")
        need(sha(prep/"assembly/assembly.json") == plan["assembly_sha256"] and a["event_sha256"] == EVENT[2]
             and a["original_ledger_sha256"] == src["ledger_sha256"] and a["environment"] == plan["implementation"]
             and a["plan_sha256"] == sha(prep/"assembly/unit-plan.json") and
             a["search_json_sha256"] == units.digest(units.encode(p["search"])) and a["unit_files_sha256"],
             "Missing or altered raw-unit assembly")
    return plan,p,members


def old_transport(root, members, read, errors):
    """Reopen original matching successful HTTP bodies, including interrupted catalog work."""
    def cached(label, url, directory, replay, records, params=None):
        full = url + ("?" + urlencode(params) if params else "")
        stem = f"{label}-{units.digest(full.encode())[:24]}.json"
        name = f"{SOURCE_CACHE}/{stem}"
        saved = members.get(name)
        target = directory/stem
        locked = False
        if saved:
            try:
                prior = fix.safe(root,saved["destination"])
                need(sha(prior) == saved["sha256"], "Original HTTP memo changed")
                q = campaign.read(prior)
                need(q["url"] == full, "Original response URL mismatch")
                if not target.exists():
                    for _, _, _, bodies in fix.request_files({"queries":{label:q}},prior.parent):
                        for body, receipt in bodies.items():
                            catalog_pack.recorded_reply(root,body,receipt,directory/body.name)
                    units.atomic(target,prior.read_bytes())
                need(sha(target) == saved["sha256"],"Conflict with pre-existing response memo")
                # Preserve a failed primary and its original successful mirror.
                locked = not q.get("ok")
            except (ValueError,OSError,KeyError,TypeError) as exc:
                errors.append(f"{label}: {exc}")
                raise
        return read(label,url,directory,replay or locked,records,params)
    return cached


def verify_rows(vet, row, cache):
    expected = row["vet"]
    fix.verify_original(expected,cache,cache.parents[1])
    need(fix.memo_matches(cache,expected,cache.parents[1]),"Stored response URL or digest does not match its memo")
    replay = vet.vet_trigger(row["trigger"],cache,replay=True)
    shape = ("status","position","matches","host","errors")
    need(all(replay.get(k) == expected.get(k) for k in shape), "Catalog selection differs on exact HTTP-byte replay")
    pins = {label:(q["url"],q["sha256"]) for label,q in expected["queries"].items() if q.get("ok")}
    fix.require_pins(replay,pins)
    need(fix.certify(replay), "Offline replay did not complete all mandatory catalog searches")


def packet_path(base, key):
    return base/"rows"/(key + ".json.gz")


def producer():
    return {"run":os.environ["GITHUB_RUN_ID"],"attempt":os.environ["GITHUB_RUN_ATTEMPT"],
            "commit":os.environ["GITHUB_SHA"],"job":os.environ["GITHUB_JOB"]}


def provenance(row):
    stamp=row["producer"]
    need(str(stamp["run"]).isdigit() and str(stamp["attempt"]).isdigit() and
         re.fullmatch("[0-9a-f]{40}",stamp["commit"]) and stamp["job"],"Unknown catalog-worker producer")


def catalogs(root, prep, group, out):
    plan,p,members = open_plan(root,prep)
    need(group in range(plan["groups"]), "Unknown catalog group")
    out.mkdir(parents=True,exist_ok=True)
    cache = out/"catalogs"/f"part-{group}"
    cache.mkdir(parents=True,exist_ok=True)
    vet = fix.load_vet(prep/"vet.after.py")
    plan_hash = sha(prep/"plan.json")
    for stream,i,event in fix.triggers(p):
        key = f"{stream}/{i}"
        desc = plan["rows"][key]
        if group != desc["group"]:
            continue
        path = packet_path(out,key)
        identity = {"plan_sha256":plan_hash,"row":key,"group":group,"raw_sha256":desc["raw_sha256"]}
        subject = {**raw_event(event), "obsid":OID, "trigger_key":key}
        if path.exists():
            record = campaign.read_obs(path)
            need(all(record[k] == v for k,v in identity.items()) and record["trigger"] == subject,
                 "Conflicting cached catalog row")
            if record.get("verified"):
                verify_rows(vet,record,cache)
            continue
        record = {**identity,"trigger":subject,"verified":False,"producer":producer()}
        integrity_errors = []
        try:
            with patch.object(vet,"read_input",old_transport(root,members,vet.read_input,integrity_errors)):
                record["vet"] = vet.vet_trigger(subject,cache)
            if fix.certify(record["vet"]):
                verify_rows(vet,record,cache)
                record["verified"] = True
        except (OSError,ValueError,KeyError,TypeError) as exc:
            record["error"] = f"{type(exc).__name__}: {exc}"
        if integrity_errors:
            save(out/"broken-inputs"/(key+".json"),{"identity":identity,"problems":integrity_errors})
            raise ValueError("Cannot complete with altered or missing original catalog responses")
        save(path,record,True)
    save(out/"group.json",{"plan_sha256":plan_hash,"group":group,"producer":producer(),"rows":
         {k:sha(packet_path(out,k)) for k,d in plan["rows"].items() if d["group"] == group},
         "source_files":cache_files(cache)})


def cache_files(cache):
    files = {}
    for path in sorted(cache.iterdir()):
        need(path.is_file() and not path.is_symlink() and path.suffix in {".bin",".json"},
             "Unexpected catalog response file")
        files[path.name] = [path.stat().st_size,sha(path)]
    return files


def load_groups(plan,plan_hash,incoming):
    found, packets = {}, {}
    for name in sorted(incoming.rglob("group.json")):
        row = campaign.read(name)
        provenance(row)
        g = row["group"]
        need(type(g) is int and g in range(plan["groups"]) and g not in found and row["plan_sha256"] == plan_hash,
             "Foreign or duplicate catalog worker")
        found[g] = name.parent
        need(cache_files(name.parent/"catalogs"/f"part-{g}") == row["source_files"],
             "Sealed response file set, size or hash changed")
        expected = {k for k,d in plan["rows"].items() if d["group"] == g}
        need(expected == set(row["rows"]),"Missing or extra catalog row in worker")
        for k,h in row["rows"].items():
            f = packet_path(name.parent,k)
            need(sha(f) == h,"Catalog row changed after worker seal")
            packet = campaign.read_obs(f)
            provenance(packet)
            need(packet["row"] == k and packet["group"] == g and packet["plan_sha256"] == plan_hash and
                 packet["raw_sha256"] == plan["rows"][k]["raw_sha256"] and
                 fix.value_sha(raw_event(packet["trigger"])) == packet["raw_sha256"] and
                 packet["trigger"].get("trigger_key") == k and packet["trigger"].get("obsid") == OID,
                 "Catalog row identity mismatch")
            packets[k] = (name.parent,f)
    required = {d["group"] for d in plan["rows"].values() if d["group"] is not None}
    need(required <= set(found),"Missing catalog work group")
    return packets,found


def finish(root,prep,incoming,repair,out,retained):
    plan,p,_ = open_plan(root,prep)
    plan_hash = sha(prep/"plan.json")
    packets,groups = load_groups(plan,plan_hash,incoming)
    audit = out/"repair_records"/str(OID)/"large-product"
    save(audit/"plan.json",plan)
    if plan["assembly_sha256"]:
        for name in ("assembly.json","unit-plan.json"):
            units.atomic(audit/name,(prep/"assembly"/name).read_bytes())
    recovery = campaign.read(repair/"plan.json")
    need(recovery["repair_sha256"] == sha(Path(fix.__file__)) and
         recovery["catalog_reader_sha256"] == sha(Path(catalog_pack.__file__)) and
         recovery["inventory_sha256"] == fix.INVENTORY,"Changed corrective plan implementation or archive census")
    need(recovery["obsids"][str(OID)]["primary_sha256"] == plan["source"]["ledger_sha256"],
         "Repair plan and source ledger disagree")
    raw_hashes = {key:d["raw_sha256"] for key,d in plan["rows"].items()}
    done = set()
    installed = set()
    def install_responses(g):
        target = out/"catalogs"/fix.namespace()/str(OID)/f"prefetch-{g}"
        if g not in installed:
            base = groups[g]
            seal = campaign.read(base/"group.json")
            for name,(size,checksum) in seal["source_files"].items():
                file = base/"catalogs"/f"part-{g}"/name
                data = file.read_bytes()
                need(len(data) == size and units.digest(data) == checksum,"Changed catalog source bytes")
                units.atomic(target/name,data)
            units.atomic(audit/f"group-{g}.json",(base/"group.json").read_bytes())
            installed.add(g)
        return target
    original = campaign.annotate
    def annotate(trigger,directory,evt,stream,i,cache,output):
        key = f"{stream}/{i}"
        need(evt.name == EVENT[0] and key not in done and
             fix.value_sha(raw_event(trigger)) == raw_hashes[key], "Campaign trigger differs from original raw row")
        done.add(key)
        pins = {}
        locked = set()
        if key in packets:
            base,packet = packets[key]
            row = campaign.read_obs(packet)
            destination = install_responses(row["group"])
            if row.get("vet",{}).get("queries"):
                if row["verified"]:
                    need(fix.certify(row["vet"]) and fix.memo_matches(destination,row["vet"],output),
                         "Prefetched catalogs do not have all verified inputs")
                    fix.verify_original(row["vet"],destination,output)
                    trigger["vet"] = deepcopy(row["vet"])
                    trigger["catalog_cache"] = str(destination.relative_to(output))
                    cache = destination
                else:
                    # Distinct per-row writable cache preserves the original failed reply.
                    cache = output/"catalogs"/fix.namespace()/str(OID)/"large-retry"/stream/str(i)
                    cache.mkdir(parents=True,exist_ok=True)
                    stamp = {"vet":row["vet"],"catalog_cache":f"catalogs/part-{row['group']}"}
                    pins = fix.seed_cache(base,{"obsid":OID},stamp,cache,live=True)
                    prior = row["vet"]
                    locked = {label for label,(mirror,_) in sys.modules["vet"].MIRRORS.items()
                              if prior["queries"].get(label) and not prior["queries"][label].get("ok") and
                              prior["queries"].get(mirror,{}).get("ok") and label not in prior.get("errors",{})}
            units.atomic(packet_path(audit,key),packet_path(base,key).read_bytes())
        if locked:
            module=sys.modules["vet"]
            read=module.read_input
            def pinned(label,url,folder,replay,records,params=None):
                return read(label,url,folder,replay or label in locked,records,params)
            with patch.object(module,"read_input",pinned):
                original(trigger,directory,evt,stream,i,cache,output)
        else:
            original(trigger,directory,evt,stream,i,cache,output)
        if pins:
            fix.require_pins(trigger["vet"],pins)
    def seeded(evt,*,seed,n_controls,n_injections):
        units.verify_evt(evt,p["input"])
        need(plan["primary_scan"] == "absent" and seed == p["search"]["seed"] and
             n_controls == campaign.N_CONTROL and n_injections == campaign.N_INJECT, "Invalid seeded-search call")
        return deepcopy(p["search"])
    with patch.object(campaign,"annotate",annotate), patch.object(search,"search_file",seeded):
        result = fix.run_one(root,repair,recovery,OID,out,retained)
    try:
        for g in groups:
            install_responses(g)
        if result["outcome"] == "complete":
            existing = audit/"completed.json"
            sealed = (existing.is_file() and campaign.read(existing).get("plan_sha256") == plan_hash and
                      campaign.read(existing).get("revised_sha256") == result["revised_sha256"])
            need(done == set(raw_hashes) or (not done and sealed),
                 "Campaign did not annotate all primary real and control rows")
            revised = campaign.read_obs(fix.old_path(out,OID))
            need(fix.raw_scan(revised["products"][0]) == plan["raw_science_sha256"] and
                 not fix.inspect(OID,revised,root,(out,fix.namespace()))[0], "Final raw search or reviewed coverage differs")
        if not (audit/"completed.json").exists():
            save(audit/"completed.json", {"plan_sha256":plan_hash,"revised_sha256":result["revised_sha256"],
                                          "repair_outcome":result["outcome"],"annotated_in_this_invocation":len(done)})
    except (OSError,KeyError,ValueError) as exc:
        path = out/"repair_records"/f"{OID}.json"
        if path.is_file():
            fix.copy(path,audit/"failed-final-validation-worker.json")
            result = {**result,"outcome":"unresolved","revised_sha256":None,"error":str(exc)}
            fix.save(path,result)
        raise
    return result


def main():
    units.repository_guard()
    a = argparse.ArgumentParser(description=__doc__)
    sub=a.add_subparsers(dest="mode",required=True)
    for cmd in ("phase","fetch","prepare","catalogs","finish"):
        s=sub.add_parser(cmd)
        s.add_argument("root",type=Path)
        if cmd == "fetch":
            s.add_argument("output",type=Path)
            s.add_argument("--retained",type=Path,action="append",default=[])
        if cmd in ("prepare","catalogs","finish"):
            s.add_argument("prep",type=Path)
        if cmd == "prepare":
            s.add_argument("--groups",type=int,default=8)
            for k in ("manifest","evt","capsules"):
                s.add_argument("--"+k,type=Path)
        if cmd == "catalogs":
            s.add_argument("group",type=int)
            s.add_argument("output",type=Path)
        if cmd == "finish":
            for k in ("incoming","repair","output"):
                s.add_argument(k,type=Path)
            s.add_argument("--retained",type=Path,action="append",default=[])
    args=a.parse_args()
    root=args.root.resolve()
    independent = [root]
    for k in ("prep","incoming","repair","output"):
        if hasattr(args,k):
            path=getattr(args,k).resolve()
            need(not any(path.is_relative_to(prev) or prev.is_relative_to(path) for prev in independent),
                 "Original evidence, prepared plan, corrective bundle and worker stages require disjoint directories")
            independent.append(path)
    if args.mode == "finish":
        need(os.getenv("CIAO_PREFIX") and os.getenv("RETAIN_EVENTS"),"Set the normal photometry and retention environment")
    if args.mode == "phase":
        row,meta,_=source(root)
        p=row["products"][0]
        print(json.dumps({"status":row["status"],"raw_scan":fix.search_state(p,OID),
                          "raw_rows":sum(1 for _ in fix.triggers(p)),"catalog_memos":meta["old_catalog_memos"],
                          "remaining":fix.inspect(OID,row,root)[0],"source":meta},sort_keys=True))
    elif args.mode == "fetch":
        row,_,_=source(root)
        args.output.mkdir(parents=True,exist_ok=True)
        attempts=[]
        try:
            fix.event_bytes(row["acquired"][0],OID,args.output,args.retained,attempts)
        finally:
            fix.save(args.output/"acquisition.json",attempts)
    elif args.mode == "prepare":
        prepare(root,args.prep.resolve(),args.groups,args.manifest,args.evt,args.capsules)
    elif args.mode == "catalogs":
        catalogs(root,args.prep.resolve(),args.group,args.output.resolve())
    else:
        finish(root,args.prep.resolve(),args.incoming.resolve(),args.repair.resolve(),args.output.resolve(),args.retained)


if __name__ == "__main__":
    main()
