"""Verify adopted upload-recovery event bytes after canonical repairs, without acquiring any data."""
import argparse
from collections import Counter, defaultdict
import hashlib
from pathlib import Path

import campaign
import followup
import heldout_original
import ingest_primary as ingest
import upload_recovery as recovery


def pins(records):
    return {(r["filename"], r["bytes"], r["sha256"], r["fits_sha256"]) for r in records}


def chain(start, end, changes):
    if start == end:
        return True
    seen, work = {start}, [start]
    while work:
        for nxt in changes.get(work.pop(), ()):
            if nxt == end:
                return True
            if nxt not in seen:
                seen.add(nxt)
                work.append(nxt)
    return False


def audit(root):
    source, coverage = heldout_original.registry(root)
    base = root / recovery.HISTORY
    paths = sorted(base.glob("*/adoption.json"))
    adopted_ids = {oid for p in paths for oid in campaign.read(p)["observations"]}
    links = defaultdict(lambda: defaultdict(set))
    if paths:
        for file in (root / followup.SCOPED).glob("*/result.json"):
            result = campaign.read(file)
            ingest.require(ingest.sha(file.parent/"plan.json") == result["plan_sha256"], "Changed corrective cycle")
            for oid, rec in result["changed_or_attempted"].items():
                if oid not in adopted_ids:
                    continue
                old, new = rec["before_sha256"], rec["after_sha256"]
                if old and new and old != new:
                    backup = file.parent/"before/observations"/f"{oid}.json.gz"
                    ingest.require(backup.is_file() and ingest.sha(backup) == old, "Changed pre-correction ledger")
                    links[oid][old].add(new)
    counts = Counter()
    observed = {}
    groups = {}
    for path in paths:
        m = campaign.read(path)
        g = m["group"]
        roots = [p for p in path.parent.glob("*/transport.json") if ingest.sha(p) == m["transport_sha256"]]
        ingest.require(len(roots) == 1, "Cannot locate unique adopted repeat transport")
        stage = roots[0].parent/"stage"
        stamp = campaign.read(roots[0])
        s = recovery.verify_stage(stage)
        i = s["identity"]
        ingest.require(stamp["seal_sha256"] == m["seal_sha256"] == ingest.sha(stage/"seal.json") and
                       (i["group"], i["repeat_run"], i["repeat_head"]) == (g,m["repeat_run"],m["repeat_code"]) and
                       source["run"]["head_sha"] == i["source_code"] == ingest.CODE and
                       set(s["observations"]) == set(m["observations"]) == set(map(str,i["assigned"])) and
                       stamp["producer_job"]["conclusion"] == "success" and
                       m["original_job"] == i["job_id"], "Adopted provenance differs from repeat transport")
        recovery.adoption_state(root,g,i["job_id"],i["assigned"],source,coverage)
        heldout = root/"held_out"/heldout_original.RECEIPT
        ingest.require(m["original_heldout_sha256"] == ingest.sha(heldout), "Original held-out provenance changed")
        retained = campaign.read(stage/"original/retained.json")
        expected = defaultdict(list)
        for row in retained["products"]:
            expected[row["obsid"]].append(row)
        statuses = Counter()
        for key, earlier in m["observations"].items():
            oid = int(key)
            ingest.require(key not in observed and str(oid) == key, "Duplicate adopted observation")
            p = root/"observations"/f"{oid}.json.gz"
            latest = ingest.ledger(p,oid) if p.is_file() else {"status":"unattempted"}
            initial = s["observations"][key]
            if earlier["selection"] == "worker":
                origin = initial["ledger_sha256"]
            else:
                ingest.require(earlier["selection"] == "checkpoint_after_interrupted_retry",
                               "Unknown recovery selection")
                origin = hashlib.sha256(ingest.blob(ingest.PLAN,ingest.PREFIX+f"/observations/{oid}.json.gz")).hexdigest()
            final_sha = ingest.sha(p) if p.is_file() else None
            ingest.require(chain(origin,final_sha,links[key]), "Canonical ledger lacks identified producer/correction provenance")
            tested = pins(latest.get("acquired", []))
            repeat_pins, witnessed = pins(initial["products"]), pins(expected[oid])
            missing = (repeat_pins | witnessed) - tested
            ingest.require(not missing or latest["status"] not in {"complete","followup_incomplete"},
                           "Completed scan lacks a witnessed/repeated product identity")
            observed[key] = {"group":g,"status":latest["status"],"canonical_sha256":final_sha,"producer_sha256":origin,
                             "original_retained_products":len(witnessed),"repeat_acquired_products":len(repeat_pins),
                             "missing_required_products":[list(t) for t in sorted(missing)]}
            statuses[latest["status"]] += 1
        counts.update(statuses)
        groups[str(g)] = {"original_job":m["original_job"],"original_upload":"missing_artifact",
                          "repeat_run":i["repeat_run"],"repeat_commit":i["repeat_head"],
                          "provider_artifact_id":stamp["artifact"]["id"],"provider_artifact_digest":stamp["artifact"]["digest"],
                          "adoption_sha256":ingest.sha(path),"statuses":dict(statuses)}
    return {"source_run":ingest.RUN,"source_attempt":ingest.ATTEMPT,"source_code":ingest.CODE,
            "assignment_commit":ingest.PLAN,"terminal_original_transport_gaps":coverage["transport_gaps"],
            "groups":groups,"statuses":dict(counts),"observations":observed,
            "scope":"Adopted repeat producers only. Original job transport gaps remain recorded in the primary ingest audit."}


if __name__ == "__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--evidence",type=Path,default=campaign.EVIDENCE)
    args=cli.parse_args()
    root=args.evidence.resolve()
    ingest.save(root/recovery.HISTORY/"audit.json",audit(root))
