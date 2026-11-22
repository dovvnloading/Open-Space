"""Evaluate one protocol fixed before opening the second half of the primary shards."""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
from unittest.mock import patch

import campaign
import early
import report
from archive import dump
from migrate import INVENTORY, verify_scan
from search import DURATIONS, _time_support

COMMIT = "0655e5b54215a3d7bbac1f6b2561379bb50c81cf"
PROTOCOL = "7c70f27929c62a1fb462485e20d7f25312ebf4eb97c132f0708faefcf86fa2f2"
PREFIX = "results/chandra_acis_fxt/"
KNOWN = {"no_counterpart","host_candidate","ambiguous","foreground_veto","known_transient"}
ACCEPT = {"no_counterpart","host_candidate"}


def sha(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f,"sha256").hexdigest()


def require(test, why):
    if not test:
        raise ValueError(why)


def git(*args):
    return subprocess.check_output(["git","-C",str(campaign.HERE.parents[1]),*args])


def freeze(root, dest):
    require(sha(campaign.HERE/"HELD_OUT.md")==PROTOCOL,"Prospectively fixed protocol changed")
    exists = subprocess.run(["git","-C",str(campaign.HERE.parents[1]),"cat-file","-e",COMMIT+"^{commit}"],
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if exists.returncode:
        git("fetch","--no-tags","--depth=1","--filter=blob:none","origin",COMMIT)
    require(git("rev-parse",COMMIT+"^{commit}").decode().strip()==COMMIT,"Assignment commit changed")
    folder=dest/"assignment_sources"
    folder.mkdir(parents=True, exist_ok=True)
    sources={}
    def blob(name, dest_file=None):
        oid=git("rev-parse",COMMIT+":"+PREFIX+name).decode().strip()
        raw=git("cat-file","blob",oid)
        if dest_file:
            (folder/dest_file).write_bytes(raw)
        sources[name]={"git_blob":oid,"sha256":hashlib.sha256(raw).hexdigest(),"bytes":len(raw)}
        return raw
    assignments=json.loads(blob("evidence/inventory/shards.json","original_shards.json"))
    inv=json.loads(blob("evidence/inventory/inventory.json"))
    require(sources["evidence/inventory/inventory.json"]["sha256"]==INVENTORY==sha(root/"inventory/inventory.json"),
            "The acquired campaign has a different parent inventory")
    blob("early.py","original_early.py.txt")
    require(sources["early.py"]["sha256"]==sha(campaign.HERE/"early.py"),"Early selector differs from the fixed revision")
    with patch.object(campaign,"EVIDENCE",root):
        early_selection=early.select()
    earlier={i for group in early_selection["groups"] for i in group}
    dump(folder/"early_selection.json",early_selection)
    pilot=json.loads(blob("evidence/pilot.json","original_pilot.json"))
    pilot_ids={r["obsid"] for r in pilot["selected"]}
    migration=json.loads(blob("evidence/history/44981778618/migration.json","original_migration.json"))
    migrated={r["obsid"] for r in migration["ledgers"]}
    require(len(assignments)==192 and len(earlier)==128 and len(pilot_ids)==3 and migrated,
            "Unexpected fixed cohort dimensions")
    all_ids=[oid for group in assignments for oid in group]
    population={r["obsid"] for r in inv}
    require(len(all_ids)==len(set(all_ids)) and set(all_ids) <= population and earlier|pilot_ids|migrated <= population,
            "Assignment or exclusion outside frozen inventory")
    exclusions={"early":earlier,"pilot":pilot_ids,"source_migration":migrated}
    selected=[{"group":g,"obsid":oid,"exclusions":[k for k,ids in exclusions.items() if oid in ids]}
              for g in range(96,192) for oid in assignments[g]]
    plan={"commit":COMMIT,"protocol_sha256":PROTOCOL,"frozen_utc":"2026-11-22T06:35:00Z","sources":sources,
          "excluded_obsids":{k:sorted(v) for k,v in exclusions.items()},"assigned":selected}
    dump(dest/"selection.json",plan)
    return selected


def num(v):
    return type(v) in (int,float) and math.isfinite(v)


def criteria(t):
    """A definitive failure dominates a missing different predicate in every stream."""
    failed, missing=[],[]
    instrument=t.get("instrument_pass")
    if instrument is False:
        failed.append("instrument")
    elif instrument is not True:
        missing.append("instrument")
    vet=t.get("vet") or {}
    status=vet.get("status")
    if status not in KNOWN or vet.get("errors") != {}:
        missing.append("catalog")
    elif status not in ACCEPT:
        failed.append("catalog")
    b=(vet.get("position") or {}).get("galactic_latitude_deg")
    if not num(b) or abs(b)>90:
        missing.append("latitude")
    elif abs(b)<30:
        failed.append("latitude")
    duration, before, after, bt, at=[t.get(k) for k in
        ("trial_duration_s","before_counts","after_counts","before_s","after_s")]
    expected=None
    if (num(bt) and bt<=0) or (num(at) and at<=0):
        failed.append("flank_exposure")
    elif not (num(duration) and duration>0 and num(before) and num(after) and min(before,after)>=0
              and num(bt) and num(at)):
        missing.append("flank_rate")
    else:
        expected=duration*(before+after)/(bt+at)
        if not math.isfinite(expected):
            missing.append("flank_rate")
        elif expected>1:
            failed.append("flank_rate")
    if not failed and not missing and not all(num(t.get(k)) for k in
         ("seed_score","ra_deg","dec_deg","start_met_s","stop_met_s","position_error_arcsec")):
        missing.append("ranking_or_deduplication")
    return ("fail" if failed else "unknown" if missing else "pass"),failed,missing,b,expected


FIELDS=("group","obsid","input_sha256","trigger_key","instrument_pass","primary_status","galactic_latitude_deg",
        "trial_duration_s","before_counts","after_counts","before_s","after_s","expected_flank_photons",
        "seed_score","secondary_status","failed","missing","real_photometry_missing")


def examine(assignment,root,writer):
    oid=assignment["obsid"]
    path=root/"observations"/f"{oid}.json.gz"
    evidence={**assignment,"ledger_sha256":None,"errors":[],"products":[],
              "decision_counts":Counter(),"real_photometry_missing":0}
    streams=defaultdict(list)
    support=defaultdict(list)
    if assignment["exclusions"]:
        evidence["disposition"]="previously_exposed"
        return evidence,streams,support
    if not path.is_file():
        evidence["disposition"]="unattempted"
        return evidence,streams,support
    evidence["ledger_sha256"]=sha(path)
    try:
        obs=campaign.read_obs(path)
        evidence["primary_status"]=obs["status"]
        evidence["primary_errors"]=obs["errors"]
        require(obs["obsid"]==oid,"Primary ledger identifies a different ObsID")
        acquired={(r["filename"],r["sha256"]) for r in obs.get("acquired",[])}
        found=[(p["input"]["filename"],p["input"]["sha256"]) for p in obs.get("products",[])]
        if not acquired or acquired!=set(found) or len(found)!=len(set(found)):
            evidence["errors"].append("missing_or_partial_product_inventory")
        for p in obs.get("products",[]):
            source=p["input"]
            evidence["products"].append({k:source.get(k) for k in ("filename","bytes","sha256","fits_sha256","url")})
            scan=p.get("search")
            if not scan:
                evidence["errors"].append(source["filename"]+" missing_primary_search")
                continue
            try:
                verify_scan(scan,source,oid)
                for d in DURATIONS:
                    support[d].extend(scan["start_support_by_duration"][str(d)])
            except (ValueError,TypeError,KeyError) as exc:
                evidence["errors"].append(source["filename"]+" "+str(exc))
            trials=[(-1,scan.get("triggers"))]+[(c.get("replica"),c.get("triggers")) for c in scan.get("controls",[])]
            for stream, events in trials:
                if not isinstance(events,list) or type(stream) is not int or not -1<=stream<24:
                    evidence["errors"].append(source["filename"]+" malformed_trigger_stream")
                    continue
                for index,t in enumerate(events):
                    decision,fails,gaps,b,expected=criteria(t)
                    key=("real" if stream==-1 else f"control-{stream}")+"/"+str(index)
                    if t.get("trigger_key") not in {None,key}:
                        evidence["errors"].append("trigger_key_disagrees")
                    vet=t.get("vet") or {}
                    absent=stream==-1 and t.get("instrument_pass") and vet.get("status") not in report.FINAL_VETO and (
                        not num((t.get("flux") or {}).get("peak_flux_erg_cm2_s")) or not (t.get("diagnostics") or {}).get("files"))
                    evidence["real_photometry_missing"]+=int(bool(absent))
                    evidence["decision_counts"][("real_" if stream==-1 else "control_")+decision]+=1
                    writer.writerow({**t,"group":assignment["group"],"obsid":oid,"input_sha256":source["sha256"],
                        "trigger_key":key,"primary_status":vet.get("status"),"galactic_latitude_deg":b,
                        "expected_flank_photons":expected,"secondary_status":decision,"failed":"|".join(fails),
                        "missing":"|".join(gaps),"real_photometry_missing":bool(absent)})
                    if decision=="pass":
                        fields=("seed_score","ra_deg","dec_deg","position_error_arcsec","start_met_s","stop_met_s",
                                "duration_s","trial_duration_s","event_mjd","before_counts","after_counts","before_s","after_s")
                        streams[stream].append({**{k:t.get(k) for k in fields},"obsid":oid,"trigger_key":key,
                            "input_sha256":source["sha256"],"input_url":source["url"],"input_name":source["filename"],
                            "input_bytes":source["bytes"],"primary_status":vet["status"],"host":vet.get("host"),
                            "galactic_latitude_deg":b,"expected_flank_photons":expected,
                            "peak_flux_erg_cm2_s":(t.get("flux") or {}).get("peak_flux_erg_cm2_s"),
                            "counting_sigma":(t.get("flux") or {}).get("flux_counting_sigma_erg_cm2_s")})
    except (OSError,EOFError,ValueError,KeyError,TypeError) as exc:
        evidence["errors"].append(f"{type(exc).__name__}: {exc}")
    ambiguous=sum(evidence["decision_counts"][s+"_unknown"] for s in ("real","control"))
    evidence["disposition"]="incomplete" if evidence["errors"] or ambiguous else "evaluable"
    return evidence,streams,support


def fap(score,maxima):
    require(len(maxima)==24,"All 24 held-out control maxima are required")
    k=sum(m is not None and m>=score for m in maxima)
    p=(k+1)/25
    return {"k":k,"replicas":24,"exceedance_fraction":k/24,"survey_fap":p,
            "exact_binomial_95":report.bounds(k,24),"two_test_bonferroni":min(1.,2*p)}


def evaluate(selection,root,dest):
    coverage,real=[],[]
    controls=[{"replica":i,"distinct_survivors":0,"maximum_score":None} for i in range(24)]
    exposure=Counter()
    temp=dest/"trigger_audit.csv.gz.tmp"
    with temp.open("wb") as raw, gzip.GzipFile(filename="",fileobj=raw,mode="wb",mtime=0) as gz:
        import io
        with io.TextIOWrapper(gz,encoding="utf-8",newline="") as fh:
            writer=csv.DictWriter(fh,fieldnames=FIELDS,extrasaction="ignore")
            writer.writeheader()
            for item in selection:
                record,streams,support=examine(item,root,writer)
                coverage.append(record)
                if record["disposition"]!="evaluable":
                    continue
                for d,spans in support.items():
                    exposure[d]+=sum((b-a)*duty for a,b,duty in _time_support(spans))
                real.extend(report.distinct(streams[-1]))
                for i,control in enumerate(controls):
                    survivors=report.distinct(streams[i])
                    control["distinct_survivors"]+=len(survivors)
                    if survivors:
                        top=max(t["seed_score"] for t in survivors)
                        m=control["maximum_score"]
                        control["maximum_score"]=top if m is None else max(m,top)
    temp.replace(dest/"trigger_audit.csv.gz")
    maxima=[c["maximum_score"] for c in controls]
    real.sort(key=lambda r:(-r["seed_score"],r["obsid"],r["trigger_key"]))
    for candidate in real:
        candidate.update(fap(candidate["seed_score"],maxima))
    tally=Counter(r["disposition"] for r in coverage)
    dump(dest/"coverage.json",coverage)
    dump(dest/"controls.json",controls)
    dump(dest/"survivors.json",real)
    dump(dest/"summary.json",{"assignment_commit":COMMIT,"protocol_sha256":PROTOCOL,
         "selection_sha256":sha(dest/"selection.json"),"implementation_sha256":sha(Path(__file__)),
         "primary_implementation":campaign.read(root/"inventory/campaign.json")["implementation"],
         "assigned_obsids":len(selection),"obsid_dispositions":tally,"evaluable_real_survivors":len(real),
         "eligible_start_support_s":{str(d):exposure[d] for d in DURATIONS},"control_replicas":24,
         "control_survivors":sum(c["distinct_survivors"] for c in controls),"minimum_adjusted_fap":0.08,
         "artifacts":{k:sha(dest/k) for k in ("coverage.json","controls.json","survivors.json","trigger_audit.csv.gz")},
         "run_id":os.environ["GITHUB_RUN_ID"],"run_attempt":os.environ["GITHUB_RUN_ATTEMPT"]})
    lines=["# Prospective held-out secondary analysis","",
           f"**{len(real)} distinct real triggers met the frozen rule in {tally['evaluable']:,} evaluable ObsIDs.**","",
           f"Groups 96–191 assigned {len(selection):,} ObsIDs. Prior exposure excluded {tally['previously_exposed']:,}. "
           f"{tally['unattempted']:,} lacked ledgers and {tally['incomplete']:,} had incomplete statistical coverage. "
           f"The 24 equal-scope control catalogs retained {sum(c['distinct_survivors'] for c in controls)} triggers in total. "
           "The minimum two-test-adjusted FAP is 0.08. Inference covers only the evaluable held-out observations.","",
           "| ObsID | Trigger | Seed score | k/24 | Survey FAP | Exact 95% interval for exceedance | Two-test adjusted |",
           "|---:|---|---:|---:|---:|---|---:|"]
    for t in real:
        lo,hi=t["exact_binomial_95"]
        lines.append(f"| {t['obsid']} | {t['trigger_key']} | {t['seed_score']:.6g} | {t['k']}/24 | "
                     f"{t['survey_fap']:.3f} | {lo:.4g}–{hi:.4g} | {t['two_test_bonferroni']:.3f} |")
    lines.extend(["","[Prospective protocol](../../HELD_OUT.md) · [Immutable sample](selection.json) · "
                 "[Coverage and exclusions](coverage.json) · [Every evaluated raw classification](trigger_audit.csv.gz) · "
                 "[Control maxima](controls.json) · [Positions, photometry, input receipts and scores](survivors.json)",""])
    (dest/"README.md").write_text("\n".join(lines))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode",choices=["freeze","evaluate"])
    args=p.parse_args()
    require(os.getenv("GITHUB_ACTIONS")=="true" and os.getenv("GITHUB_REPOSITORY")=="dovvnloading/Open-Space"
            and os.getenv("GITHUB_REF")=="refs/heads/research/chandra-acis-fast-transients",
            "Use the Chandra research workflow")
    root=campaign.EVIDENCE
    dest=root/"held_out"
    selected=freeze(root,dest)
    if args.mode=="evaluate":
        evaluate(selected,root,dest)
