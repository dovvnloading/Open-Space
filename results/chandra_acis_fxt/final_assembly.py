"""Combine explicit individual reviews against the corrected primary and supplementary ledgers."""
from collections import defaultdict
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import campaign
import publication as pub

DEFAULTS=("EARLY_DECISIONS.json","PILOT_DECISIONS.json","PRIMARY_EVEN_DECISIONS.json",
          "PRIMARY_ODD_DECISIONS.json","SUPPLEMENT_DECISIONS.json")


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(",",":"),allow_nan=False).encode()).hexdigest()


def event_state(event):
    raw={k:v for k,v in event.items() if k not in {"obsid","trigger_key","catalog_cache","vet","diagnostics","flux"}}
    vet=event.get("vet") or {}
    counterpart={k:v for k,v in vet.items() if k!="queries"}
    counterpart["queries"]={route:{k:q.get(k) for k in ("url","ok","sha256","size_bytes","row_count","query_error")}
                            for route,q in vet.get("queries",{}).items()}
    flux=event.get("flux") or {}
    calibration={k:flux.get(k) for k in ("event_sha256","calibrated_fits_sha256","peak_flux_erg_cm2_s",
                                        "flux_counting_sigma_erg_cm2_s","error")}
    return digest([raw,counterpart,calibration,event.get("diagnostics",{})])


@lru_cache(maxsize=2)
def read_ledger(path):
    return campaign.read_obs(path)


def locate(root,ident):
    oid,h,tag=ident
    path=root/"observations"/f"{oid}.json.gz"
    if not path.is_file():
        return None
    row=read_ledger(path)
    pub.need(row["obsid"]==oid,"Mismatched manual-review observation")
    for product in row["products"]:
        if product["input"]["sha256"]==h:
            i=int(tag.split("/")[1])
            events=product.get("search",{}).get("triggers",[])
            if i<len(events):
                event=events[i]
                pub.need(event.get("trigger_key",tag)==tag,"Inconsistent historical raw trigger")
                return event
    return None


def earliest_primary(root,ident):
    history=root/"history/ned-untyped-v1"
    revisions=sorted((p for p in history.glob("*") if p.is_dir()),key=lambda p:tuple(map(int,p.name.split("-"))))
    for previous in revisions:
        if (previous/"before/observations"/f"{ident[0]}.json.gz").is_file():
            # A first snapshot lacking this raw event means legacy review cannot be authenticated.
            return locate(previous/"before",ident)
    return locate(root,ident)


def choices(ident,reviewed,rows,prior_unchanged):
    """Return a single witnessed manual decision, or a reason for individual re-review."""
    if not rows:
        return None,"missing_manual_review"
    exact=[r for r in rows if r[1].get("review_sha256")==reviewed]
    if exact:
        adopted=exact
    elif any("review_sha256" in row for _,row in rows):
        return None,"manual_review_fingerprint_changed"
    elif not prior_unchanged:
        return None,"reviewed_event_or_catalog_changed"
    else:
        adopted=rows
    dispositions={r["decision"] for _,r in adopted}
    if len(dispositions)!=1:
        return None,"conflicting_current_manual_decisions"
    state=next(iter(dispositions))
    if state=="candidate" and not exact:
        return None,"candidate_needs_final_measurement_review"
    if any(not isinstance(r.get("reason"),str) or not r["reason"].strip() or
           not isinstance(r.get("evidence"),list) or not r["evidence"] for _,r in adopted):
        return None,"review_missing_reason_or_evidence"
    reasons=list(dict.fromkeys(r["reason"].strip() for _,r in adopted))
    evidence=list(dict.fromkeys(x for _,r in adopted for x in r["evidence"]))
    return {"obsid":ident[0],"input_sha256":ident[1],"trigger_key":ident[2],"decision":state,
            "reason":" ".join(reasons),"evidence":evidence,"review_sha256":reviewed,
            "sources":[name for name,_ in adopted]},None


def assemble(supplement,out,extra=()):
    import csv
    root=campaign.EVIDENCE
    groups,leads,sig,inputs=pub.load(supplement)
    files={p for p in (campaign.HERE/name for name in DEFAULTS) if p.is_file()}
    files.update((campaign.HERE/"manual_decisions").glob("*.json"))
    files.update(pub.local(name) for name in extra)
    by_key=defaultdict(list)
    sources={}
    for path in sorted(files):
        name=path.relative_to(campaign.HERE).as_posix()
        sources[name]=pub.file_sha(path)
        doc=campaign.read(path)
        records=doc if isinstance(doc,list) else doc["decisions"]
        seen=set()
        for r in records:
            k=pub.key(r)
            pub.need(k not in seen and r.get("decision") in pub.ALLOWED,"Duplicate or unsupported disposition in "+name)
            seen.add(k)
            by_key[k].append((name,r))
    decisions,requests,not_retained,unrecognized=[],[],[],[]
    optical_index={}
    op_csv=root/"lead_optical/index.csv"
    if op_csv.is_file():
        with op_csv.open(newline="") as stream:
            for row in csv.DictReader(stream):
                optical_index[pub.key(row)]=row
    for k,r in sorted(leads.items()):
        group=r["group"]
        state=event_state(r["event"])
        optical=group["root"]/"lead_optical/observations"/str(k[0])/k[1]/(k[2]+".json.gz")
        optical_hash=None
        if optical.is_file():
            op=campaign.read_obs(optical)
            pub.need(pub.key(op["lead"])==k and op.get("schema")==1,"Optical investigation points to a different trigger")
            optical_hash=pub.file_sha(optical)
        flag=optical_index.get(k,{}).get("gaia_variables_inside_r90")
        flagged=flag not in (None,"","0","0.0","False","false","[]")
        review_hash=digest({"identity":k,"cohort":group["name"],"event_sha256":state,"optical_sha256":optical_hash,
                            "receipt":r["receipt"],"fap":r["fap"],"k":r["n"],"interval":[r["lo"],r["hi"]]})
        prior=earliest_primary(root,k) if group["name"]==pub.PRIMARY else locate(group["root"],k)
        selected,why=choices(k,review_hash,by_key[k],prior is not None and event_state(prior)==state and not flagged)
        if selected:
            try:
                diagnostic=group["root"]/"photometry"/str(k[0])/r["lead"]["input_name"]/k[2]/"diagnostics.json"
                if diagnostic.is_file():
                    recorded=campaign.read(diagnostic)
                    if recorded==r["event"].get("diagnostics") and recorded.get("input_sha256")==k[1]:
                        ref=diagnostic.relative_to(campaign.HERE).as_posix()
                        if ref not in selected["evidence"]:
                            selected["evidence"].append(ref)
                if optical_hash:
                    ref=optical.relative_to(campaign.HERE).as_posix()
                    if ref not in selected["evidence"]:
                        selected["evidence"].append(ref)
                pub.adjudicate({"reviewed_inputs_sha256":sig,"decisions":[selected]},{k:r},sig)
            except (OSError,ValueError,KeyError,TypeError) as exc:
                why="invalid_specific_evidence: "+str(exc)
            else:
                decisions.append(selected)
        requests.append({"obsid":k[0],"input_sha256":k[1],"trigger_key":k[2],"cohort":group["name"],
                         "review_sha256":review_hash,"review_event_sha256":state,"prior_event_sha256":event_state(prior) if prior else None,
                         "catalog_status":r["event"]["vet"]["status"],"optical_record_sha256":optical_hash,
                         "gaia_variables_inside_r90":flag,"survey_fap":r["fap"],"control_exceedances":r["n"],
                         "source_decisions":[{"file":name,"decision":row["decision"]} for name,row in by_key[k]],
                         "status":why or "adjudicated"})
    for k,items in sorted(by_key.items()):
        if k in leads:
            continue
        owners=[g for g in groups if k[0] in g["ids"]]
        raw=locate(owners[0]["root"],k) if len(owners)==1 else None
        if raw is None and len(owners)==1 and owners[0]["name"]==pub.PRIMARY:
            raw=earliest_primary(root,k)
        row={"obsid":k[0],"input_sha256":k[1],"trigger_key":k[2],"sources":[n for n,_ in items],
             "manual_decisions":[d for _,d in items],"recorded_final_or_historical_status":(raw or {}).get("vet",{}).get("status")}
        (not_retained if raw is not None else unrecognized).append(row)
    pending=sum(r["status"]!="adjudicated" for r in requests)+len(unrecognized)
    from ingest_primary import save
    out.mkdir(parents=True,exist_ok=True)
    save(out/"review_queue.json",{"reviewed_inputs_sha256":sig,"review_source_sha256":sources,
         "measurement_files_sha256":inputs,"retained_real_events":requests,"previously_reviewed_not_retained":not_retained,
         "unrecognized_manual_identities":unrecognized,"pending":pending})
    dest=out/"decisions.json"
    if pending:
        # A previously valid release must never be publishable with the new pending queue.
        if dest.exists():
            stale=out/"previous"/(pub.file_sha(dest)+".json")
            stale.parent.mkdir(exist_ok=True)
            dest.replace(stale)
        raise ValueError(f"{pending} individual manual reviews are missing, changed or conflicting; see {out/'review_queue.json'}")
    document={"reviewed_inputs_sha256":sig,"decisions":decisions,"review_source_sha256":sources}
    pub.adjudicate(document,leads,sig)
    save(dest,document)
    return dest
