#!/usr/bin/env python3
"""Independent result-file referee: deterministic consistency, not an independent astrophysical confirmation."""
import csv, hashlib, json
from pathlib import Path
HERE=Path(__file__).resolve().parent/"output"
def read(name):
    p=HERE/name
    return json.loads(p.read_text()) if p.exists() else None
def check(name, ok, details=""):
    tests.append({"check":name,"pass":bool(ok),"details":details})
tests=[]
r=read("receipts.json"); c=read("checks.json"); p=read("protocol.json")
check("receipts_and_checks_exist",r is not None and c is not None and p is not None)
check("protocol_frozen",bool(r and p and r.get("protocol")==p))
if r:
    cats=r.get("catalogues",[])
    check("official_main_hard_receipts",len(cats)==2 and all(len(x.get("archive_sha256",""))==64 and len(x.get("fits_sha256",""))==64 for x in cats))
    check("input_catalogue_denominators",len(cats)==2 and [x.get("rows") for x in cats]==[930203,5466])
if c:
    x=c.get("counts",{})
    keys=["all_main","quality_finite","flux_det_counts","pointlike","soft","upstream_flags_clean","without_strong_hard_match"]
    check("monotonic_frozen_gates",all(k in x for k in keys) and all(x[keys[i]]>=x[keys[i+1]] for i in range(len(keys)-1)))
    check("followup_no_larger_than_screen",x.get("followup_top_n",1e9)<=min(60,x.get("without_strong_hard_match",0)))
    check("zero_qualified_is_NOT_null",x.get("qualified_novel_XDINS") in (0,None) and any("not" in s.lower() for s in c.get("claims",[])))
if (HERE/"candidates.csv").exists():
    with (HERE/"candidates.csv").open() as f: rows=list(csv.DictReader(f))
    check("candidate_rows_bounded",len(rows)<=60)
    check("no_unvalidated_discovery_status",all(row["status"]!="NEW_XDINS" for row in rows))
    check("candidate_counts_match",bool(c and len(rows)==c.get("counts",{}).get("followup_top_n",-1)))
    check("ranked_deterministically",all(float(rows[i]["ml_counts"])*float(rows[i]["soft_fraction"]) >= float(rows[i+1]["ml_counts"])*float(rows[i+1]["soft_fraction"])-0.04 for i in range(len(rows)-1)))
else: check("candidate_rows_exist",False,"Catalogue stage did not complete; not a null result")
status="passed_output_consistency" if all(t["pass"] for t in tests) else "failed_or_pending"
result={"status":status,"scope":"output-only consistency; not a blind image/spectral/literature referee",
        "checks":tests,"pass":sum(t["pass"] for t in tests),"total":len(tests),
        "not_certified":["completeness","non-variability","optical/IR upper limits","spectral kT","pulsations","novelty"]}
(HERE/"referee.json").write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
print(json.dumps({"referee":status,"pass":result["pass"],"total":len(tests)}))
