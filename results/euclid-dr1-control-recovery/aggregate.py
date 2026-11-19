#!/usr/bin/env python3
"""Receipt-audited aggregation of all selected MER source rows.
Each partition retains its local top 200, hence union contains exact global top
200 (and no claim about omitted lower-ranked individual eligible objects).
"""
from pathlib import Path
import gzip, hashlib, json, math
from collections import Counter
ROOT=Path(__file__).resolve().parent
OUT=ROOT/"output"
PUBLISHED=888152165
def sha(x):return hashlib.sha256(x).hexdigest()
def save(name,x):
    (OUT/name).write_text(json.dumps(x,indent=2,sort_keys=True,allow_nan=False)+"\n")
def key(r):return (-r["rank_score"],r["id"])
def main():
    frozen=json.loads((ROOT/"freeze.json").read_text())
    cal=json.loads((OUT/"calibration.json").read_text())
    if cal["matched_quasar_recovery"]["numerator"]!=21:
        raise RuntimeError("Frozen calibration mismatch")
    observed=set()
    counts=Counter()
    candidates=[]
    for b in range(32):
        manifest=OUT/f"batches/batch_{b:03d}_manifest.json"
        register=OUT/f"batches/batch_{b:03d}_top.jsonl.gz"
        if not manifest.exists() or not register.exists():
            raise RuntimeError("Full 32 batch source scan not yet complete")
        m=json.loads(manifest.read_text())
        if m["status"]!="COMPLETE" or m["batch"]!=b or m["selection_id"]!=frozen["selection_id"]:
            raise RuntimeError("Wrong batch/selection identity")
        if m["index_sha256"]!=frozen["source_index_sha256"] or m["top_sha256"]!=sha(register.read_bytes()):
            raise RuntimeError("Wrong index or corrupt source register in batch "+str(b))
        rows=[]
        with gzip.open(register,"rt",encoding="utf8") as f:
            for line in f:rows.append(json.loads(line))
        if len(rows)!=m["retained_top_rows"]:raise AssertionError("Truncated batch top list")
        allowed={(p["Norder"],p["Npix"]) for p in m["partitions"]}
        by_part=Counter((r["Norder"],r["Npix"]) for r in rows)
        if not set(by_part).issubset(allowed):raise AssertionError("Rows not in manifest partitions")
        if any(n>200 for n in by_part.values()):raise AssertionError("Per-partition rank cap exceeded")
        for p in m["partitions"]:
            ident=(p["Norder"],p["Npix"])
            if ident in observed:raise AssertionError("Partition repeated")
            observed.add(ident)
            if by_part[ident]!=p["retained_top_rows"]:
                raise AssertionError("Per-partition truncated or duplicated top list")
            counts["rows"]+=p["read_rows"]
            counts["eligible"]+=p["eligible_rows"]
            counts["unflagged_eligible"]+=p["eligible_rows"]-p["eligible_flagged_JH"]
            counts["not_optical_gt3"]+=p["eligible_optical_not_gt3"]
            counts["partition_top_saved"]+=p["retained_top_rows"]
        if sum(p["eligible_rows"] for p in m["partitions"])!=m["eligible_rows"]:
            raise AssertionError("Batch selection counts differ")
        candidates+=rows
    if len(observed)!=768 or counts["rows"]!=PUBLISHED:
        raise AssertionError("Source data not exhaustively covered")
    if len(candidates)!=counts["partition_top_saved"]:
        raise AssertionError("Ranked register is incomplete")
    identifiers=[r["id"] for r in candidates]
    if len(set(identifiers))!=len(identifiers):
        raise AssertionError("Duplicate source object ID")
    control=json.loads((OUT/"control_audit.json").read_text())
    excluded={str(c["object_id"]) for c in control if "object_id" in c}
    ranked=sorted((r for r in candidates if r["id"] not in excluded),key=key)
    if len(ranked)<200:raise AssertionError("Too few new photometric rank rows")
    for i,r in enumerate(ranked):
        r["rank"]=i+1
        r["class"]="PHOTOMETRIC_UNVERIFIED_NOT_SPECTROSCOPIC_DISCOVERY"
        r["new_spectroscopic_redshift"]=None
        r["WISE_color_status"]="NOT_YET_VETTED"
    with gzip.open(OUT/"ranked_candidates.jsonl.gz","wt",encoding="utf8",compresslevel=6) as f:
        for r in ranked[:2000]:
            f.write(json.dumps(r,sort_keys=True,allow_nan=False)+"\n")
    save("top_candidates.json",ranked[:100])
    summary={"status":"FULL_HATS_RESELECTION_COMPLETE_UNVERIFIED_PHOTOMETRIC_CANDIDATES",
      "selection_id":frozen["selection_id"],
      "scanned_partitions":len(observed),"scanned_published_catalogue_rows":counts["rows"],
      "eligible_catalogue_rows":counts["eligible"],
      "eligible_rows_without_ground_optical_gt3":counts["not_optical_gt3"],
      "eligible_rows_without_JH_flag":counts["unflagged_eligible"],
      "rank_top_per_partition":200,
      "rank_register_preserved_per_partition":counts["partition_top_saved"],
      "global_rank_top_200_exact_from_partition_tops":True,
      "ranked_without_38_literature_matched_control_ids":len(ranked),
      "quasar_control_recovery":cal["matched_quasar_recovery"],
      "dwarf_control_retention":cal["matched_dwarf_retention"],
      "z_gt7_control_recovery":cal["new_z_gt7_matched_quasars"],
      "source_selection_tuned_on_same_controls":True,
      "footprint_effective_area_deg2":None,"field_completeness_calibrated":False,
      "quasar_purity_calibrated":False,"astrophysical_quasar_posterior":None,
      "new_confirmed_spectroscopic_quasars":0,
      "discovery_level":False,
      "prior_QLF_expectation_complete_1900deg2":"~8-19 z>7 at earlier survey assumptions; not a predicted screened yield"}
    save("summary.json",summary)
    print("FULL_HATS_FINAL",json.dumps(summary),flush=True)
if __name__=="__main__":main()
