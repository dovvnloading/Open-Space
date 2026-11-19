#!/usr/bin/env python3
"""Network-free, independent validation of saved visit and area receipts."""
import hashlib,json
from pathlib import Path
P=Path(__file__).resolve().parent
O=P/"output"
def load(p):return json.loads(p.read_text())
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
f=load(P/"freeze.json");m=load(O/"frozen_manifest.json");s=load(O/"summary.json")
assert m["freeze_file_sha256"]==digest(P/"freeze.json")
assert {x["name"] for x in f["fields"]}=={"north","south","fornax"}
assert f["patch_pixels"]==1200
chosen={}
for v in m["visits"]:
    assert v["band"] in ("VIS","Y","J","H")
    assert v["field"] in ("north","south","fornax")
    assert v["obs_id"] not in chosen
    chosen[v["obs_id"]]=v
    assert v["products"]["science"].startswith("dr1/MER/")
statuses=[];raw_count=0
for field in f["fields"]:
    p=load(O/(field["name"]+".json"))
    assert p["field"]==field["name"]
    assert 0<=p["nominal_square_deg"]<=.001112
    raw_count+=len(p["raw_triggers"])
    for b in p["band_results"]:
        assert b["band"] in ("VIS","Y","J","H")
        if b["status"]=="EMPIRICALLY_SCREENED_UNCONFIRMED":
            v1,v2=[chosen[k] for k in b["obs_ids"]]
            l1,l2=v1["layer"],v2["layer"]
            assert l1["status"]=="LAYER_POINT_COVERAGE"
            assert l2["status"]=="LAYER_POINT_COVERAGE"
            assert set(l1["exposure_ids"]).isdisjoint(set(l2["exposure_ids"]))
            assert abs(l1["median_mjd"]-l2["median_mjd"])>1
            assert 0<b["detection_area_deg2"]<=b["joint_area_deg2"]<=p["nominal_square_deg"]*1.001
            assert b["raw_trigger_count"]==len([x for x in p["raw_triggers"]
                if x["band"]==b["band"] and x["pair_obs_ids"]==b["obs_ids"]])
            assert b["calibration"]["status"]=="EMPIRICAL_COMPACT_SOURCE_PSF_AND_FLUX_PROXY"
        statuses.append(b["status"])
assert s["independent_calibrated_visit_pairs"]==statuses.count("EMPIRICALLY_SCREENED_UNCONFIRMED")
assert s["raw_difference_triggers"]==raw_count
assert s["validated_astrophysical_candidates"]==0
report={"status":"BOUNDED_REPLAY_PASS_NO_VALIDATED_DISCOVERY" if s["independent_calibrated_visit_pairs"]
    else "INCONCLUSIVE_RECEIPT_REPLAY_PASS_NO_SCREENED_PAIR",
    "checked_selected_visits":len(chosen),
    "verified_difference_pairs":s["independent_calibrated_visit_pairs"],
    "raw_unvetted_triggers":raw_count,"rate_limit_validated":False,"remarkable_threshold_met":False}
(O/"referee.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
print(json.dumps(report),flush=True)
