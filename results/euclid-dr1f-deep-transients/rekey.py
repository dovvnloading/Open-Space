#!/usr/bin/env python3
"""Repair duplicate display IDs in pre-pixel v2 S3 manifest, NOT the data selection.

All science keys, layer identities, masks, S3 timestamps, independent exposure
sets and chosen observation order remain unchanged; only unique identifiers
and pair-reference keys change. Archives the uncorrected manifest for audit.
"""
import hashlib,json,sys
from pathlib import Path
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P))
import run as R
O=P/"output";M=O/"frozen_manifest.json";DONE=O/"rekey_v3.json"
if DONE.exists():
    print("REKEY_V3_ALREADY_FROZEN",flush=True)
    raise SystemExit(0)
manifest=json.loads(M.read_text())
original_bytes=M.read_bytes()
# Independent of SIA nonresponses; fixed metadata used before this key repair.
if manifest.get("method_version")!="2026-11-19-direct-S3-v2-no-image-inspection":
    raise SystemExit("Unexpected source manifest version: stop rather than silently update.")
if (O/"direct_v2_status.json").exists():
    expected=json.loads((O/"direct_v2_status.json").read_text())["manifest_sha256"]
    if hashlib.sha256(original_bytes).hexdigest()!=expected:
        raise SystemExit("v2 manifest hash changed; cannot do transparent ID repair")
(O/"frozen_manifest_v2_duplicate_id_archive.json").write_bytes(original_bytes)
for v in manifest["visits"]:
    if v["obs_id"].startswith("102158276_"):
        continue
    key=v["products"]["science"]
    unique=hashlib.sha256(key.encode("utf-8")).hexdigest()[:14]
    v["pre_rekey_display_id"]=v["obs_id"]
    v["obs_id"]=f"{v['field']}_{v['band']}_{unique}"
assert len({v["obs_id"] for v in manifest["visits"]})==len(manifest["visits"])
for field in R.FROZEN["fields"]:
    for band in R.BANDS:
        vv=[v for v in manifest["visits"] if v["field"]==field["name"] and v["band"]==band]
        for j in range(len(vv)):
            v2=vv[j]
            v2.pop("pair_independence",None)
            for i in range(j):
                v1=vv[i]
                a,b=v1["layer"],v2["layer"]
                qualifies=(v1["product_association_cleared"] and v2["product_association_cleared"] and
                    a["status"]=="LAYER_POINT_COVERAGE" and b["status"]=="LAYER_POINT_COVERAGE" and
                    bool(a["exposure_ids"]) and bool(b["exposure_ids"]) and
                    set(a["exposure_ids"]).isdisjoint(set(b["exposure_ids"])) and
                    abs(a["median_mjd"]-b["median_mjd"])>R.FROZEN["minimum_independent_time_days"])
                v2.setdefault("pair_independence",{})[v1["obs_id"]]=bool(qualifies)
# Preserve early outputs as provisional if current runner already wrote them;
# no provisional result is a discovery and none is silently mislabelled final.
for p in (O/"summary.json",O/"referee.json",O/"north.json",O/"south.json",
          O/"fornax.json",O/"external_receipts.json"):
    if p.exists():
        p.replace(O/("pre_rekey_"+p.name))
manifest["method_version"]="2026-11-19-direct-S3-v3-unique-object-ids"
manifest["fixed_issue"]="v2 S3 observation display identifiers used a non-unique name segment. v3 adds SHA256(S3 science key) display IDs, no visit/image change."
R.dump(M,manifest)
R.dump(DONE,{"status":"REKEYED_NO_CHANGE_TO_PUBLISHER_SOURCES_OR_SELECTION",
   "v2_original_sha256":hashlib.sha256(original_bytes).hexdigest(),
   "v3_manifest_sha256":hashlib.sha256(M.read_bytes()).hexdigest(),
   "observations":len(manifest["visits"]),
   "eligible_independent_pairs":sum(sum(v.get("pair_independence",{}).values()) for v in manifest["visits"])})
print("REKEY_V3_COMPLETE",len(manifest["visits"]),flush=True)
