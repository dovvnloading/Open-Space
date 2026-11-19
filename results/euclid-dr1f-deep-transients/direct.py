#!/usr/bin/env python3
"""Second, prospectively pinned acquisition route after IRSA SIA 502.

Uses ONLY the already deposited metadata S3 probe, the older North SIA
response, and before-pixel deterministic selections. Treats uncertain
science/RMS/flag layering as an uncalibrated eligibility failure.
"""
import csv,hashlib,json,re,sys
from datetime import datetime,timezone
from pathlib import Path
from xml.etree import ElementTree as ET
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P))
import run as R
O=P/"output"
PROBE=O/"s3_probe.json"
OLD=P.parent/"euclid-foundation-eligibility"/"inputs"/"sia-vis.xml"
MARKER=O/"direct_v2_status.json"
def dt(key):
    m=re.search(r"_(\d{8}T\d{6})[.0-9]*Z_",key.rsplit("/",1)[-1])
    if not m: raise ValueError("Publisher filename has no release timestamp: "+key)
    return datetime.strptime(m.group(1),"%Y%m%dT%H%M%S").replace(tzinfo=timezone.utc)
def original_sia_vis():
    doc=ET.parse(OLD).getroot()
    table=doc.find(".//{*}TABLE")
    cols=[f.attrib["name"] for f in table.findall("./{*}FIELD")]
    result={}
    for tr in table.findall(".//{*}TR"):
        r=dict(zip(cols,[(c.text or "").strip() for c in tr.findall("./{*}TD")]))
        obs=r.get("obs_id","")
        if obs not in ("102158276_013663548_VIS","102158276_028106832_VIS",
            "102158276_343226670_VIS"):continue
        url=r.get("access_url","")
        try:key=R.key_from_url(url)
        except ValueError:continue
        stem=key.rsplit("/",1)[-1]
        kind=next((name for name,token in (("science","BGSUB-MOSAIC"),
            ("rms","-RMS_"),("flag","-FLAG_"),("layer","DETECTOR-LAYERING"))
            if token in stem),None)
        if kind:result.setdefault(obs,{})[kind]=key
    return result
def merge_files(tile,band):
    names={}
    for kind in ("science","rms","flag","layer"):
        names[kind]=list(tile.get("files",{}).get(band+"_"+kind,[]))
    return names
def nearest_association(science,files):
    out={"science":science};proof={"selection":"FIRST_LEXICAL_PUBLISHED_S3_SCIENCE_KEY",
        "science_s3_key":science,"product_mapping":"NEAREST_PUBLISHER_TIMESTAMP_NOT_A_SIA_JOIN",
        "separation_hours":{},"second_nearest_margin_hours":{}}
    t=dt(science)
    for kind in ("rms","flag","layer"):
        opts=sorted(((abs((dt(k)-t).total_seconds())/3600,k) for k in files[kind]))
        if not opts:continue
        distance,k=opts[0]
        out[kind]=k
        proof["separation_hours"][kind]=round(distance,3)
        proof["second_nearest_margin_hours"][kind]=round(opts[1][0]-distance,3) if len(opts)>1 else None
    # At least 20-hour separation from every competing mosaic for each product.
    # A single instrument exposure per tile is unambiguous (even if a new
    # cross-tile overlap still must pass the detector-footprint test).
    separation=proof["second_nearest_margin_hours"]
    unambiguous=all(kind in out and
        (separation.get(kind) is None or separation[kind]>=20)
        for kind in ("rms","flag","layer"))
    proof["association_provenance"]=("UNAMBIGUOUS_DISTINCT_PUBLICATION_EPOCHS"
        if unambiguous else "AMBIGUOUS_ASSOCIATION_NOT_USABLE_FOR_PHOTOMETRY")
    return out,proof
def main():
    if MARKER.exists():
        print("DIRECT_V2_MANIFEST_ALREADY_FROZEN",flush=True);return
    probe=json.loads(PROBE.read_text())
    if probe["status"]!="S3_PROBE_DONE":raise RuntimeError("S3 probe not complete")
    if not OLD.exists():raise RuntimeError("North prior SIA receipts not in checkout")
    legacy=original_sia_vis()
    fields={f["field"]:f for f in probe["fields"]}
    manifest={"method_version":"2026-11-19-direct-S3-v2-no-image-inspection",
        "predecessor":"failed IRSA SIA 502 for all three fields; retained in output/sia_v1_failure.json",
        "freeze_file_sha256":R.sha((P/"freeze.json").read_bytes()),
        "probe_file_sha256":R.sha(PROBE.read_bytes()),"created_utc":R.utc(),
        "release":"IRSA DR1-Foundation (2026-11-12), public deep HATS/S3",
        "fields":R.FROZEN["fields"],"visits":[],"receipts":[]}
    # Freeze from the metadata probe, not transient detection counts. No data
    # pixels have yet been read by the v2 workflow.
    for field in R.FROZEN["fields"]:
        name=field["name"]
        p=fields[name]
        if p["status"]!="S3_DISCOVERY_PROBE_COMPLETE":
            manifest["receipts"].append({"field":name,"status":"S3_METADATA_UNAVAILABLE",
                "failure":p.get("error")});continue
        tiles=p.get("tiles",[])
        if name=="north":
            tiles=[t for t in tiles if t["tileid"]=="102158276"]
        if name=="fornax":
            tiles=[t for t in tiles if t["tileid"]=="101375805"]
        choice={"field":name,"status":"S3_PUBLISHER_METADATA_SELECTED",
            "hats_sources_near_center":p.get("hats_sources_within_3arcmin"),
            "hats_tile_counts":p.get("tile_counts"),"tileids":[t["tileid"] for t in tiles],
            "selected_observation_keys":[],"selection":"pre-pixel S3 prefix and prior North exact SIA IDs"}
        manifest["receipts"].append(choice)
        for band in R.BANDS:
            selected=[]
            if name=="north" and band=="VIS":
                for obs in ("102158276_013663548_VIS","102158276_028106832_VIS",
                    "102158276_343226670_VIS"):
                    pths=legacy.get(obs,{})
                    if "science" not in pths:continue
                    selected.append((obs,pths,{"association_provenance":"EXACT_OLDER_SIA_OBS_ID",
                        "selection":"FROZEN_OLD_NORTH_VIS_IDS","source_sia_xml":str(OLD.relative_to(P.parent.parent))}))
            else:
                for tile in tiles:
                    sources=merge_files(tile,band)
                    for key in sorted(sources["science"]):
                        products,provenance=nearest_association(key,sources)
                        obs=f"{name}_{band}_{tile['tileid']}_{key.rsplit('/',1)[-1].split('_')[1]}"
                        selected.append((obs,products,provenance))
                # Sort by science key; only first max_visits_per_band eligible
                selected=sorted(selected,key=lambda x:x[1]["science"])[:R.FROZEN["max_visits_per_band"]]
            for obs,paths,provenance in selected:
                visit={"obs_id":obs,"band":band,"field":name,
                    "center":[field["ra_deg"],field["dec_deg"]],
                    "products":paths,"provenance":provenance}
                # An unambiguous date group or exact SIA association is needed
                # before a pixel pair may be counted as calibrated.
                visit["product_association_cleared"]=provenance["association_provenance"] in (
                    "EXACT_OLDER_SIA_OBS_ID","UNAMBIGUOUS_DISTINCT_PUBLICATION_EPOCHS")
                layer=R.layer_receipt(visit,*visit["center"])
                visit["layer"]=layer
                manifest["visits"].append(visit)
                choice["selected_observation_keys"].append({"obs_id":obs,"band":band,
                    "science":paths["science"],"association":provenance["association_provenance"],
                    "layer_status":layer["status"],"mjd":layer.get("median_mjd")})
                print("S3_SELECTED",name,band,obs,"assoc",visit["product_association_cleared"],
                      "layer",layer["status"],flush=True)
    for field in R.FROZEN["fields"]:
        for band in R.BANDS:
            vv=[v for v in manifest["visits"] if v["field"]==field["name"] and v["band"]==band]
            for j in range(len(vv)):
                for i in range(j):
                    a,b=vv[i],vv[j]
                    la,lb=a["layer"],b["layer"]
                    independence=(a["product_association_cleared"] and b["product_association_cleared"] and
                        la["status"]=="LAYER_POINT_COVERAGE" and lb["status"]=="LAYER_POINT_COVERAGE" and
                        bool(la["exposure_ids"]) and bool(lb["exposure_ids"]) and
                        set(la["exposure_ids"]).isdisjoint(set(lb["exposure_ids"])) and
                        abs(la["median_mjd"]-lb["median_mjd"])>1.0)
                    b.setdefault("pair_independence",{})[a["obs_id"]]=bool(independence)
    old=O/"frozen_manifest.json"
    if old.exists() and not (O/"sia_v1_failure.json").exists():
        (O/"sia_v1_failure.json").write_bytes(old.read_bytes())
    for path in [old,O/"summary.json",O/"referee.json",
                 O/"north.json",O/"south.json",O/"fornax.json",
                 O/"external_receipts.json"]:
        path.unlink(missing_ok=True)
    R.dump(old,manifest)
    R.dump(MARKER,{"method_version":manifest["method_version"],"status":"MANIFEST_PINNED_BEFORE_DIFF_PIXELS",
        "manifest_sha256":R.sha(old.read_bytes()),"s3_probe_sha256":manifest["probe_file_sha256"],
        "visits_selected":len(manifest["visits"]),
        "independently_eligible_pairs":sum(sum(bool(x) for x in v.get("pair_independence",{}).values())
             for v in manifest["visits"])})
    print("S3_MANIFEST_FROZEN",len(manifest["visits"]),flush=True)
if __name__=="__main__":main()
