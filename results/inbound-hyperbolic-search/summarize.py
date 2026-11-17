"""Join the frozen catalogues and independent fits by JPL SPK identifier."""
import argparse
import collections
import csv
import gzip
import io
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/"catalog"))
import audit


def read_catalogue():
    raw = gzip.decompress((ROOT/"catalog/sbdb.json.gz").read_bytes())
    dataset = json.loads(raw)
    return [audit.process(dict(zip(dataset["fields"], row, strict=True))) for row in dataset["data"]]


def original_fields(fit):
    """Require the requested reference frame, inbound motion and a far boundary."""
    if not fit or fit.get("status") != "fit returned":
        return None
    if fit["center"].lower() not in {"object_-1", "barycenter", "solar system barycenter", "solar system barycentre"} or str(fit["inbound"]) != "True":
        return None
    if float(fit["r_au"]) < 100 or abs(float(fit["epoch"])-float(fit["requested_epoch"])) > .01:
        return None
    inv = float(fit["inv_a_au"])
    sigma = fit.get("sigma_e")
    z = (float(fit["e"])-1)/float(sigma) if sigma and float(sigma)>0 else None
    return {"inbound_radius_au":round(float(fit["r_au"]),5), "inv_a_au_inv":inv,
            "z_e_one_model":z,"v_kms":math.sqrt(-inv*audit.GM_SUN/audit.AU_KM*1.001341838645975) if inv<0 else None,
            "n_used":fit["used"],"n_total":fit["total"],"rms_arcsec":fit["rms_arcsec"],
            "used_arc_days":float(fit["last_used_jd"])-float(fit["first_used_jd"])}


def read_softened_shards():
    paths=sorted((ROOT/"batch/evidence/unpenalized").glob("*.csv.gz"))
    assert len(paths)==12, "Expected twelve independent slices"
    result={}
    for path in paths:
        with gzip.open(path,"rt") as stream:
            for row in csv.DictReader(stream):
                key=row["spkid"]
                assert key not in result
                result[key]=row
    return result


def main(out):
    rows = read_catalogue()
    relaxed=read_softened_shards()
    assert {str(r["spkid"]) for r in rows}==set(relaxed)
    default_path=ROOT/"batch/evidence/refits.csv.gz"
    partial_path=ROOT/"batch/evidence/standard_partial.csv.gz"
    if default_path.exists() or partial_path.exists():
        with gzip.open(default_path if default_path.exists() else partial_path,"rt") as f:
            recorded={r["spkid"]:r for r in csv.DictReader(f)}
        if default_path.exists():
            assert set(recorded)==set(relaxed)
    else:
        recorded={}
    batch={k:recorded.get(k,{"status":"standard-pass outcome unavailable"}) for k in relaxed}
    recent_rows = json.loads(gzip.decompress((ROOT/"batch/evidence/latest_mpc_refits.json.gz").read_bytes()))
    initial_rows=json.loads(gzip.decompress((ROOT/"batch/evidence/initial_mpc_queries.json.gz").read_bytes()))
    first_by_id={str(r["spkid"]):r.get("fit",{}) for r in initial_rows}
    for r in recent_rows:
        first=first_by_id.get(str(r["spkid"]),{})
        second=r.get("fit",{})
        if second.get("status")!="fit returned" and first.get("status")=="fit returned" and first.get("astrometry_sha256")==second.get("astrometry_sha256"):
            r["fit"]=first
    recent = {str(r["spkid"]):r for r in recent_rows if r["selection"]=="catalogue refresh"}
    with (ROOT/"catalog/code-original.csv").open() as f:
        published = {r["spkid"]:r for r in csv.DictReader(f)}
    focused = json.loads((ROOT/"orbits/data/refits.json").read_text())
    extra=ROOT/"orbits/slow/evidence/results.json"
    if extra.exists():
        focused.update(json.loads(extra.read_text())["objects"])
    # Only the following two identifications include a checked pre-encounter bound state.
    ejections={"1980 E1":"Jupiter","2024 L5":"Saturn"}
    model_checks={"2025 V1","1997 P2","2008 J4","1877 G2"}
    items,states=[],collections.Counter()
    for row in rows:
        key=str(row["spkid"])
        old=batch[key]
        free=relaxed[key]
        update=recent.get(key,{})
        new=update.get("fit",{})
        old_orig=original_fields(old)
        new_orig=original_fields(new)
        alternate=original_fields(free)
        best = new_orig or alternate or old_orig or {}
        code=published.get(key,{})
        pdes=row["pdes"]
        focused_id=row["known_iso"] or pdes.replace(" ","-")
        focus=focused.get(focused_id,{}).get("models",{})
        g250=(focus.get("gravity",{}).get("boundaries",{}).get("250au") or {})
        ng250=(focus.get("marsden_outward",{}).get("boundaries",{}).get("250au") or {})
        disposition="no accepted original-energy determination"
        if row["known_iso"]:
            disposition="known interstellar control"
        elif pdes in ejections:
            disposition="known planetary ejection; checked pre-encounter"
        elif pdes=="2024 U2":
            disposition="independent bound incoming state at 250 au"
        elif pdes in model_checks:
            disposition="see individual astrometric models and uncertainty"
        elif pdes in {"1995 Y1","2007 W1"}:
            disposition="published alternate models or arcs change original energy"
        elif pdes in {"2011 M1", "2020 W5", "2012 C2"}:
            disposition="published short-arc or model-dependent original candidate; see catalogue notes"
        elif code:
            val, sigma=map(float,[code["published_original_1_over_a_micro_au_inv"],code["sigma"]])
            if val>0:
                disposition="published preferred original nominally bound"
            elif val+5*sigma<0:
                disposition="published preferred original negative beyond 5 formal sigma; review"
            else:
                disposition="published preferred original negative within 5 formal sigma"
        elif best:
            disposition="independent gravity-only inbound estimate; origin unvalidated"
        states[disposition]+=1
        items.append({"spkid":key,"designation":pdes,"name":row["full_name"],"known_iso":row["known_iso"],
                      "osculating_vinf_km_s":row["nominal_osculating_v_km_s"],
                      "catalogue_e":row["e"],"catalogue_sigma_e":row["sigma_e"],
                      "catalogue_ng_fields":row["nonzero_model_fields"],"status":disposition,
                      "batch_fit":old["status"],"latest_fit":new.get("status",update.get("error","")),
                      "unpenalized_fit":free["status"],
                      "default_inverse_a":(old_orig or {}).get("inv_a_au_inv"),
                      "unpenalized_inverse_a":(alternate or {}).get("inv_a_au_inv"),
                      "live_inverse_a":(new_orig or {}).get("inv_a_au_inv"),
                      "default_z":(old_orig or {}).get("z_e_one_model"),
                      "unpenalized_z":(alternate or {}).get("z_e_one_model"),
                      "batch_sign_disagreement":bool(alternate and old_orig and alternate["inv_a_au_inv"]*old_orig["inv_a_au_inv"]<0),
                      "inbound_gravity_source":"live MPC" if new_orig else "bulk unpenalized" if alternate else "bulk default" if old_orig else "",
                      **{k:best.get(k) for k in ("inbound_radius_au","inv_a_au_inv","z_e_one_model","v_kms","n_used","n_total","rms_arcsec","used_arc_days")},
                      "CODE_model":code.get("model"),"CODE_original_micro_au_inv":code.get("published_original_1_over_a_micro_au_inv"),
                      "CODE_sigma_micro_au_inv":code.get("sigma"),"focused_models":bool(focus),
                      "independent_gravity_inv_a_250au":g250.get("inverse_a_au_inv"),
                      "independent_gravity_sigma_250au":g250.get("formal_sigma_inverse_a"),
                      "independent_marsden_outward_inv_a_250au":ng250.get("inverse_a_au_inv"),
                      "independent_marsden_outward_sigma_250au":ng250.get("formal_sigma_inverse_a")})
    assert len(items) == len(batch) == len({r["spkid"] for r in items})
    metrics = {
        "catalogue_rows":len(items),
        "bulk_status":dict(collections.Counter(r["batch_fit"] for r in items)),
        "unpenalized_status":dict(collections.Counter(r["unpenalized_fit"] for r in items)),
        "batch_sign_disagreements":sum(r["batch_sign_disagreement"] for r in items),
        "independent_far_inbound":sum(bool(r["inbound_gravity_source"]) for r in items),
        "independent_far_inbound_negative":sum((r["inv_a_au_inv"] or 0)<0 for r in items),
        "inputs_for_inbound":dict(collections.Counter(r["inbound_gravity_source"] for r in items if r["inbound_gravity_source"])),
        "published_CODE":sum(bool(r["CODE_model"]) for r in items),
        "individual_fit_objects":sum(bool(r["focused_models"]) for r in items),
        "interpretation":dict(states)
    }
    out.mkdir(exist_ok=True,parents=True)
    buf=io.StringIO()
    w=csv.DictWriter(buf,fieldnames=list(items[0]),lineterminator="\n")
    w.writeheader()
    w.writerows(items)
    (out/"dispositions.csv.gz").write_bytes(gzip.compress(buf.getvalue().encode(),mtime=0))
    (out/"counts.json").write_text(json.dumps(metrics,indent=2)+"\n")
    print(json.dumps(metrics))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out",type=Path)
    main(parser.parse_args().out)
