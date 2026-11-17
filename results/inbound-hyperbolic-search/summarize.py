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
    if int(fit["used"])<3:
        return None
    inv = float(fit["inv_a_au"])
    if not math.isfinite(inv):
        return None
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
    sampled_path=ROOT/"batch/evidence/recent_stratified.json.gz"
    sampled={str(r["spkid"]):r for r in json.loads(gzip.decompress(sampled_path.read_bytes()))} if sampled_path.exists() else {}
    with (ROOT/"catalog/code-original.csv").open() as f:
        published = {r["spkid"]:r for r in csv.DictReader(f)}
    focused = json.loads((ROOT/"orbits/data/refits.json").read_text())
    extra=ROOT/"orbits/slow/evidence/results.json"
    if extra.exists():
        focused.update(json.loads(extra.read_text())["objects"])
    fast={}
    for path in sorted((ROOT/"batch/fast-priority/evidence").glob("*shard-*.json")):
        for rec in json.loads(path.read_text())["records"]:
            fast[str(rec["spkid"])]=rec
    review_path=ROOT/"batch/fast-priority/evidence/station-review.json"
    if review_path.exists():
        review=json.loads(review_path.read_text())
        rec=fast[review["spkid"]]
        assert rec["raw_sha256"] == review["input_sha256"]
        rec.update(status=review["reviewed_status"],far_inbound_100au=review["reviewed_far_inbound_100au"],
                   far_inbound_250au=review["reviewed_far_inbound_250au"])
    # Native alias grouping succeeds for controls with permanent and provisional MPC labels.
    for path in sorted((ROOT/"batch/fast-priority/evidence").glob("control-shard-*.json")):
        for rec in json.loads(path.read_text())["records"]:
            control_review_path=ROOT/"batch/fast-priority/evidence/control-station-review.json"
            if rec["status"]!="identity and observation epochs verified" and control_review_path.exists():
                control_review=json.loads(control_review_path.read_text())
                if str(rec["spkid"])==str(control_review.get("spkid")):
                    assert rec["raw_sha256"]==control_review["input_sha256"]
                    rec.update(status=control_review["reviewed_status"],far_inbound_100au=control_review["reviewed_far_inbound_100au"],
                               far_inbound_250au=control_review["reviewed_far_inbound_250au"])
            if rec["far_inbound_100au"]:
                fast[str(rec["spkid"])]=rec
    event={}
    for group in ("census","census_short"):
        path=ROOT/f"orbits/{group}/evidence/summary.csv"
        with path.open() as f:
            for rec in csv.DictReader(f):
                key=str(rec["spkid"])
                assert key not in event
                event[key]=rec
    with (ROOT/"orbits/robust/evidence/dispositions.csv").open() as f:
        robustness={r["spkid"]:r for r in csv.DictReader(f)}
    current_event_path=ROOT/"orbits/recent/evidence/results.json"
    current_event=json.loads(current_event_path.read_text())["objects"] if current_event_path.exists() else {}
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
        sample=sampled.get(key,{}).get("fit",{})
        sampled_orig=original_fields(sample)
        best = new_orig or sampled_orig or alternate or old_orig or {}
        code=published.get(key,{})
        pdes=row["pdes"]
        focused_id=row["known_iso"] or pdes.replace(" ","-")
        focus=focused.get(focused_id,{}).get("models",{})
        g250=(focus.get("gravity",{}).get("boundaries",{}).get("250au") or {})
        ng250=(focus.get("marsden_outward",{}).get("boundaries",{}).get("250au") or {})
        independent_current=current_event.get(key,{}).get("models",{}).get("gravity",{}).get("boundaries",{}).get("250au") or {}
        priority=fast.get(key,{})
        ground=event.get(key,{})
        robustness_row=robustness.get(key,{})
        boundary=independent_current or g250
        exact_inv=boundary.get("inverse_a_au_inv")
        exact_sig=boundary.get("formal_sigma_inverse_a")
        event_source="recent MPC sample" if independent_current else "targeted object fit" if g250 else ""
        if ground.get("inverse_a_au_inv") not in ("",None):
            assert ground["status"] == "inbound 250 au" and float(ground["inbound_radial_au_d"])<0 and abs(float(ground["radius_au"])-250)<0.1
        if exact_inv is None and ground.get("inverse_a_au_inv") not in ("",None):
            exact_inv=float(ground["inverse_a_au_inv"])
            exact_sig=float(ground["formal_sigma_inverse_a"]) if ground.get("formal_sigma_inverse_a") else None
            event_source="seeded ground-based census"
        disposition="no accepted original-energy determination"
        if row["known_iso"]:
            disposition="known interstellar control"
        elif pdes in ejections:
            disposition="known planetary ejection; checked pre-encounter"
        elif pdes in {"2005 D3","2006 S8","2004 D4","2006 S9"}:
            disposition="SOHO sungrazer on a short spacecraft arc; see parabolic and free fits"
        elif pdes=="1996 J1-A":
            disposition="post-breakup fragment; intact parent is compatible with parabolic entry"
        elif pdes in {"2025 K1-D","2025 K1-B","2019 Y4-A","2019 Y4-C"}:
            disposition="documented comet fragment; component orbit is not an independent incoming body"
        elif pdes=="2024 U2":
            disposition="independent bound incoming state at 250 au"
        elif robustness_row:
            disposition=robustness_row["status"]
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
        elif exact_inv is not None:
            disposition="independent 250 au gravity solution; physical origin requires model checks"
        elif priority.get("far_inbound_100au"):
            disposition="longer-budget independent far-inbound fit; see retained arc and quality"
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
                      "station_night_sample_fit":sample.get("status"),
                      "station_night_sample_inverse_a":(sampled_orig or {}).get("inv_a_au_inv"),
                      "default_z":(old_orig or {}).get("z_e_one_model"),
                      "unpenalized_z":(alternate or {}).get("z_e_one_model"),
                      "batch_sign_disagreement":bool(alternate and old_orig and alternate["inv_a_au_inv"]*old_orig["inv_a_au_inv"]<0),
                      "inbound_gravity_source":"live MPC" if new_orig else "live MPC station/night sample" if sampled_orig else "bulk unpenalized" if alternate else "bulk default" if old_orig else "",
                      **{k:best.get(k) for k in ("inbound_radius_au","inv_a_au_inv","z_e_one_model","v_kms","n_used","n_total","rms_arcsec","used_arc_days")},
                      "long_budget_status":priority.get("status"),"long_budget_far_100au":priority.get("far_inbound_100au"),
                      "long_budget_inv_a":priority.get("inv_a") if priority.get("far_inbound_100au") else None,
                      "long_budget_formal_e_z":priority.get("formal_e_excess_sigma") if priority.get("far_inbound_100au") else None,
                      "long_budget_used":priority.get("used"),"long_budget_input_days":priority.get("input_days"),"long_budget_used_days":priority.get("used_days"),
                      "activity_or_holdout_status":robustness_row.get("status"),"activity_or_holdout_note":robustness_row.get("note"),
                      "ground_census_status":ground.get("status"),"ground_census_reason":ground.get("reason"),
                      "event_250au_source":event_source,"event_250au_inv_a":exact_inv,"event_250au_sigma_inv_a":exact_sig,
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
        "records_with_actual_inbound_250au":sum(r["event_250au_inv_a"] is not None for r in items),
        "records_with_any_far_inbound_independent_fit":sum(bool(r["inbound_gravity_source"] or r["long_budget_far_100au"] or r["event_250au_inv_a"] is not None) for r in items),
        "records_with_250au_signed_uncertainty":sum(r["event_250au_sigma_inv_a"] is not None for r in items),
        "independent_far_inbound_long_budget":sum(bool(r["long_budget_far_100au"]) for r in items),
        "independent_far_or_published_CODE":sum(bool(r["inbound_gravity_source"] or r["long_budget_far_100au"] or r["event_250au_inv_a"] is not None or r["CODE_model"]) for r in items),
        "records_with_negative_250au_gravity":sum((r["event_250au_inv_a"] or 0)<0 for r in items),
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
