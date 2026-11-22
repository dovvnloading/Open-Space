"""Inspect Gaia photometry, Pan-STARRS morphology and AllWISE at retained Chandra events."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re

import campaign
from report import FINAL_VETO
import target_optical as optical
from vet import CAMPAIGN_REF


SERVICES = ("gaia_cone", "ps1_mean", "ps1_stack", "allwise")
FIELDS = ("obsid", "input_sha256", "trigger_key", "ra_deg", "dec_deg", "position_error_arcsec", "event_mjd")
SCHEMA = 1


def identity():
    here = Path(__file__).parent
    return {name: hashlib.sha256((here/name).read_bytes()).hexdigest()
            for name in ("lead_optical.py", "target_optical.py", "vet.py")}


def inputs(evidence, selection):
    """Validate the lead table and exact product against its observation ledger."""
    wanted = {}
    if selection:
        for row in json.loads(selection.read_bytes()):
            ident = row["obsid"], row["input_sha256"], row["trigger_key"]
            if ident in wanted:
                raise ValueError("Repeated lead identity")
            wanted[ident] = row
        files = [evidence/"observations"/f"{oid}.json.gz" for oid in sorted({k[0] for k in wanted})]
    else:
        files = sorted((evidence/"observations").glob("*.json.gz"))
    found = set()
    for path in files:
        ledger = campaign.read_obs(path)
        if ledger["status"] not in {"complete", "followup_incomplete"}:
            continue
        oid = ledger["obsid"]
        for product in ledger["products"]:
            if not product.get("complete"):
                continue
            sha = product["input"]["sha256"]
            for event in product["search"]["triggers"]:
                if not event["instrument_pass"] or event["vet"]["status"] in FINAL_VETO:
                    continue
                key = event["trigger_key"]
                ident = oid, sha, key
                if selection and ident not in wanted:
                    continue
                if (not re.fullmatch(r"real/(0|[1-9][0-9]*)", key) or
                        not re.fullmatch(r"[0-9a-f]{64}", sha) or oid != int(path.name.split(".")[0])):
                    raise ValueError("Invalid retained trigger identity")
                lead = {"obsid": oid, "input_sha256": sha, **{k: event[k] for k in FIELDS[2:]}}
                if selection and any(wanted[ident][k] != y for k, y in lead.items()):
                    raise ValueError("Retained lead disagrees with its individual trigger")
                matches = event["vet"].get("matches", [])
                ids = {str(m["name"]) for m in matches if m["catalog"] == "gaia"}
                if any(not re.fullmatch(r"[0-9]{1,20}", x) for x in ids):
                    raise ValueError("Invalid Gaia identifier")
                found.add(ident)
                yield lead, sorted(ids)
    if selection and set(wanted) != found:
        raise ValueError("Some retained events have no corresponding primary ledger")


def run_catalog(name, lead, cache, replay):
    out = {"status":"incomplete", "rows":[], "receipts":{}}
    try:
        if name.startswith("ps1_"):
            optical.ps1_service(name[4:], lead, cache, replay, out)
        elif name.startswith("gaia_id_"):
            prior = name[len("gaia_id_"):]
            optical.tap_service(name, lead, cache, replay, out, known=prior, prior_gaia_id=prior)
        else:
            optical.tap_service(name, lead, cache, replay, out)
        out["rows"].sort(key=lambda r: r["offset_arcsec"] if r["offset_arcsec"] is not None else math.inf)
        out["status"] = "rows_returned" if out["rows"] else "empty_validated_cone"
    except Exception as exc:
        out["error"] = f"{type(exc).__name__}: {exc}"
    return out


def save_index(out, index, failed):
    path = out/"index.csv"
    temp = path.with_suffix(".tmp")
    fields = (*FIELDS, "catalog_status","gaia_ids_matched_by_primary","record",
              "gaia_sources_inside_r90","gaia_variables_inside_r90","ps1_stack_rows_inside_r90",
              "allwise_rows_inside_r90")
    with temp.open("w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fields)
        writer.writeheader()
        writer.writerows(index)
    temp.replace(path)
    campaign.dump(out/"summary.json", {"lead_count":len(index), "incomplete":failed, "schema":SCHEMA,
                  "index_sha256":hashlib.sha256(path.read_bytes()).hexdigest(),
                  "note":"Human-review context for real leads. The primary false-alarm and injection selections "
                         "use the registered pipeline. Gaia photometric flags, optical morphology and static "
                         "catalog associations alone cannot identify an X-ray transient."})


def run(evidence, out, selection, replay):
    if (os.getenv("GITHUB_ACTIONS") != "true" or os.getenv("GITHUB_REF") != CAMPAIGN_REF or
            os.getenv("GITHUB_EVENT_NAME") not in {"push", "workflow_dispatch"}):
        raise RuntimeError("Use the Chandra research workflow")
    cache = out/"source"
    cache.mkdir(parents=True, exist_ok=True)
    source_code = identity()
    index = []
    failed = []
    for lead, prior in inputs(evidence, selection):
        oid, sha, key = (lead[n] for n in FIELDS[:3])
        path = out/"observations"/str(oid)/sha/(key+".json.gz")
        obj = campaign.read_obs(path) if path.is_file() else {}
        if not (obj.get("lead") == lead and obj.get("prior_gaia_ids") == prior
                and obj.get("implementation_sha256") == source_code and obj.get("complete") and not replay):
            obj = {"lead":lead, "prior_gaia_ids":prior, "epoch_coordinate_system":"ICRS",
                   "schema":SCHEMA, "implementation_sha256":source_code, "search_radius_arcsec":optical.CONE,
                   "catalogs":{}, "complete":True}
            for name in (*SERVICES, *("gaia_id_"+n for n in prior)):
                catalog = run_catalog(name, lead, cache, replay)
                obj["catalogs"][name] = catalog
                if catalog["status"] == "incomplete":
                    obj["complete"] = False
            path.parent.mkdir(parents=True, exist_ok=True)
            campaign.save_obs(path, obj)
        catalogs = obj["catalogs"]
        if not obj["complete"]:
            failed.append({"obsid":oid,"input_sha256":sha,"trigger_key":key,
                           "catalogs":[k for k,v in catalogs.items() if v["status"]=="incomplete"]})
        def inside(name):
            return [v for v in catalogs[name]["rows"] if v.get("within_xray_r90") is True]
        near = {str(v["catalog_values"]["source_id"]):v["catalog_values"] for name in catalogs
                if name.startswith("gaia_") for v in inside(name)}
        index.append({**lead, "catalog_status":"complete" if obj["complete"] else "incomplete",
                      "gaia_ids_matched_by_primary":" ".join(prior), "record":path.relative_to(out).as_posix(),
                      "gaia_sources_inside_r90":" ".join(sorted(near)),
                      "gaia_variables_inside_r90":" ".join(sorted(k for k,v in near.items()
                                                                  if v["phot_variable_flag"]=="VARIABLE")),
                      "ps1_stack_rows_inside_r90":len(inside("ps1_stack")),
                      "allwise_rows_inside_r90":len(inside("allwise"))})
        if len(index) % 25 == 0:
            print(f"Reviewed {len(index)} fixed optical positions", flush=True)
            save_index(out, index, failed)
    save_index(out, index, failed)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--from-ledgers", action="store_true", help="Use deposited pre-deduplication primary rows")
    parser.add_argument("--replay", action="store_true", help="Verify or reconstruct using the saved response bytes")
    args = parser.parse_args()
    run(args.evidence, args.out, None if args.from_ledgers else args.evidence/"leads.json", args.replay)
