"""Audit field denominators using only the committed evidence, with no network."""
import argparse
import csv
import gzip
import json
import math
from collections import Counter
from pathlib import Path


def yes(x):
    return str(x).lower() == "true"


def read(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as stream:
        yield from csv.DictReader(stream)


def rate(k, n):
    if not n:
        return {"hits":k, "trials":n, "fraction":None, "wilson95":None}
    z = 1.959963984540054
    p = k / n
    center = (p + z*z/(2*n))/(1+z*z/n)
    half = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n))/(1+z*z/n)
    return {"hits":k,"trials":n,"fraction":p,"wilson95":[max(0,center-half),min(1,center+half)]}


def audit(out, shared=None):
    shared = shared or (out.parent if out.name == "wide" else out)
    selections = {row["system_id"]:row for row in read(shared/"system_dispositions.csv") if ("W" in row["system_id"][:2]) == (out.name=="wide")}
    groups = {k:Counter() for k in ("science","control")}
    selected = set()
    incomplete = set()
    for field in read(out/"gaia_fields.csv.gz"):
        ident,stage = field["system_id"],field["stage"]
        group = groups[field["cohort"]]
        group["fields"] += 1
        if yes(field.get("dr3_incomplete",False)):
            incomplete.add(ident)
            group["dr3_service_incomplete"] += 1
        if stage in ("fewer_than_two_clean_fpr_components_in_aperture","pair_geometry_or_fpr_flux_contrast"):
            continue
        group["geometry_eligible"] += 1
        if stage != "passed_gaia":
            continue
        selected.add(ident)
        item = selections[ident]
        group["gaia"] += 1
        if item["imaging_status"] in ("no_robust_extended_light","extended_light"):
            group["image_test_conclusive"] += 1
        if yes(item["full_morphology"]):
            group["morphology"] += 1
        if yes(item["plausible_before_novelty"]):
            group["anchored"] += 1
        if not yes(item["prior"]) and yes(item["plausible_before_novelty"]):
            group["unlisted_for_review"] += 1
    if set(selections) != selected:
        raise ValueError("Survivor IDs differ between field and system tables")
    output = {}
    for cohort,g in groups.items():
        output[cohort] = dict(counts=dict(g),
            gaia_per_aperture=rate(g["gaia"],g["fields"]),
            gaia_per_geometric_multiplet=rate(g["gaia"],g["geometry_eligible"]),
            light_plus_model_per_gaia=rate(g["morphology"],g["gaia"]),
            light_plus_model_per_conclusive_image=rate(g["morphology"],g["image_test_conclusive"]),
            selected_per_aperture=rate(g["anchored"],g["fields"]))
    b = list(read(out/"known_recovery.csv"))
    eligible = [r for r in b if yes(r["in_eligible_parent"])]
    unknown = [r["name"] for r in eligible if set(r["fields"].split()) & incomplete]
    output["benchmark"] = dict(total=len(b),eligible=len(eligible), partly_unqueried=unknown,
        gaia_recovery_eligible=rate(sum(yes(r["passed_gaia"]) for r in eligible),len(eligible)),
        light_plus_model_recovery_eligible=rate(sum(yes(r["passed_light_and_model"]) for r in eligible),len(eligible)),
        end_to_end_recovery_eligible=rate(sum(yes(r["passed_with_quasar_anchor"]) for r in eligible),len(eligible)),
        end_to_end_all=rate(sum(yes(r["passed_with_quasar_anchor"]) for r in b),len(b)))
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=Path)
    parser.add_argument("--write",action="store_true")
    args = parser.parse_args()
    result = audit(args.directory)
    packed = json.dumps(result,indent=2,sort_keys=True)+"\n"
    dest = args.directory/"audit.json"
    if args.write:
        dest.write_text(packed)
    elif dest.exists() and json.loads(dest.read_text()) != result:
        raise SystemExit("The deposited audit and the reconstructed audit disagree")
    print(packed)
