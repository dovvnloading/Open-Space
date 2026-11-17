"""Add measured Gaia distance intervals and isotropic-equivalent line powers to flux bounds."""
import argparse
import csv
import math
from pathlib import Path


def power(flux17, distance_pc):
    return 4*math.pi*(float(distance_pc)*3.085677581491367e16)**2*float(flux17)*1e-20


def convert(limits, audits, output):
    lookup = {r["targetid"]:r for p in audits for r in csv.DictReader(p.open())}
    with limits.open() as source, output.open("w",newline="") as dest:
        reader = csv.DictReader(source)
        base = ["targetid","arm","source_id","r_med_geo","r_lo_geo","r_hi_geo"]
        levels = [c for c in reader.fieldnames if c.startswith("flux") and c.endswith("_17")]
        fields = base + [c for c in reader.fieldnames if c not in base and c not in levels]
        for col in levels:
            stem = col.removeprefix("flux").removesuffix("_17")
            fields.extend([col]+[f"eirp{stem}_{q}_W" for q in ("med","lo","hi")])
        writer = csv.DictWriter(dest,fieldnames=fields)
        writer.writeheader()
        for row in reader:
            obj = lookup[row["targetid"]]
            if obj["all_assigned"]!="True" or obj["star_class_consistent"]!="True" or obj["ambiguous_gaia"]!="False":
                raise ValueError(f"non-audited exposure set for {row['targetid']}")
            row.update({k:obj[k] for k in base[2:]})
            for col in levels:
                if not row[col] or not row["r_med_geo"]:
                    continue
                stem = col.removeprefix("flux").removesuffix("_17")
                for q in ("med","lo","hi"):
                    row[f"eirp{stem}_{q}_W"] = format(power(row[col],row[f"r_{q}_geo"]),".8g")
            writer.writerow(row)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("limits",type=Path)
    p.add_argument("output",type=Path)
    p.add_argument("audits",nargs="+",type=Path)
    a=p.parse_args()
    convert(a.limits,a.audits,a.output)
