"""Join the signed-off source-by-source reading of primary papers to the measurements."""
import argparse
import csv
import json
import math
from pathlib import Path


def join(folder,destination):
    evidence=folder/"evidence"
    decisions=json.loads((folder/"decisions.json").read_text())
    with (evidence/"stable-periods.csv").open() as stream:
        periods=list(csv.DictReader(stream))
    used=set()
    published=[]
    for p in periods:
        match=[(i,x) for i,x in enumerate(decisions) if x["source_id"]==p["source_id"] and
               math.isclose(float(x["period_days"]),float(p["period_days"]),rel_tol=1e-7)]
        if len(match)!=1 or match[0][0] in used:
            raise ValueError(f"Expected one reviewed decision for {p['source_id']} at {p['period_days']} d")
        i,x=match[0]
        used.add(i)
        if not x["explanation"] or not x["references"]:
            raise ValueError(f"Missing evidence for {p['source_id']}")
        published.append(dict(p, decision=x["decision"],explanation=x["explanation"],
                              references=" | ".join(r["url"] for r in x["references"])))
    if len(used)!=len(decisions):
        raise ValueError("The decisions file has rows absent from the frozen survey")
    destination.mkdir(exist_ok=True,parents=True)
    if published:
        with (destination/"period-decisions.csv").open("w",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=list(published[0]))
            writer.writeheader()
            writer.writerows(published)
    doc=["# Review of independent BAT and MAXI periods","",
         "BAT and MAXI each passed the blind red-noise and window checks. Each row also passed the period and phase stability checks in the two halves of each light curve. The decision below states what the period can be attributed to after examination of the source, nearby objects, aliases and published periods.",""]
    for row in published:
        match=next(x for x in decisions if x["source_id"]==row["source_id"] and math.isclose(float(x["period_days"]),float(row["period_days"]),rel_tol=1e-7))
        doc.extend([f'## {row["source"]} · {float(row["period_days"]):.8g} days',"",
                    f'**{row["decision"]}.** {row["explanation"]}',"",
                    "| Independent measurement | BAT (15–50 keV) | MAXI (2–20 keV) |",
                    "|---|---:|---:|",
                    "| Period (days, 1σ) | "+" | ".join(f'{float(row[b+"_period_days"]):.9g} ± {float(row[b+"_uncertainty_days"]):.2g}' for b in ("bat","maxi"))+" |",
                    "| Semiamplitude (native units, 1σ) | "+" | ".join(f'{float(row[b+"_semiamplitude"]):.5g} ± {float(row[b+"_semiamplitude_se"]):.2g}' for b in ("bat","maxi"))+" |",
                    "| Full-grid Monte Carlo FAP | "+" | ".join(f'{float(row[b+"_source_global_mc_fap"]):.6g}' for b in ("bat","maxi"))+" |",
                    "| Model-conditional survey FAP bound | "+" | ".join(f'{float(row[b+"_survey_fap_bound"]):.4g}' for b in ("bat","maxi"))+" |","",
                    " · ".join(f'[{r["label"]}]({r["url"]})' for r in match["references"]),""])
    doc.extend(["BAT semiamplitudes are counts cm⁻² s⁻¹ and MAXI semiamplitudes are photons cm⁻² s⁻¹. The Monte Carlo rank resolution is 1/256. The smaller analytic bounds are conditional on the fitted continuum. See the [phase folds](folds/index.md) and [numerical alias audit](index.md).",""])
    (destination/"period-decisions.md").write_text("\n".join(doc))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("folder",type=Path)
    p.add_argument("destination",type=Path)
    args=p.parse_args()
    join(args.folder,args.destination)
