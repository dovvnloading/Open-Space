"""Recompute the strict blind-search funnel and the predeclared held-out lithium/control test."""
import csv
import json
from collections import defaultdict
from pathlib import Path
import numpy as np

ROOT = Path(__file__).parent


def summary():
    lookup = {r["targetid"]:r for folder in ("catalog","catalog_holdout")
              for r in csv.DictReader((ROOT/folder/"audit.csv").open())}
    results = defaultdict(list)
    cutflow = {}
    for p in sorted((ROOT/"scans").glob("scan-*.jsonl")):
        key = p.stem.removeprefix("scan-")
        counts = defaultdict(int)
        for line in p.open():
            row = json.loads(line)
            m = row["meta"]
            t = str(m["targetid"])
            a = lookup[t]
            complete = a["all_assigned"] == "True" and a["star_class_consistent"] == "True" and a["ambiguous_gaia"]=="False"
            for prefix in (["local","strict"] if complete else ["local"]):
                counts[prefix+"_stars"] += 1
                counts[prefix+"_observations"] += len(m["exposures"])
                counts[prefix+"_distance"] += bool(a["r_med_geo"])
                for arm in row["arms"].values():
                    counts[prefix+"_triggers"] += sum(len(h) for h in arm["exposure_hits"].values())
                    counts[prefix+"_persistent"] += len(arm["candidates"])
        cutflow[key] = dict(counts)
    for r in csv.DictReader((ROOT/"measurements/lithium-main-bright-49134.csv").open()):
        if r["half_window_px"] == "15" and r["observed_snr"]:
            results[(r["targetid"],r["rest_A"])].append(r)
    stats = defaultdict(list)
    for (_, wave), rows in results.items():
        sig = np.array([float(r["flux_sigma17"]) for r in rows])
        amplitude = np.array([float(r["observed_flux17"])-float(r["template_flux17"]) for r in rows])
        signed_ew = [float(r["observed_ew_A"])-float(r["template_ew_A"]) for r in rows]
        sn = np.sum(amplitude/sig**2)/np.sqrt(np.sum(1/sig**2))
        stats[wave].append((sn, np.median(signed_ew), len(rows)))
    comparison = {}
    for w, values in stats.items():
        vals = np.asarray(values)
        comparison[w] = dict(stars=len(vals), three_visits=int((vals[:,2]==3).sum()),
                             ge_3sigma=int((vals[:,0]>=3).sum()), ge_7sigma=int((vals[:,0]>=7).sum()),
                             le_minus3sigma=int((vals[:,0]<=-3).sum()),
                             median_ew_difference_A=float(np.median(vals[:,1])))
    return dict(shards=cutflow, heldout=comparison)


if __name__=="__main__":
    print(json.dumps(summary(),indent=2))
