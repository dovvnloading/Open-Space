"""Create the coverage, cut-flow, peak and injection-recovery tables on the runner."""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import json
from pathlib import Path

import numpy as np


def write_csv(path, rows, fields):
    with path.open("w", newline="") as stream:
        w = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def main(folder):
    inv = json.loads((folder / "inventory.json").read_text())
    roster = {r["id"]:r for r in inv["sources"]}
    srcs, peaks, injections, cases, failures = [], [], [], [], []
    by = defaultdict(list)
    counts = Counter()
    for path in sorted(folder.glob("search-*.jsonl.gz")):
        with gzip.open(path, "rt") as stream:
            for line in stream:
                s = json.loads(line)
                sid = s["roster_id"]
                original = roster[sid]
                d = dict(source_id=sid, source=original["name"], source_type=original.get("source_type"),
                         maxi_id=original.get("maxi_id"), match=original["match"], status=s["status"])
                counts[s["status"]] += 1
                if s["status"] == "analysis_error":
                    failures.append([sid,s.get("traceback")])
                for band in ("bat","maxi"):
                    b = s.get("instruments", {}).get(band, {})
                    for k in ("status","samples","span_days","min_period_days","max_period_days","model_ok",
                              "local_maxima","above_mc_threshold","weighted_rms"):
                        d[band+"_"+k] = b.get(k)
                    for p in b.get("peaks",[]):
                        peaks.append(dict(source_id=sid, instrument=band, **p))
                        counts[band+"_significant_before_windows"] += 1
                        if not p["window_reason"] and p["model_ok"] and p["survey_fap_bound"] <= .05:
                            counts[band+"_eligible_peaks"] += 1
                d["matching_pairs"] = len(s.get("pairs", []))
                d["stable_pairs"] = len(s.get("candidates", []))
                counts["matching_pairs"] += d["matching_pairs"]
                counts["stable_pairs"] += d["stable_pairs"]
                srcs.append(d)
                for c in s.get("candidates",[]):
                    row = dict(source_id=sid, source=original["name"], period_days=c["period_days"],
                               uncertainty_days=c["uncertainty_days"], survey_fap_bound=c["survey_fap_bound"])
                    for band in ("bat","maxi"):
                        for k in ("period_days","uncertainty_days","semiamplitude","semiamplitude_se","source_global_mc_fap","survey_fap_bound"):
                            row[band+"_"+k]=c[band][k]
                    cases.append(row)
                for inj in s.get("injections",[]):
                    z=dict(source_id=sid,period_days=inj["period_days"],relative_semiamplitude=inj["relative_semiamplitude"],
                           trials=inj["trials"], joint=inj["joint_recovered"])
                    for band in ("bat","maxi"):
                        z[band+"_semiamplitude"]=inj["rate_semiamplitudes"].get(band)
                        z[band+"_recovered"]=inj["recovered_by_band"].get(band)
                    injections.append(z)
                    for group, numerator in dict(joint=z["joint"],bat=z["bat_recovered"],maxi=z["maxi_recovered"]).items():
                        if numerator is not None:
                            by[(group,z["period_days"],z["relative_semiamplitude"])].append((z["trials"],numerator,z))
    srcs.sort(key=lambda row:row["source_id"])
    write_csv(folder/"sources.csv", srcs, list(srcs[0]))
    write_csv(folder/"threshold-peaks.csv", peaks, ["source_id","instrument","period_days","score","source_global_mc_fap","survey_fap_bound","window_reason","model_ok"])
    write_csv(folder/"injections.csv", injections, ["source_id","period_days","relative_semiamplitude","trials","bat_semiamplitude","maxi_semiamplitude","bat_recovered","maxi_recovered","joint"])
    write_csv(folder/"stable-periods.csv", cases, ["source_id","source","period_days","uncertainty_days","survey_fap_bound"] + [
        band+"_"+k for band in ("bat","maxi") for k in ("period_days","uncertainty_days","semiamplitude","semiamplitude_se","source_global_mc_fap","survey_fap_bound")])
    sensitivity=[]
    rng=np.random.default_rng(730321)
    for (group,period,amp), rows in sorted(by.items()):
        trials,rec = np.array([[r[0],r[1]] for r in rows]).T
        # Resample independent source records, keeping each source's injected phases together.
        idx=rng.integers(len(rows),size=(2000,len(rows)))
        lo,hi=np.quantile(rec[idx].sum(1)/trials[idx].sum(1),[.025,.975])
        row=dict(group=group,period_days=period,relative_semiamplitude=amp,sources=len(rows),
                 trials=int(trials.sum()), recovered=int(rec.sum()), fraction=float(rec.sum()/trials.sum()),
                 source_bootstrap_lo=float(lo),source_bootstrap_hi=float(hi))
        for band in ("bat","maxi"):
            amps=[r[2][band+"_semiamplitude"] for r in rows if r[2][band+"_semiamplitude"] is not None]
            row[band+"_median_semiamplitude"]=float(np.median(amps)) if amps else None
        sensitivity.append(row)
    if sensitivity:
        write_csv(folder/"sensitivity.csv", sensitivity, list(sensitivity[0]))
    summary={"counts":counts, "analysis_errors":failures, "joint_injection_recovery":[r for r in sensitivity if r["group"]=="joint"],
             "numerical_survivors":cases}
    (folder/"search-summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    if failures:
        raise ValueError(f"{len(failures)} sources encountered analysis errors; see search-summary.json")


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder",type=Path)
    main(parser.parse_args().folder)
