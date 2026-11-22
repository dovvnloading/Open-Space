"""Recovery of periodic 10%-duty outburst trains on frozen dual-instrument sources."""
import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import time

import numpy as np

from field_checks import load_archives, read_rows
from search import ALPHA, INJ_PHASES, Search, match, prepare, rng_for, survey_trials

DUTY = .10
PEAK_RMS = (1., 3., 10.)


def daily_pulse(t, frequency, phase):
    """Mean of a unit-height box train over a one-day bin centered at each t."""
    center = t*frequency + phase
    def primitive(cycles):
        lower = np.floor(cycles)
        return lower*DUTY + np.minimum(cycles-lower,DUTY)
    return np.clip((primitive(center+.5*frequency) - primitive(center-.5*frequency))/frequency, 0., 1.)


def replay_source(source, band, curve, saved, roster_size):
    data, meta = prepare(curve)
    if data is None:
        raise ValueError(f"{band} frozen curve is unavailable or ineligible")
    s = Search(source,band,data,meta,roster_size)
    ref = np.array(saved["red_variances"])
    # Replaying the original inputs, pinned versions and deterministic seed streams
    # provides the complete original thresholds, including the per-frequency training
    # means and full-scan Monte Carlo maxima, which were not exported as arrays.
    if not (survey_trials(len(s.freq),roster_size) == saved["survey_trial_factor"] and
            s.model_ok == saved["model_ok"] and len(s.freq) == saved["frequencies"] and
            s.noise.method == saved["red_fit_method"] and
            np.allclose(s.noise.variances,ref,rtol=1e-6,atol=max(1e-30,abs(ref).max()*1e-8)) and
            np.isclose(s.scale,saved["weighted_rms"],rtol=1e-8,atol=1e-30)):
        raise ValueError(f"{band} original continuum, sampling or normalization failed replay")
    rows, counts = s.peaks(s.base,s.y)
    if [p["index"] for p in rows] != [p["index"] for p in saved["peaks"]] or any(
            n != saved[k] for k,n in counts.items()):
        raise ValueError(f"{band} original blind peak list failed replay")
    for actual, expected in zip(rows,saved["peaks"]):
        if actual["source_global_mc_fap"] != expected["source_global_mc_fap"] or not np.isclose(
                actual["score"], expected["score"], rtol=1e-6):
            raise ValueError(f"{band} original Monte Carlo calibration failed replay")
    return s


def fundamental(engine, sums, y, frequency):
    """Apply the primary blind gates, then identify recoveries near the planted fundamental."""
    scores = engine.scores(sums)
    nearby = np.abs(engine.freq-frequency) <= 1/engine.span
    # These screens only skip calculations for sources with zero possible recoveries.
    if not engine.model_ok or not nearby.any() or scores[nearby].max() <= engine.threshold:
        return []
    recovered = []
    for k, flag in engine.ranked(scores)[0]:
        if flag is not None or not nearby[k]:
            continue
        _, bound = engine.analytic(k,y)
        if bound <= ALPHA and engine.halves(engine.freq[k],y)["pass"]:
            recovered.append({"period_days":float(1/engine.freq[k])})
    return recovered


def pulse_experiment(source, record, curves, roster_size):
    assigned = {r["period_days"] for r in record["injections"]}
    if not assigned:
        return {"source_id":source,"status":"no eligible original period stratum","experiments":[]}
    if len(assigned) != 1:
        raise ValueError("Original injections did not use one common period")
    period = assigned.pop()
    frequency = 1/period
    bands = {name: replay_source(source,name,curves.get(f"{source}.{name}"),record["instruments"][name],roster_size)
             for name in ("bat","maxi")}
    origin = min(s.origin for s in bands.values())
    phases = rng_for(source,"10-percent-duty-pulse-injections").uniform(0,1,INJ_PHASES)
    runs = {a:{"peak_to_baseline_rms":a, "phase_trials":[], "recovered":{"bat":0,"maxi":0,"joint":0},
               "peak_rate_excess":{name:a*s.scale for name,s in bands.items()}} for a in PEAK_RMS}
    for phase in phases:
        pulses = {}
        for name,s in bands.items():
            pulse = daily_pulse((s.origin-origin)+s.t,frequency,phase)
            pulses[name] = (pulse,s.transforms(pulse)[0])
        for a,row in runs.items():
            found = {}
            for name,s in bands.items():
                level = row["peak_rate_excess"][name]
                pulse, psums = pulses[name]
                found[name] = fundamental(s,s.base+level*psums,s.y+level*pulse,frequency)
            joint = any(match(bands["bat"],bands["maxi"],b,m) for b in found["bat"] for m in found["maxi"])
            success = {name:bool(f) for name,f in found.items()}
            success["joint"] = joint
            for name,passed in success.items():
                row["recovered"][name] += int(passed)
            row["phase_trials"].append({"phase_cycles":float(phase),**success})
    return {"source_id":source,"status":"tested","period_days":period,"duty_cycle":DUTY,
            "phase_reference_mjd":origin,"phase_draws":len(phases),"experiments":list(runs.values())}


def shard(args):
    begin, cpu = time.monotonic(), time.process_time()
    args.out.mkdir(parents=True,exist_ok=True)
    inv = json.loads((args.evidence/"inventory.json").read_text())
    ledger = args.evidence/f"search-{args.part:02}.jsonl.gz"
    snapshot = args.snapshot or args.evidence
    records = {r["roster_id"]:r for r in read_rows(ledger) if r["status"] == "dual_search"}
    wanted = {f"{s}.{b}" for s in records for b in ("bat","maxi")}
    curves = load_archives([snapshot/f"curves-{args.part:02}.tar.gz"],wanted)
    statuses = Counter()
    with (args.out/f"pulse-{args.part:02}.jsonl").open("w") as out:
        for i,(sid,record) in enumerate(records.items(),1):
            try:
                row = pulse_experiment(sid,record,curves,len(inv["sources"]))
            except Exception as exc:
                row = {"source_id":sid,"status":"replay error","error":f"{type(exc).__name__}: {exc}","experiments":[]}
            finally:
                Search.precision.cache_clear()
            statuses[row["status"]] += 1
            out.write(json.dumps(row,allow_nan=False,separators=(",",":"))+"\n")
            out.flush()
            if i%5 == 0:
                print(f"Shard {args.part}: {i}/{len(records)} dual sources",flush=True)
    timing = {"shard":args.part,"dual_sources":len(records),"statuses":statuses,
              "wall_seconds":time.monotonic()-begin,"process_cpu_seconds":time.process_time()-cpu}
    (args.out/f"pulse-timing-{args.part:02}.json").write_text(json.dumps(timing,indent=2)+"\n")
    print(json.dumps(timing))
    if statuses["replay error"]:
        raise RuntimeError("One or more frozen pairs could not be replayed")


def summarize(args):
    args.out.mkdir(parents=True,exist_ok=True)
    buckets, statuses, seen = defaultdict(list),Counter(),set()
    for path in sorted(args.input.glob("pulse-[0-9][0-9].jsonl")):
        for r in read_rows(path):
            sid = r["source_id"]
            if sid in seen:
                raise ValueError(f"Duplicate source in pulse sensitivity ledger: {sid}")
            seen.add(sid)
            statuses[r["status"]] += 1
            for ex in r["experiments"]:
                key = (r["period_days"],ex["peak_to_baseline_rms"])
                buckets[key].append((r["phase_draws"],ex))
    if not seen or statuses["replay error"]:
        raise ValueError("Pulse-shard records are missing or contain replay failures")
    if args.evidence:
        expected = {r["roster_id"] for p in args.evidence.glob("search-*.jsonl.gz")
                    for r in read_rows(p) if r["status"]=="dual_search"}
        if seen != expected:
            raise ValueError(f"Incomplete pulse census: {len(seen)} source records, {len(expected)} expected")
    rows = []
    rng = rng_for("survey","pulse-source-bootstrap")
    for (period,a),items in sorted(buckets.items()):
        trials = np.array([n for n,_ in items])
        indices = rng.integers(len(items),size=(2000,len(items)))
        for band in ("bat","maxi","joint"):
            hits = np.array([ex["recovered"][band] for _,ex in items])
            lo,hi = np.quantile(hits[indices].sum(axis=1)/trials[indices].sum(axis=1),[.025,.975])
            rows.append({"instrument":band,"period_days":period,"duty_cycle":DUTY,"peak_to_baseline_rms":a,
                "sources":len(items),"phase_trials":int(trials.sum()),"fundamentals_recovered":int(hits.sum()),
                "recovery_fraction":float(hits.sum()/trials.sum()),
                "source_bootstrap_025":float(lo),"source_bootstrap_975":float(hi),
                **{f"{b}_median_peak_rate_excess":float(np.median([ex["peak_rate_excess"][b] for _,ex in items]))
                   for b in ("bat","maxi")}})
    if rows:
        with (args.out/"pulse-sensitivity.csv").open("w",newline="") as out:
            writer = csv.DictWriter(out,fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    timing = [json.loads(p.read_text()) for p in sorted(args.input.glob("pulse-timing-*.json"))]
    summary = {"source_statuses":statuses,"shard_timing":timing,
               "total_process_cpu_seconds":sum(r["process_cpu_seconds"] for r in timing),
               "joint_recovery":[r for r in rows if r["instrument"]=="joint"]}
    (args.out/"pulse-summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command",required=True)
    run = commands.add_parser("shard")
    run.add_argument("--evidence",type=Path,required=True)
    run.add_argument("--snapshot",type=Path)
    run.add_argument("--part",type=int,required=True)
    run.add_argument("--out",type=Path,required=True)
    summary = commands.add_parser("summarize")
    summary.add_argument("--input",type=Path,required=True)
    summary.add_argument("--evidence",type=Path)
    summary.add_argument("--out",type=Path,required=True)
    args = parser.parse_args()
    (shard if args.command=="shard" else summarize)(args)


if __name__=="__main__":
    main()
