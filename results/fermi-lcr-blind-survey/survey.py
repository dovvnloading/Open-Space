"""Run the frozen tests once per catalog source, then calibrate across the roster."""
import argparse
import csv
import gzip
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np

import acquire
import periods
import states


def clean(value):
    if isinstance(value, dict):
        return {k: clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [clean(v) for v in value]
    if isinstance(value, (np.bool_, np.integer, np.floating)):
        return clean(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(clean(value), indent=2, allow_nan=False) + "\n")


def table(path, rows):
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v) if isinstance(v, (list, dict)) else v
                             for k, v in clean(row).items()})


def shard(args):
    data = Path(args.data)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    catalog = json.loads(Path(args.index).read_text())
    rows, nulls = [], {}
    for src in catalog[args.shard::args.shards]:
        name = src["name"]
        curves = {}
        for cadence in ("monthly", "weekly", "3day"):
            path = data / "curves" / acquire.slug(name) / (cadence + ".npz")
            if not path.is_file():
                continue
            with np.load(path, allow_pickle=False) as series:
                curve = dict(series)
                curve["time"], curve["upper"] = curve["t"], curve["ul"]
                curves[cadence] = curve
        row, null, injection = periods.analyze(name, curves)
        change = states.scan_state(curves.get("monthly"), curves.get("weekly"), curves.get("3day"))
        rows.append({"source": name, "row": row, "injection": injection, "state": change})
        nulls[name] = np.asarray(null, dtype=float) if null is not None else np.array([])
    with gzip.open(out / "calibration.json.gz", "wt") as stream:
        json.dump(clean(rows), stream, allow_nan=False, separators=(",", ":"))
    np.savez_compressed(out / "null.npz", **nulls)


def reduce(args):
    work, out = Path(args.work), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    catalog = json.loads(Path(args.index).read_text())
    records, state_rows, seen, receipts = [], [], set(), []
    source_by_name = {s["name"]: s for s in catalog}
    for part in sorted(work.glob("part-*")):
        with gzip.open(part / "calibration.json.gz", "rt") as stream:
            payload = json.load(stream)
        with np.load(part / "null.npz", allow_pickle=False) as nulls:
            for entry in payload:
                name = entry["source"]
                if name in seen:
                    raise ValueError(f"Repeated source {name}")
                seen.add(name)
                records.append((entry["row"], nulls[name], entry["injection"]))
                state_rows.append({"source": name, **entry["state"]})
        with (part / "manifest.jsonl").open() as src:
            receipts.extend(json.loads(line) for line in src if line.strip())
    if seen != set(source_by_name):
        raise ValueError("Shards do not partition the archive roster")
    if len(receipts) != 3 * len(seen):
        raise ValueError("Manifest must account for all three cadences per source")
    periodic, injections, summary = periods.finalize(records)
    state_rows.sort(key=lambda r: r["source"])
    periodic.sort(key=lambda r: r["source"])
    dump(out / "period-summary.json", summary)
    table(out / "periods.csv", periodic)
    table(out / "injections.csv", injections)
    table(out / "states.csv", [{k:v for k,v in r.items() if k != "transitions"} for r in state_rows])
    transitions = [{"source": r["source"], **t} for r in state_rows for t in r.get("transitions", [])]
    table(out / "transitions.csv", transitions or [])
    with (out / "inputs.jsonl").open("w") as stream:
        for r in receipts:
            stream.write(json.dumps(clean(r), separators=(",", ":")) + "\n")
    leads = []
    for row in periodic:
        if row.get("passes"):
            src = source_by_name[row["source"]]
            leads.append({**src, "kind": "period", "measurement": row})
    for row in state_rows:
        if row.get("transitions"):
            src = source_by_name[row["source"]]
            leads.append({**src, "kind": "state", "measurement": row})
    dump(out / "leads.json", leads)
    dump(out / "coverage.json", {"roster": len(catalog), "cadence_receipts": len(receipts),
         "inputs": dict(Counter(str(r.get("status")) for r in receipts)),
         "periods": dict(Counter(r["status"] for r in periodic)),
         "states": dict(Counter(r["status"] for r in state_rows)),
         "period_pass": sum(bool(r.get("passes")) for r in periodic),
         "state_sources": sum(bool(r.get("transitions")) for r in state_rows),
         "transitions": len(transitions)})


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    sub = cli.add_subparsers(dest="command", required=True)
    s = sub.add_parser("shard")
    s.add_argument("--data", required=True)
    s.add_argument("--shard", type=int, required=True)
    s.add_argument("--shards", type=int, required=True)
    r = sub.add_parser("reduce")
    r.add_argument("--work", required=True)
    for p in (s, r):
        p.add_argument("--out", required=True)
        p.add_argument("--index", required=True)
    opt = cli.parse_args()
    (shard if opt.command == "shard" else reduce)(opt)
