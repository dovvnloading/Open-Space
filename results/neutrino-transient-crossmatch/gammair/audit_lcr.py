"""Preserve every returned bin's selection fields and audit all preregistered quality decisions offline."""
import argparse
import csv
import gzip
import json
import itertools
import lzma
from pathlib import Path
from fetch_lcr import HERE
from freeze_lcr import extract, keyed

def compact(name, raw):
    # MET is the common key in all four LCR arrays. Repeated week-length
    # differences compress efficiently and preserve arbitrary missing weeks.
    ts=keyed(raw["ts"],2)
    times=list(ts)
    positions={t:i for i,t in enumerate(times)}
    flux=keyed(raw["flux"],2)
    errors=keyed(raw["flux_error"],3)
    codes=keyed(raw["fit_convergence"],2)
    if set(flux)!=set(errors) or set(codes)!=set(ts) or set(flux)-set(ts):
        raise ValueError(f"inconsistent NASA time keys for {name}")
    return dict(name=name,met_delta=times[:1]+[b-a for a,b in zip(times,times[1:])],
                ts=[ts[t][0] for t in times], convergence=[codes[t][0] for t in times],
                detected_at=[positions[t] for t in flux],flux=[flux[t][0] for t in flux],
                lower=[errors[t][0] for t in flux],upper=[errors[t][1] for t in flux])

def expand(record):
    times=list(itertools.accumulate(record["met_delta"]))
    detected=[times[i] for i in record["detected_at"]]
    return dict(ts=list(zip(times,record["ts"],strict=True)),
                fit_convergence=list(zip(times,record["convergence"],strict=True)),
                flux=list(zip(detected,record["flux"],strict=True)),
                flux_error=list(zip(detected,record["lower"],record["upper"],strict=True)))

def pack(cache):
    folder=HERE/"all_fit_inputs"
    folder.mkdir(exist_ok=True)
    paths=sorted(cache.glob("4FGL_J*.json.gz"))
    for i in range(0,len(paths),50):
        records=[]
        for f in paths[i:i+50]:
            raw=json.loads(gzip.decompress(f.read_bytes()))
            name=f.name.split(".json",1)[0].replace("_"," ")
            records.append(compact(name,raw))
        payload=json.dumps(records,separators=(",",":")).encode()
        (folder/f"fits_{i//50:02}.json.xz").write_bytes(lzma.compress(payload))

def audit():
    with (HERE/"lcr_source_counts.csv").open() as f:
        counts={row["source_id"]:int(row["bins_before_cutoff"]) for row in csv.DictReader(f)}
    bins={name:[] for name in counts}
    for file in sorted((HERE/"weekly_inputs").glob("*.gz")):
        with gzip.open(file,"rt") as f:
            for r in csv.DictReader(f):
                bins[r["source_id"]].append((int(r["met_mid"]),float(r["ts"]),float(r["flux"]),float(r["lo"]),float(r["hi"])))
    seen=set()
    for file in sorted((HERE/"all_fit_inputs").glob("*.xz")):
        for row in json.loads(lzma.decompress(file.read_bytes())):
            ident=row["name"]
            if ident in seen:
                raise ValueError(f"duplicate source {ident}")
            seen.add(ident)
            all_bins, surviving_bins=extract(expand(row))
            if counts[ident]!=len(all_bins) or bins[ident]!=surviving_bins:
                raise ValueError(f"quality audit failed for {ident}")
    if seen!=set(counts):
        raise ValueError("missing source")
    print(f"Verified {len(seen)} sources, {sum(counts.values())} weekly fits, and {sum(map(len,bins.values()))} qualifying fits.")

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache",type=Path,help="pack new full fit-selection inputs before auditing")
    args=p.parse_args()
    if args.cache:
        pack(args.cache)
    audit()
