"""Bundle the requested wide-sample follow-up in checked, size-bounded per-cap files."""
import argparse
from collections import defaultdict
import csv
import gzip
import hashlib
import json
from pathlib import Path

BASE=Path(__file__).resolve().parent
CONTROL="J085715.68+560405.8"
MAX_BYTES=850000

def line(value):
    return (json.dumps(value,separators=(",",":"))+"\n").encode()

def dump_parts(lines, stem, dest, index):
    def write_group(group):
        payload=gzip.compress(b"".join(group),mtime=0)
        if len(payload)>MAX_BYTES:
            if len(group)<2:raise ValueError(f"Single oversized record in {stem}")
            mid=len(group)//2
            write_group(group[:mid]);write_group(group[mid:])
            return
        key=f"{stem}-{len(parts)+1:02d}.jsonl.gz"
        (dest/key).write_bytes(payload)
        parts.append(dict(file=key,compressed_bytes=len(payload),sha256=hashlib.sha256(payload).hexdigest()))
    parts=[];write_group(lines);index[stem]=parts

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--selected",type=Path,required=True);ap.add_argument("--matches",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True);ap.add_argument("--partial",action="store_true")
    args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    known={r["source_name"] for r in csv.DictReader(args.matches.open())}
    targets=[r for r in csv.DictReader(args.selected.open()) if r["source_name"] not in known or r["source_name"]==CONTROL]
    missing=[];by_cap=defaultdict(lambda:{"records":[],"neighbours":[]});summary=[]
    for target in targets:
        name=target["source_name"];record=BASE/"retained"/(name+".json");path=BASE/"retained"/(name+".csv.gz")
        if not record.exists() or not path.exists():
            missing.append(name);continue
        z=json.loads(record.read_text());data=gzip.decompress(path.read_bytes())
        if z["retained_csv_sha256"]!=hashlib.sha256(data).hexdigest():raise ValueError(f"Changed subset {name}")
        if (z["ra"],z["dec"],z["pmra"],z["pmdec"])!=tuple(float(target[k]) for k in ["ra2015_deg","dec2015_deg","pmra","pmdec"]):
            raise ValueError(f"Stale trajectory {name}")
        cap=by_cap[target["field"]]
        cap["records"].append(line(z))
        # A rare very large individual CSV may contribute multiple string pieces.
        # Concatenating by name and piece index exactly restores its hashed CSV bytes.
        sz=600000
        for i in range(0,len(data),sz):
            cap["neighbours"].append(line([name,i//sz,data[i:i+sz].decode("ascii")]))
        summary.append(z["summary"])
    if missing and not args.partial:raise ValueError(f"{len(missing)} missing of {len(targets)}")
    with (args.out/"summary.csv").open("w") as f:
        w=csv.DictWriter(f,fieldnames=summary[0]);w.writeheader();w.writerows(summary)
    index={}
    for cap,kinds in sorted(by_cap.items()):
        for label,lines in kinds.items():dump_parts(lines,f"{cap}-{label}",args.out,index)
    payload={"requested":len(targets),"completed":len(summary),"missing":missing,"heldout_mjd":58484.,
             "selection_sha256":hashlib.sha256(args.selected.read_bytes()).hexdigest(),"matches_sha256":hashlib.sha256(args.matches.read_bytes()).hexdigest(),
             "description":"Per-cap gzip JSON lines. Neighbour records are [source, zero-based piece, ASCII CSV text]. Join by source and piece to recover exact hashed subset bytes.","parts":index}
    (args.out/"coverage.json").write_text(json.dumps(payload,indent=2)+"\n")
    print(payload["completed"],"/",payload["requested"])
if __name__=="__main__":main()
