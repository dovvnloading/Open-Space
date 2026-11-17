"""Run epoch-correct Gaia queries on complete caps; results remain valid as new caps arrive."""
import csv
import hashlib
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from collections import defaultdict
from gaia_batch import run
ROOT=Path(__file__).resolve().parent
def main():
    inputs=ROOT.parent/"wide"/"selected.csv";out=ROOT/"gaia_wide";out.mkdir(exist_ok=True)
    groups=defaultdict(list)
    for r in csv.DictReader(inputs.open()):groups[r["field"]].append(r)
    def field(item):
        name,data=item;path=out/name;path.mkdir(exist_ok=True)
        infile=path/"targets.csv"
        keys=sorted((r["source_name"],r["ra2015_deg"],r["dec2015_deg"],r["pmra"],r["pmdec"]) for r in data)
        signature=hashlib.sha256(json.dumps(keys).encode()).hexdigest()
        with infile.open("w") as f:
            w=csv.DictWriter(f,fieldnames=data[0]);w.writeheader();w.writerows(data)
        if (path/"selection_sha256").exists() and (path/"selection_sha256").read_text().strip()!=signature:
            raise ValueError("Selection changed for "+name+"; replay in a fresh field directory")
        if not (path/"manifest.json").exists():run(infile,path)
        (path/"selection_sha256").write_text(signature+"\n")
        return list(csv.DictReader((path/"matches.csv").open())),json.loads((path/"manifest.json").read_text())
    allrows=[];allmeta=[]
    with ThreadPoolExecutor(2) as pool:
        for matches,meta in pool.map(field,sorted(groups.items())):allrows.extend(matches);allmeta.extend(meta)
    with (out/"matches.csv").open("w") as f:
        w=csv.DictWriter(f,fieldnames=["source_name","separation_arcsec","source_id","ra","dec","pmra","pmdec","pmra_error","pmdec_error","parallax","phot_g_mean_mag"]);w.writeheader();w.writerows(allrows)
    (out/"manifest.json").write_text(json.dumps(allmeta,indent=1)+"\n")
    print(len(groups),sum(map(len,groups.values())),len(allrows),flush=True)
if __name__=="__main__":main()
