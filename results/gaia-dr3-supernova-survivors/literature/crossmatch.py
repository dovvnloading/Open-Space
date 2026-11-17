"""Compare each screened source with published objects using full Gaia identifiers."""
import csv
import hashlib
import json
import re
import subprocess
import tempfile
from pathlib import Path
from urllib.request import urlopen

HERE = Path(__file__).resolve().parent
PRIOR = HERE.parents[1]/"dr1-dr19-unbound-stars"/"literature"/"index.csv"
URL = "https://arxiv.org/pdf/2209.09915"
DIGEST = "1e9f68f399578035a5e5fbc94a50e5c7eeb53156c2b4d8baba66fae7a21db58b"


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def main():
    data=urlopen(URL,timeout=80).read()
    if hashlib.sha256(data).hexdigest()!=DIGEST:
        raise ValueError("Igoshev et al. document version changed")
    with tempfile.TemporaryDirectory() as tmp:
        p=Path(tmp)/"paper.pdf";p.write_bytes(data)
        body=subprocess.check_output(["pdftotext","-layout",str(p),"-"],text=True)
    ids=set(re.findall(r"(?<!\d)\d{17,19}(?!\d)",body))
    controls={r["source_id"]:r["name"] for r in read(HERE.parent/"control"/"controls.csv")}
    targets={r["source_id"]:r["verdict"] for r in read(HERE.parent/"control"/"published_targets.csv")}
    prior=read(PRIOR)
    inputs={r["source_id"] for r in read(HERE.parent/"search"/"parent.csv")}
    inputs.update(r["source_id"] for r in read(HERE.parent/"desi"/"tail.csv"))
    inputs.update(r["source_id"] for r in read(HERE.parent/"desi"/"control_matches.csv"))
    out=[]
    for sid in sorted(inputs,key=int):
        sources=sorted({r["catalogue"] for r in prior if r["release"]=="DR3" and r["source_id"]==sid})
        out.append(dict(source_id=sid,ElBadry26_control=controls.get(sid,""),
            ElBadry26_target=targets.get(sid,""),Igoshev22_id_in_text=sid in ids,
            prior_index_dr3_matches="|".join(sources)))
    with (HERE/"matches.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=out[0]);w.writeheader();w.writerows(out)
    receipt=dict(igoshev_url=URL,igoshev_sha256=DIGEST,igoshev_distinct_numeric_ids=len(ids),
                 prior_index_file=str(PRIOR.relative_to(PRIOR.parents[2])),prior_rows=len(prior),
                 prior_sha256=hashlib.sha256(PRIOR.read_bytes()).hexdigest())
    (HERE/"receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")


if __name__=="__main__":
    main()
