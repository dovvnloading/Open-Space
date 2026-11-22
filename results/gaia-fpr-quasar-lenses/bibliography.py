"""Retrieve object-linked references and aliases for individual novelty review."""
import csv
import hashlib
import io
import json
import re
import sys
from pathlib import Path

import requests

HOSTS = ["https://simbad.cfa.harvard.edu/simbad/sim-tap/sync",
         "https://simbad.cds.unistra.fr/simbad/sim-tap/sync",
         "https://simbad.idia.ac.za/simbad/sim-tap/sync"]
TOPIC = re.compile(r"lens|quasar pair|dual (?:quasar|AGN)|double quasar|binary quasar|multiplet|varstrometr", re.I)


def query(statement, columns):
    errors = []
    for service in HOSTS:
        try:
            res = requests.get(service, params={"REQUEST":"doQuery","LANG":"ADQL","FORMAT":"csv",
                                                "MAXREC":20000,"QUERY":statement},timeout=(8,45))
            res.raise_for_status()
            reader = csv.DictReader(io.StringIO(res.text))
            if reader.fieldnames is None or set(reader.fieldnames) != set(columns):
                raise ValueError("Unexpected response schema: "+str(reader.fieldnames)[:100])
            rows = list(reader)
            if len(rows)>=20000:
                raise ValueError("Output may be truncated")
            return rows,dict(service=service,query=statement,bytes=len(res.content),
                             sha256=hashlib.sha256(res.content).hexdigest(),failures=errors)
        except (ValueError,requests.RequestException) as exc:
            errors.append(dict(service=service,error=str(exc)[:200]))
    return [],dict(query=statement,failures=errors)


def inspect(ident, ra, dec):
    region = f" WHERE CONTAINS(POINT('ICRS',b.ra,b.dec),CIRCLE('ICRS',{ra:.10f},{dec:.10f},0.00277777778))=1"
    aliases,raudit = query("SELECT b.main_id,b.otype,b.ra,b.dec,i.id FROM basic AS b JOIN ident AS i ON b.oid=i.oidref"+region,
                           ["main_id","otype","ra","dec","id"])
    references,baudit = query("SELECT DISTINCT b.main_id,r.bibcode,r.title,r.doi FROM basic AS b "
                         "JOIN has_ref AS h ON b.oid=h.oidref JOIN ref AS r ON h.oidbibref=r.oidbib"+region,
                         ["main_id","bibcode","title","doi"])
    return dict(system_id=ident,ra=ra,dec=dec,aliases=aliases,references=references,
                topic_titles=[r for r in references if TOPIC.search(r["title"])],receipts=[raudit,baudit])


def run(evidence, dest):
    dest.mkdir(exist_ok=True,parents=True)
    with (evidence/"system_dispositions.csv").open() as f:
        dispositions={r["system_id"]:r for r in csv.DictReader(f)}
    # Examine all unpublished highly ranked systems under the main selection.
    positions = {}
    for line in (evidence/"review_packets.jsonl").read_text().splitlines():
        packet = json.loads(line)
        row = dispositions[packet["system_id"]]
        if row["prior"]=="False" and row["full_morphology"]=="True":
            ident = row["representative"]
            positions[ident]=packet["imaging"]["lens_model"]["center_radec_deg"]
    # Additional spectroscopic follow-up can request the same object-linked check.
    targets = dest.parent/"bibliography_targets.csv"
    if targets.exists():
        with targets.open() as f:
            for row in csv.DictReader(f):
                positions[row["system_id"]]=[float(row["ra"]), float(row["dec"])]
    result = [inspect(i,*pos) for i,pos in sorted(positions.items())]
    (dest/"individual_literature.json").write_text(json.dumps(result,indent=2)+"\n")


if __name__=="__main__":
    run(Path(sys.argv[1]),Path(sys.argv[2]))
