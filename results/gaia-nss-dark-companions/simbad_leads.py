"""Resolve all explicitly discussed Gaia NSS leads by exact DR3 ID in SIMBAD.

Does not infer novelty from a missing SIMBAD row. Stores the literal response and
SHA-256 receipt in a chosen output directory, never rewriting committed inputs.
Run: python simbad_leads.py --out /tmp/simbad-nss-referee
"""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

LEADS = [
    # The twelve repeat-observed above-three-solar-mass screening entries.
    "182029303939094784", "275582835635329024", "422991301559494528",
    "460498563758448512", "878555832642451968", "1594082579405063168",
    "2224573745414713984", "2951241529971223808", "3372818683473351936",
    "3664684869697065984", "3998734261904478848", "5713565792393224960",
    # Other explicitly discussed follow-ups, luminous negative controls and benchmarks.
    "3378588057203660160", "1830677037678583680", "4661290764764683776",
    "4685960537144695552", "902085660651149440", "2077800002648034944",
    "2130638481326861824", "4658457941738020736",
    "4373465352415301632", "5870569352746779008", "6328149636482597888",
]
TAP = "https://simbad.cfa.harvard.edu/simbad/sim-tap/sync?"
QUERY = ("SELECT a.id,b.main_id,b.otype,b.sp_type,b.ra,b.dec "
         "FROM ident AS a JOIN basic AS b ON b.oid=a.oidref "
         "WHERE a.id IN ("+",".join("'Gaia DR3 "+sid+"'" for sid in LEADS)+")")
URL = TAP+urlencode({"REQUEST":"doQuery","LANG":"ADQL","FORMAT":"csv","QUERY":QUERY})


def main(out):
    out.mkdir(parents=True,exist_ok=True)
    with urlopen(URL,timeout=180) as response:
        raw=response.read()
    text=raw.decode("utf-8-sig")
    rows=list(csv.DictReader(io.StringIO(text)))
    assert "main_id" in text.splitlines()[0] and "id" in text.splitlines()[0], "Unexpected SIMBAD CSV"
    requested={"Gaia DR3 "+x for x in LEADS}
    assert all(r["id"] in requested for r in rows), "SIMBAD returned an unrequested source"
    assert len({r["id"] for r in rows})==len(rows), "Duplicate SIMBAD ID"
    (out/"lead_simbad.csv").write_bytes(raw)
    manifest={"query":QUERY,"url":URL,"bytes":len(raw),
              "sha256":hashlib.sha256(raw).hexdigest(),
              "requested_source_ids":LEADS,
              "returned_matches":len(rows),
              "unmatched_exact_ids":[sid for sid in LEADS if "Gaia DR3 "+sid not in {r["id"] for r in rows}]}
    (out/"lead_simbad_receipt.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print("SIMBAD QUERY RECEIPT:",json.dumps({k:manifest[k] for k in ("url","bytes","sha256","returned_matches","unmatched_exact_ids")},indent=2))
    for row in rows:
        print("SIMBAD EXACT-ID:",json.dumps(row,ensure_ascii=False))
    print("No exact-ID match is not an assertion of novelty or of source absence.")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",required=True,type=Path)
    main(p.parse_args().out)
