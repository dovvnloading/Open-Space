from follow_fetch import get, ROOT, degrees
from urllib.parse import urlencode
import csv

r = next(
    l.split()
    for l in (ROOT / "inputs/smc4_ident.dat").read_text().splitlines()
    if l.startswith("OGLE-SMC-CEP-2332 ")
)
ra = 15 * degrees(r[2])
dec = degrees(r[3])


def query(sql):
    return "https://simbad.cfa.harvard.edu/simbad/sim-tap/sync?" + urlencode(
        dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=sql)
    )


sql = f"SELECT oid,main_id,ra,dec,otype,nbref FROM basic WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{dec},0.00055556))"
receipts = [get(("2332_simbad.csv", query(sql)))]
objs = list(csv.DictReader((ROOT / "inputs/2332_simbad.csv").open()))
for obj in objs:
    oid = int(obj["oid"])
    receipts.append(
        get(
            (
                "2332_references.csv",
                query(
                    f"SELECT bibcode,title,doi FROM has_ref JOIN ref ON has_ref.oidbibref=ref.oidbib WHERE has_ref.oidref={oid}"
                ),
            )
        )
    )
    receipts.append(
        get(("2332_identifiers.csv", query(f"SELECT id FROM ident WHERE oidref={oid}")))
    )
with (ROOT / "literature_manifest.csv").open("w") as f:
    w = csv.writer(f)
    w.writerow(["file", "sha256", "url", "bytes"])
    w.writerows(receipts)
print(objs)
