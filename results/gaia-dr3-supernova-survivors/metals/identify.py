"""Resolve selected identifiers and check DR3 neighbours of four DR2-only rows."""
import csv
import hashlib
import io
import json
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent


def request(stem, endpoint, query, expected):
    url = endpoint + "?" + urllib.parse.urlencode(
        dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=query, MAXREC=200))
    with urllib.request.urlopen(url, timeout=90) as response:
        data, status = response.read(), response.status
    rows = list(csv.DictReader(io.StringIO(data.decode())))
    assert rows and set(rows[0]) == set(expected)
    (HERE/f"{stem}.adql").write_text(query)
    (HERE/f"{stem}.csv").write_bytes(data)
    (HERE/f"{stem}_receipt.json").write_text(json.dumps(dict(
        url=url, utc=datetime.now(timezone.utc).isoformat(), http_status=status,
        rows=len(rows), sha256=hashlib.sha256(data).hexdigest(),
        query_sha256=hashlib.sha256(query.encode()).hexdigest()), indent=2) + "\n")
    return rows


targets = list(csv.DictReader((HERE/"targets.csv").open()))
aliases = ["Gaia DR3 "+r["source_id"] for r in targets]
aliases += ["LP 93-21","SDSS J124043.01+671034.6","WD J110907.96+000132.92"]
query = ("SELECT i.id, b.main_id, b.otype, b.sp_type, b.nbref FROM ident AS i "
         "JOIN basic AS b ON i.oidref=b.oid WHERE i.id IN ("
         + ",".join("'"+a+"'" for a in aliases) + ")")
resolved = request("simbad", "https://simbad.cfa.harvard.edu/simbad/sim-tap/sync", query,
                   ["id","main_id","otype","sp_type","nbref"])
assert {a for a in aliases if a.startswith("Gaia DR3 ")} <= {r["id"] for r in resolved}

dr2 = list(csv.DictReader((HERE/"dr2_only.csv").open()))
ids = [r["designation"].split()[-1] for r in dr2]
query = f"""SELECT n.dr2_source_id, n.dr3_source_id, n.angular_distance, n.magnitude_difference, g.parallax, g.parallax_error, g.pm, g.ruwe
FROM gaiadr3.dr2_neighbourhood AS n JOIN gaiadr3.gaia_source AS g ON n.dr3_source_id = g.source_id
WHERE n.dr2_source_id IN ({", ".join(ids)})"""
rows = request("dr2_check","https://gea.esac.esa.int/tap-server/tap/sync",query,
               ["dr2_source_id","dr3_source_id","angular_distance","magnitude_difference","parallax","parallax_error","pm","ruwe"])
assert {r["dr2_source_id"] for r in rows} == set(ids)
