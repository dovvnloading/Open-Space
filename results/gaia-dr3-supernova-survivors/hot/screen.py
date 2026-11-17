"""Screen Culpan et al. known hot subdwarfs by survey RV with and without error margins."""
import csv
import hashlib
import json
import pathlib
import os
import sys
import urllib.parse
import urllib.request
from datetime import datetime, timezone

HERE = pathlib.Path(__file__).resolve().parent
CACHE = pathlib.Path(os.environ.get("SURVIVOR_CACHE", "cache"))
VIZ = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv"
CUTOFF = 700.
Z = 1.6448536269514722
FIELDS = "recno,GaiaEDR3,Name,RA_ICRS,DE_ICRS,SpClass,RVSDSS,e_RVSDSS,RVLAMOST,e_RVLAMOST"


def download():
    params = {"-source":"J/A+A/662/A40/knownhsd", "-out":FIELDS, "-out.max":"unlimited"}
    url = VIZ + "?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout=120) as response:
        data = response.read()
        report = dict(url=url, retrieved_utc=datetime.now(timezone.utc).isoformat(), status=response.status,
                      server_date=response.getheader("Date"), nbytes=len(data), sha256=hashlib.sha256(data).hexdigest())
    CACHE.mkdir(parents=True, exist_ok=True)
    (CACHE / "culpan_knownhsd.tsv").write_bytes(data)
    (HERE / "catalog_receipt.json").write_text(json.dumps(report, indent=2)+"\n")


def catalogue():
    output = (CACHE/"culpan_knownhsd.tsv").read_text().splitlines()
    i = next(i for i,s in enumerate(output) if s.startswith("recno\t"))
    lines = [output[i], *[r for r in output[i+3:] if r.strip() and not r.startswith("#")]]
    return [{k:v.strip() for k,v in r.items()} for r in csv.DictReader(lines, delimiter="\t")]


def analyse():
    sky = catalogue()
    n_all = len(sky)
    assert {int(a["recno"]) for a in sky} == set(range(1,6617))
    rv = [a for a in sky if a["RVSDSS"] or a["RVLAMOST"]]
    n_rv = len(rv)
    scored = []
    for r in rv:
        signed = [float(r[k]) for k in ("RVSDSS", "RVLAMOST") if r[k]]
        lower = [abs(float(r[k]))-Z*float(r["e_"+k]) for k in ("RVSDSS","RVLAMOST") if r[k] and r["e_"+k] and float(r["e_"+k])>0]
        scored.append({**r, "max_abs_rv":max(map(abs,signed)), "max_abs_rv_minus_1645err":max(lower,default=None)})
    scored.sort(key=lambda r: -r["max_abs_rv"])
    with (HERE / "tail.csv").open("w") as f:
        w=csv.DictWriter(f, fieldnames=scored[0], lineterminator="\n");w.writeheader()
        w.writerows(r for r in scored if r["max_abs_rv"]>=500)
    report = dict(catalogue=n_all, with_rv=n_rv, nonnull_sdss=sum(bool(r["RVSDSS"]) for r in rv),
                  nonnull_lamost=sum(bool(r["RVLAMOST"]) for r in rv), both=sum(bool(r["RVLAMOST"] and r["RVSDSS"]) for r in rv),
                  rows_with_gaia_id=sum(bool(r["GaiaEDR3"]) for r in sky),
                  distinct_gaia_ids=len({r["GaiaEDR3"] for r in sky if r["GaiaEDR3"]}),
                  rv_rows_without_gaia_id=sum(not r["GaiaEDR3"] for r in rv), threshold_kms=CUTOFF,
                  point_ge_700=sum(r["max_abs_rv"]>=CUTOFF for r in scored),
                  lower_ge_700=sum((r["max_abs_rv_minus_1645err"] or -1)>=CUTOFF for r in scored),
                  point_ge_500=sum(r["max_abs_rv"]>=500 for r in scored))
    (HERE / "summary.json").write_text(json.dumps(report, indent=2)+"\n")
    print(report)

if __name__ == "__main__":
    if "--download" in sys.argv: download()
    analyse()
