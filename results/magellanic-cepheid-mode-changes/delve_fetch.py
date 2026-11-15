"""Obtain DELVE DR3 MC objects and detections nearest six selected Cepheids."""

from urllib.parse import urlencode
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import csv, math
from follow_fetch import ROOT, TARGETS, get, degrees


def tap(query):
    return "https://datalab.noirlab.edu/tap/sync?" + urlencode(
        dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=query)
    )


def retrieve(target):
    ident = next(
        l
        for l in (ROOT / "inputs/smc4_ident.dat").read_text().splitlines()
        if l.startswith("OGLE-SMC-CEP-" + target + " ")
    ).split()
    ra, dec = 15 * degrees(ident[2]), degrees(ident[3])
    bounds = f"ra BETWEEN {ra-.006} AND {ra+.006} AND dec BETWEEN {dec-.0018} AND {dec+.0018}"
    fname = f"{target}_delve_cone.csv"
    cat = get(
        (
            fname,
            tap(
                f"SELECT TOP 100 objid,ra,dec,gmag,rmag,imag FROM delve_dr3.mc_object WHERE {bounds}"
            ),
        )
    )
    objs = list(csv.DictReader((ROOT / "inputs" / fname).open()))

    def sep(r):
        return 3600 * math.hypot(
            (float(r["ra"]) - ra) * math.cos(math.radians(dec)), float(r["dec"]) - dec
        )

    objs.sort(key=sep)
    if not objs or sep(objs[0]) > 1:
        return [cat]
    id_ = objs[0]["objid"]
    meas = get(
        (
            f"{target}_delve_phot.csv",
            tap(
                f"SELECT mjd,mag,err,filter,exposure,chi,sharp,ra,dec,objid FROM delve_dr3.mc_meas WHERE objid='{id_}'"
            ),
        )
    )
    print(target, id_, round(sep(objs[0]), 2), meas[-1], flush=True)
    return [cat, meas]


if __name__ == "__main__":
    with ThreadPoolExecutor(4) as pool:
        out = [row for result in pool.map(retrieve, TARGETS) for row in result]
    with (ROOT / "delve_manifest.csv").open("w") as f:
        w = csv.writer(f)
        w.writerow(["file", "sha256", "url", "bytes"])
        w.writerows(out)
