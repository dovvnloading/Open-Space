"""Fetch Gaia DR3 and Euclid Q1 neighbours near the expected epoch-dependent positions."""
import argparse
import csv
import hashlib
import io
import json
import math
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parent
TABLES = {
    "gaia": ("gaia_dr3_source",2016.0,8,"source_id,ra,dec,pmra,pmdec,pmra_error,pmdec_error,parallax,parallax_error,phot_g_mean_mag"),
    "euclid": ("euclid_q1_mer_catalogue",2024.25,12,"object_id,ra,dec,point_like_prob,extended_prob,spurious_prob," + ",".join(f"{prefix}_{band}_2fwhm_aper" for band in ("vis","y","j","h") for prefix in ("flux","fluxerr"))),
}
def fetch(job):
    row,kind=job
    table,epoch,radius,cols=TABLES[kind]
    ra0,dec0=float(row["ra2015_deg"]),float(row["dec2015_deg"])
    ra=ra0+(epoch-2015.405)*float(row["pmra"])/(3600*math.cos(math.radians(dec0)))
    dec=dec0+(epoch-2015.405)*float(row["pmdec"])/3600
    query=f"SELECT {cols} FROM {table} WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra:.7f},{dec:.7f},{radius/3600:.9f}))"
    url="https://irsa.ipac.caltech.edu/TAP/sync?"+urlencode(dict(REQUEST="doQuery",LANG="ADQL",FORMAT="csv",QUERY=query,MAXREC=10000))
    path=ROOT/"inputs"/(row["source_name"]+"_"+kind+".csv")
    if not path.exists():
        with urlopen(url,timeout=120) as response:b=response.read()
        if not b.startswith(cols.split(",")[0].encode()+b","):raise ValueError(b[-450:])
        path.write_bytes(b)
    b=path.read_bytes()
    found=list(csv.DictReader(io.StringIO(b.decode())))
    if len(found)>=10000:raise ValueError("row limit")
    for r in found:
        r["separation"]=3600*math.hypot((float(r["ra"])-ra)*math.cos(math.radians(dec)),float(r["dec"])-dec)
    found.sort(key=lambda r:r["separation"])
    print(row["source_name"],kind,len(found),round(found[0]["separation"],2) if found else None,flush=True)
    return dict(source_name=row["source_name"],table=table,epoch=epoch,ra_deg=ra,dec_deg=dec,radius_arcsec=radius,sha256=hashlib.sha256(b).hexdigest(),n=len(found),closest=found[0] if found else None,query=query)
if __name__ == "__main__":
    ap=argparse.ArgumentParser();ap.add_argument("rows",type=Path);ap.add_argument("output",type=Path);ap.add_argument("--gaia-only",action="store_true");args=ap.parse_args()
    (ROOT/"inputs").mkdir(exist_ok=True)
    rows=list(csv.DictReader(args.rows.open()))
    with ThreadPoolExecutor(4) as pool:
        results=list(pool.map(fetch,((r,k) for r in rows for k in (["gaia"] if args.gaia_only else TABLES))))
    args.output.write_text(json.dumps(results,indent=2)+"\n")
