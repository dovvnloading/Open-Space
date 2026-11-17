"""Retrieve Gaia DR3 neighbours around all 2016.0-predicted positions, using chunks of ADQL cones."""
import csv
import hashlib
import io
import json
import math
import sys
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen,Request
from concurrent.futures import ThreadPoolExecutor

ROOT=Path(__file__).resolve().parent
URL="https://gea.esac.esa.int/tap-server/tap/sync"
def run(input_path,output_dir):
    rows=list(csv.DictReader(Path(input_path).open()))
    out=Path(output_dir);out.mkdir(exist_ok=True,parents=True)
    def part(job):
        idx,group=job
        centers=[]
        for r in group:
            d=float(r["dec2015_deg"]);a=float(r["ra2015_deg"])
            a+=0.595*float(r["pmra"])/(3600*math.cos(math.radians(d)))
            d+=0.595*float(r["pmdec"])/3600
            centers.append((r["source_name"],a,d))
        cond=" OR ".join(f"1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{a:.7f},{d:.7f},0.0022222222))" for name,a,d in centers)
        q="SELECT source_id,ra,dec,pmra,pmdec,pmra_error,pmdec_error,parallax,phot_g_mean_mag FROM gaiadr3.gaia_source WHERE "+cond+" ORDER BY source_id"
        url=URL+"?"+urlencode(dict(REQUEST="doQuery",LANG="ADQL",FORMAT="csv",MAXREC=20000,QUERY=q))
        path=out/f"batch{idx:03d}.csv"
        if not path.exists():
            for endpoint in [URL,"https://gaia.aip.de/tap/sync","https://gaia.ari.uni-heidelberg.de/tap/sync",URL]:
                try:
                    data=urlopen(endpoint+"?"+url.split("?",1)[1],timeout=55).read()
                    if data.lstrip(b'\x22').startswith(b"source_id"):break
                except Exception as exc:
                    data=str(exc).encode()
            if not data.lstrip(b'"').startswith(b"source_id"): raise ValueError(data[-550:])
            path.write_bytes(data)
        b=path.read_bytes();source=list(csv.DictReader(io.StringIO(b.decode())))
        if len(source)>=20000:raise ValueError("row limit")
        matches=[]
        for name,a,d in centers:
            for z in source:
                s=3600*math.hypot(((float(z["ra"])-a+180)%360-180)*math.cos(math.radians(d)),float(z["dec"])-d)
                if s<=8:matches.append(dict(source_name=name,separation_arcsec=round(s,4),**z))
        return matches,dict(file=path.name,sha256=hashlib.sha256(b).hexdigest(),query=q,targets=[v[0] for v in centers],n_rows=len(source))
    with ThreadPoolExecutor(3) as ex:data=list(ex.map(part,enumerate([rows[i:i+18] for i in range(0,len(rows),18)])))
    matches=[m for ms,_ in data for m in ms]
    cols=["source_name","separation_arcsec","source_id","ra","dec","pmra","pmdec","pmra_error","pmdec_error","parallax","phot_g_mean_mag"]
    with (out/"matches.csv").open("w") as f:
        w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows(matches)
    (out/"manifest.json").write_text(json.dumps([v for _,v in data],indent=2)+"\n")
    print(len(rows),len(matches))
if __name__=="__main__":run(*sys.argv[1:])
