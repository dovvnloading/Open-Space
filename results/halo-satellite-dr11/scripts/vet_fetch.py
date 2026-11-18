"""Freeze independent astrometry, variables, and known-object/contaminant lists."""
import argparse
import concurrent.futures
import csv
import datetime
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import re
import time
from urllib.parse import urlencode
from urllib.request import Request, urlopen

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "vet"
FIELDS = {"bootes": (209,212,13.5,15.5), "coma": (185.75,188.75,22.9,24.9)}
VIZ = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?"
GAIA = "https://gaia.ari.uni-heidelberg.de/tap/sync?"
LVDB = "https://github.com/apace7/local_volume_database/releases/download/v1.1.2/comb_all.csv"
HARRIS = "https://www.physics.mcmaster.ca/~harris/mwgc.dat"

def sha(b):
    return hashlib.sha256(b).hexdigest()

def tabular(raw):
    lines = [x for x in raw.decode().splitlines() if x.strip() and not x.startswith("#")]
    if len(lines) < 3 or not lines[2].startswith("---"):
        raise ValueError("VizieR did not return a valid TSV table")
    return [{k:v.strip() for k,v in r.items()} for r in csv.DictReader([lines[0]]+lines[3:], delimiter="\t")]

def local_sky(ra,dec):
    # Retain a two-degree coordinate margin for possible extended neighbours.
    return any(lo-2<=ra<=hi+2 and south-2<=dec<=north+2 for lo,hi,south,north in FIELDS.values())

def reference_rows(key, raw):
    if key=="lvdb":
        rows = list(csv.DictReader(io.StringIO(raw.decode())))
        assert len(rows)>1000 and all(k in rows[0] for k in ('ra','dec','rhalf','table'))
        cols = ("key","ra","dec","name","host","table","confirmed_real","confirmed_galaxy","rhalf","distance_modulus","pmra","pmdec","ref_proper_motion")
        return [dict((k,r[k]) for k in cols) for r in rows if r["ra"] and r["dec"] and local_sky(float(r["ra"]),float(r["dec"]))],len(rows),cols
    if key=="harris":
        text = raw.decode().split("Part I:")[1].split("Part II:")[0]
        result = []
        pat = r'^.{25}(\d\d) (\d\d) ([\d.]+)\s+([+-])(\d\d) (\d\d) ([\d.]+)'
        count = 0
        for line in text.splitlines():
            m = re.search(pat,line)
            if not m: continue
            h,mn,s = map(float, m.groups()[:3])
            d,am,asec = map(float, m.groups()[4:])
            ra,dec = 15*(h+mn/60+s/3600), (d+am/60+asec/3600)*(1 if m[4]=="+" else -1)
            count += 1
            if local_sky(ra,dec): result.append(dict(name=line[:12].strip(),ra=ra,dec=dec))
        assert count==157, f"parsed {count} of 157 Harris positions"
        return result,count,("name","ra","dec")
    rows = tabular(raw)
    cols = ("Name","RA_ICRS","DE_ICRS","r50","dist50")
    assert len(rows)>5000 and all(k in rows[0] for k in cols)
    return [{k:r[k] for k in cols} for r in rows if local_sky(float(r["RA_ICRS"]), float(r["DE_ICRS"]))],len(rows),cols

def fetch(key, url, kind):
    rec_path = OUT / "receipts" / (key+".json")
    out = OUT/(key+".csv.gz")
    if rec_path.exists() and out.exists() and sha(out.read_bytes())==json.loads(rec_path.read_text())["output_sha256"]:
        return
    for i in range(4):
        try:
            with urlopen(Request(url,headers={'User-Agent':'Open-Space astronomy research'}),timeout=180) as response:
                raw = response.read()
                server_date = response.headers.get('Date')
            if kind=="gaia":
                assert raw.startswith(b"source_id,"), raw[:300]
                n=len(raw.splitlines())-1
                assert 0<=n<150000
                output, nret = raw,n
            else:
                if kind=='reference':
                    rows,n,cols = reference_rows(key,raw)
                elif kind=='rr':
                    names={"Source":"source_id","RA_ICRS":"ra","DE_ICRS":"dec","Class":"best_classification","Gmagavg":"int_average_g","AG":"g_absorption","PF":"pf","P1O":"p1_o"}
                    rows=[{v:r[k] for k,v in names.items()} for r in tabular(raw)]
                    n,cols=len(rows),tuple(names.values())
                else:
                    rows = tabular(raw)
                    cols = tuple(rows[0]) if rows else ("id","ra","dec")
                    n=len(rows)
                stream=io.StringIO()
                w=csv.DictWriter(stream,fieldnames=cols,lineterminator="\n")
                w.writeheader()
                w.writerows(rows)
                nret,output=len(rows),stream.getvalue().encode()
            enc = gzip.compress(output,mtime=0)
            temp = out.with_suffix(".tmp")
            temp.write_bytes(enc)
            temp.replace(out)
            rec = dict(url=url, checked_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),response_date=server_date,response_bytes=len(raw),response_sha256=sha(raw),response_rows=n,retained_rows=nret,output_sha256=sha(enc),output=out.name)
            rec_path.write_text(json.dumps(rec,indent=2)+"\n")
            return
        except Exception:
            if i==3:raise
            time.sleep(2**i)

def sources():
    yield "lvdb",LVDB,"reference"
    yield "harris",HARRIS,"reference"
    for y,table in [(2024,"686/A42"),(2023,"673/A114")]:
        p={"-source":"J/A+A/"+table+"/clusters","-out":"Name,RA_ICRS,DE_ICRS,r50,dist50","-out.max":"unlimited"}
        yield f"hunt_reffert_{y}",VIZ+urlencode(p),"reference"
    for field,(lo,hi,south,north) in FIELDS.items():
        # Small padding also covers nearest bright stars outside footprint.
        bounds=(f"1=CONTAINS(POINT('ICRS',g.ra,g.dec), CIRCLE('ICRS',{(lo+hi)/2},{(south+north)/2},2.5)) AND " +
                f"g.ra BETWEEN {lo-.3} AND {hi+.3} AND g.dec BETWEEN {south-.3} AND {north+.3}")
        cols="g.source_id, g.ra, g.dec, g.pmra, g.pmdec, g.pmra_error, g.pmdec_error, g.pmra_pmdec_corr, g.parallax, g.parallax_error, g.ruwe, g.phot_g_mean_mag, g.bp_rp"
        query=f"SELECT {cols} FROM gaiadr3.gaia_source AS g WHERE {bounds}"
        yield f"gaia_{field}",GAIA+urlencode(dict(REQUEST="doQuery",LANG="ADQL",FORMAT="csv",MAXREC=150000,QUERY=query)),"gaia"
        center=f"{(lo+hi)/2} {(south+north)/2}"
        params={"-source":"I/358/vrrlyr","-out":"Source,RA_ICRS,DE_ICRS,Class,Gmagavg,AG,PF,P1O","-c":center,"-c.rd":3,"-out.max":"unlimited"}
        yield f"rrlyrae_{field}",VIZ+urlencode(params),"rr"
        for table,cols in [("table3","WHL,RAJ2000,DEJ2000,zph,zsp"),("table4","WH,RAJ2000,DEJ2000,zsp")]:
            params={"-source":"J/ApJ/807/178/"+table,"-out":cols,"-c":f"{(lo+hi)/2} {(south+north)/2}","-c.rd":3,"-out.max":"unlimited"}
            yield f"whl_{field}_{table}",VIZ+urlencode(params),"whl"

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("keys",nargs="*")
    args=ap.parse_args()
    (OUT/"receipts").mkdir(parents=True,exist_ok=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        jobs={pool.submit(fetch,*s):s[0] for s in sources() if not args.keys or s[0] in args.keys}
        for job in concurrent.futures.as_completed(jobs):
            try: job.result()
            except Exception as err:
                print(jobs[job],repr(err),flush=True)
                raise
