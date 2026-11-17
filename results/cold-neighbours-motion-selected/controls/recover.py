"""Measure CatWISE recovery of the full UltracoolSheet nearby spectroscopic-Y subset."""
import csv
import hashlib
import math
import sys
from pathlib import Path
from urllib.request import urlopen
from urllib.parse import urlencode
from concurrent.futures import ThreadPoolExecutor
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/"pilot"))
sys.path.insert(0,str(ROOT.parent/"catalogues"))
from fetch_catwise import COLUMNS,ENDPOINT
from select_catwise import select,prefilter,load_corrections
from crosscheck import acquire

def number(v):
    return float(v or "nan")
def main():
    catalogue=acquire()/"UCS"
    controls=[r for r in csv.DictReader(catalogue.open()) if r["spt_ir"].strip().startswith("Y") and number(r["plx_formula"])>=50]
    points=ROOT/"inputs";points.mkdir(exist_ok=True)
    offsets=load_corrections()
    def check(item):
        i,r=item
        dec=number(r["dec_j2000_formula"]);ra=number(r["ra_j2000_formula"])
        pmra=number(r["pmra_formula"])/1000;pmdec=number(r["pmdec_formula"])/1000
        ra+=15.405*pmra/(3600*math.cos(math.radians(dec)));dec+=15.405*pmdec/3600
        cone=f"1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra:.7f},{dec:.7f},0.018))"
        query=f"SELECT {COLUMNS} FROM catwise_2020 WHERE {cone} ORDER BY source_id"
        url=ENDPOINT+"?"+urlencode(dict(REQUEST="doQuery",LANG="ADQL",FORMAT="csv",MAXREC=10000,QUERY=query))
        file=points/f"ctrl{i}.csv"
        if not file.exists():
            b=urlopen(url,timeout=120).read()
            if not b.startswith(b"source_name,"):raise ValueError(b[-350:])
            file.write_bytes(b)
        rows=list(csv.DictReader(file.open()))
        if len(rows)>=10000:raise ValueError("row limit")
        hits=[]
        for row in rows:
            raw_dec=number(row["dec_pm"])
            da0,dd0,mx0,my0=offsets[row["source_id"][:8]]
            d=raw_dec+dd0;a=number(row["ra_pm"])+da0/math.cos(math.radians(d))
            sep=3600*math.hypot(((a-ra+180)%360-180)*math.cos(math.radians(dec)),d-dec)
            # Protect the control identity against a fragment or unrelated crossing trajectory.
            discrep=math.hypot(number(row["pmra"])+mx0-pmra,number(row["pmdec"])+my0-pmdec)
            if sep<=2.5 and discrep<=max(0.20,0.3*math.hypot(pmra,pmdec)):
                hits.append((sep,row))
        hits.sort(key=lambda z:z[0])
        passed=[(d,row,select(row,offsets)) for d,row in hits if select(row,offsets)]
        bright=any(s["bright_tier"] for _,_,s in passed)
        rep=passed[0][:2] if passed else hits[0] if hits else (None,{})
        return dict(control=r["name"],spt_ir=r["spt_ir"],parallax_mas=r["plx_formula"],ra2015=ra,dec2015=dec,
                    pmra=pmra,pmdec=pmdec,parent_rows=len(rows),trajectory_match=int(bool(hits)),
                    query_prefilter=int(any(prefilter(x) for _,x in hits)),screen=int(bool(passed)),bright=int(bright),
                    representative=rep[1].get("source_name",""),separation_arcsec=round(rep[0],3) if rep[0] is not None else "",
                    w2=rep[1].get("w2mpro_pm",""),file=file.name,sha256=hashlib.sha256(file.read_bytes()).hexdigest(),query=query)
    with ThreadPoolExecutor(3) as pool:results=list(pool.map(check,enumerate(controls)))
    with (ROOT/"recovery.csv").open("w") as f:
        w=csv.DictWriter(f,fieldnames=results[0]);w.writeheader();w.writerows(results)
    print(len(results),{k:sum(x[k] for x in results) for k in ("trajectory_match","query_prefilter","screen","bright")})
if __name__=="__main__":main()
