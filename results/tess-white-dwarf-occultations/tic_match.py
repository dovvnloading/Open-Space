"""Crossmatch the frozen EDR3 catalogue to TIC v8 through epoch-2000 positions."""
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

def load_wds(path):
    lines=[x for x in path.read_text().splitlines() if x and not x.startswith("#")]
    return list(csv.DictReader([lines[0]]+lines[3:],delimiter="\t"))

def crossmatch(row, cache):
    ra,dec,pmra,pmde = [float(row[k]) for k in ("RA_ICRS","DE_ICRS","pmRA","pmDE")]
    cosd = math.cos(math.radians(dec))
    ra0 = (ra - 16*pmra/(3600000*cosd)) % 360
    dec0 = dec - 16*pmde/3600000
    source = row["GaiaEDR3"].strip()
    query = dict(service="Mast.Catalogs.Filtered.Tic.Position.Rows",format="json",
                 params=dict(columns="ID,ra,dec,GAIA,Tmag,Gaiamag,disposition,pmRA,pmDEC",
                             filters=[],ra=ra0,dec=dec0,radius=10/3600))
    url = "https://mast.stsci.edu/api/v0/invoke?"+urlencode({"request":json.dumps(query,separators=(",",":"))})
    path = cache/(source+".json")
    base = dict(gaia=source,name=row["WDJname"].strip(),gmag=float(row["Gmag"]),ra_2016=ra,dec_2016=dec,ra_2000=ra0,dec_2000=dec0)
    rec = dict(gaia=source,url=url)
    try:
        if not path.exists():
            with urlopen(url, timeout=45) as f:
                raw = f.read()
            doc=json.loads(raw)
            if doc.get("status")!="COMPLETE":raise ValueError(str(doc.get("msg") or doc.get("status")))
            path.write_bytes(raw)
        raw=path.read_bytes()
        data=json.loads(raw)["data"]
        contenders=[]
        for t in data:
            sep=3600*math.hypot(((t["ra"]-ra0+180)%360-180)*cosd,t["dec"]-dec0)
            if (t["disposition"] or "").upper() in ("DUPLICATE","ARTIFACT","SPLIT"): continue
            if sep>2: continue
            same=str(t["GAIA"])==source
            if same or (t["Gaiamag"] is not None and abs(t["Gaiamag"]-float(row["Gmag"]))<0.5):
                contenders.append((not same, sep, t))
        contenders.sort(key=lambda t:t[:2])
        matches=[t for t in contenders if t[0]==contenders[0][0]] if contenders else []
        rec.update(bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest(),n_in_cone=len(data))
        if len(matches)==1:
            fallback,sep,t=matches[0]
            return dict(**base,tic=t["ID"],tmag=t["Tmag"],sep_arcsec=round(sep,4),
                        status="position_and_G" if fallback else "source_id_and_position"),rec
        return dict(**base,status="ambiguous" if matches else "unmatched"),rec
    except Exception as exc:
        return dict(**base,status="query_error"),dict(**rec,error=str(exc))

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("parent",type=Path)
    p.add_argument("--cache",type=Path,default=Path("cache/tic"))
    p.add_argument("--out",type=Path,default=Path("record"))
    a=p.parse_args()
    a.cache.mkdir(exist_ok=True,parents=True);a.out.mkdir(exist_ok=True,parents=True)
    with ThreadPoolExecutor(6) as pool:
        results=list(pool.map(lambda r:crossmatch(r,a.cache),load_wds(a.parent)))
    with (a.out/"tic_match.csv").open("w") as f:
        w=csv.DictWriter(f,fieldnames="gaia name gmag ra_2016 dec_2016 ra_2000 dec_2000 tic tmag sep_arcsec status".split());w.writeheader();w.writerows(r for r,_ in results)
    (a.out/"tic_receipts.json").write_text(json.dumps([s for _,s in results],indent=1)+"\n")
    from collections import Counter
    print(Counter(r["status"] for r,_ in results))
