"""Obtain a bright all-sky Gaia EDR3 roster, geometric sector footprints, and archive inputs."""
import argparse
import csv
import hashlib
import io
import json
import time
import zipfile
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from tess_stars2px import tess_stars2px_function_entry

VIZIER = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?"
TESSCUT = "https://mast.stsci.edu/tesscut/api/v0.1/astrocut?"
FIELDS = "WDJname,GaiaEDR3,RA_ICRS,DE_ICRS,pmRA,pmDE,Gmag,BPmag,RPmag,Pwd,MassH,loggH"
GAIA_FIELDS = "Source,RA_ICRS,DE_ICRS,pmRA,pmDE,Gmag,BPmag,RPmag"
TIC_FIELDS = "TIC,RAJ2000,DEJ2000,Tmag,GAIA,pmRA,pmDE,Disp"

def vizier(table, cols, **extra):
    return VIZIER + urlencode({"-source":table, "-out":cols, "-out.max":"unlimited", **extra})

def tsv_rows(path):
    lines = [s for s in path.read_text().splitlines() if s.strip() and not s.startswith("#")]
    if len(lines)<3:
        return []
    table = list(csv.reader(lines,delimiter="\t"))
    return [dict(zip(table[0],map(str.strip,row))) for row in table[3:]]

def propagated(row, year):
    dec = float(row["DE_ICRS"])
    dt = (year - 2016.0)/3.6e6
    ra = float(row["RA_ICRS"]) + dt*float(row["pmRA"] or 0)/np.cos(np.radians(dec))
    return ra % 360, dec + dt*float(row["pmDE"] or 0)

def download(path, url, force=False):
    receipt = path.with_suffix(path.suffix+".json")
    if path.is_file() and receipt.is_file() and not force:
        return
    meta = dict(url=url, retrieved_utc=datetime.now(timezone.utc).isoformat())
    for attempt in range(3):
        try:
            with urlopen(url,timeout=170) as response:
                raw=response.read()
                meta["http_status"]=response.status
            if path.suffix==".fits":
                if raw[:1]==b"{":
                    meta["server_response"]=json.loads(raw)
                    raise ValueError("TESScut returned a status message instead of FITS")
                with zipfile.ZipFile(io.BytesIO(raw)) as archive:
                    names=[n for n in archive.namelist() if n.lower().endswith(".fits")]
                    if len(names)!=1:
                        raise ValueError(f"Expected one FITS member, got {names}")
                    raw=archive.read(names[0])
                    meta["member"]=names[0]
            elif path.suffix==".tsv" and (b"#INFO" not in raw or b"#INFO\tError" in raw):
                raise ValueError("VizieR error or incomplete response")
            meta.update(status="ok",bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest())
            meta.pop("error",None)
            path.write_bytes(raw)
            break
        except Exception as exc:
            meta.update(status="error",error=str(exc))
            if attempt<2:
                time.sleep(2**attempt)
    receipt.write_text(json.dumps(meta,indent=2)+"\n")
    print(path.name,meta["status"],meta.get("bytes",meta.get("error")),flush=True)

def footprint(stars,sector):
    # Explicit fixed coordinate epoch for Oct-Nov 2026, the two selected sectors.
    coords=np.array([propagated(row,2026.8) for row in stars])
    out=tess_stars2px_function_entry(np.arange(len(stars)),*coords.T,trySector=sector)
    rows=[]
    for index,camera,ccd,x,y in zip(out[0],*out[4:8]):
        if index<0:
            continue
        rows.append(dict(stars[index],sector=sector,ra_obs=coords[index,0],
                         dec_obs=coords[index,1],camera=camera,ccd=ccd,column=x,row=y))
    return rows

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out",type=Path)
    parser.add_argument("--sectors",nargs="+",type=int,default=[109,110],choices=[109,110])
    parser.add_argument("--workers",type=int,default=4)
    args=parser.parse_args()
    args.out.mkdir(exist_ok=True,parents=True)
    cat=args.out/"bright.tsv"
    download(cat,vizier("J/MNRAS/508/3877/maincat",FIELDS,Pwd=">0.75",Gmag="<14.0",**{"-sort":"Gmag"}))
    jobs={}
    for sector in args.sectors:
        rows=footprint(tsv_rows(cat),sector)
        with (args.out/f"s{sector}.csv").open("w",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=[*FIELDS.split(","),"sector","ra_obs","dec_obs","camera","ccd","column","row"])
            writer.writeheader()
            writer.writerows(rows)
        for row in rows:
            source=row["GaiaEDR3"]
            jobs[args.out/f"{source}-s{sector}.fits"]=TESSCUT+urlencode(dict(ra=f"{row['ra_obs']:.8f}",dec=f"{row['dec_obs']:.8f}",sector=sector,x=9,y=9,units="px"))
            coord=f"{row['RA_ICRS']} {row['DE_ICRS']}"
            jobs[args.out/f"{source}-gaia.tsv"]=vizier("I/355/gaiadr3", GAIA_FIELDS,
                                 **{"-c":coord, "-c.rs":"180", "-c.eq":"J2000"})
            ra,dec=propagated(row,2000.0)
            jobs[args.out/f"{source}-tic82.tsv"]=vizier("IV/39/tic82",TIC_FIELDS,
                                 **{"-c":f"{ra:.9f} {dec:+.9f}", "-c.rs":"30"})
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        list(pool.map(lambda item:download(*item), jobs.items()))

if __name__=="__main__":
    main()
