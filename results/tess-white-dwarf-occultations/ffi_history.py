"""Retrieve all earlier Cycle 5+ TESScut visits available for the frozen S109 accepted stars."""
import argparse
import csv
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode

from ffi_acquire import TESSCUT,download,propagated,footprint

def main():
    from tess_stars2px import TESS_Spacecraft_Pointing_Data
    from astropy.time import Time
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input",type=Path)
    parser.add_argument("visits",type=Path)
    args=parser.parse_args()
    accepted={r["target"] for r in csv.DictReader(args.visits.open())}
    yield_rows=[]
    jobs={}
    for row in csv.DictReader((args.input/"s109.csv").open()):
        source=row["GaiaEDR3"]
        if source not in accepted: continue
        url="https://mast.stsci.edu/tesscut/api/v0.1/sector?"+urlencode(dict(ra=row["ra_obs"],dec=row["dec_obs"]))
        path=args.input/f"{source}-sectors.json"
        download(path,url)
        sections=json.loads(path.read_text())["results"]
        for hit in sections:
            sector=int(hit["sector"])
            if not 56<=sector<109:continue
            epoch=float(Time(TESS_Spacecraft_Pointing_Data(trySector=sector).midtimes[0],format="jd").jyear)
            ra,dec=propagated(row,epoch)
            item=dict(row,sector=sector,camera=hit["camera"],ccd=hit["ccd"],ra_obs=ra,dec_obs=dec,
                      column="",row="",coordinate_epoch=epoch)
            yield_rows.append(item)
            jobs[args.input/f"{source}-s{sector}.fits"]=TESSCUT+urlencode(
                 dict(ra=f"{ra:.8f}",dec=f"{dec:.8f}",sector=sector,x=9,y=9,units="px"))
    with (args.input/"history.csv").open("w",newline="") as file:
        writer=csv.DictWriter(file,fieldnames=yield_rows[0].keys())
        writer.writeheader()
        writer.writerows(yield_rows)
    with ThreadPoolExecutor(4) as workers:
        list(workers.map(lambda x:download(*x),jobs.items()))

if __name__=="__main__":
    main()
