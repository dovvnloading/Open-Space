"""Fetch/replay the complete public Level 3 F1500W mosaics of two repeated fields."""
import argparse
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import sys
from urllib.parse import urlencode
from urllib.request import urlopen

HERE = Path(__file__).resolve().parent
CAMPAIGN = HERE.parents[1] / "jwst-repeat-imaging-transients" / "repeat-catalogs"
sys.path.insert(0, str(CAMPAIGN))
from fetch_catalogs import product_url

ENDPOINT = "https://mast.stsci.edu/api/v0/invoke?"
COLS = "obs_id,target_name,filters,instrument_name,t_min,t_max,proposal_id,dataRights,t_obs_release,obsid,calib_level,dataURL,s_ra,s_dec"
FILTERS = {"obs_collection":["JWST"], "filters":["F1500W"], "calib_level":[3],
           "proposal_id":["1911","4857"], "dataRights":["PUBLIC"]}
QUERY = {"service":"Mast.Caom.Filtered", "params":{"columns":COLS, "filters":[
    {"paramName":k,"values":v} for k,v in FILTERS.items()]},
         "format":"json","page":1,"pagesize":2000}

def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()

def download(row, cache):
    dest = cache / row["filename"]
    if not dest.exists():
        temp = dest.with_suffix(".part")
        with urlopen(row["url"], timeout=120) as response, temp.open("wb") as f:
            while block := response.read(1024*1024):
                f.write(block)
        temp.replace(dest)
    sha = digest(dest)
    if row.get("sha256", sha) != sha:
        raise ValueError(f"Frozen input mismatch: {dest}")
    return dict(row, bytes=dest.stat().st_size, sha256=sha)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache", type=Path)
    parser.add_argument("--freeze", action="store_true", help="Create a new selection, without overwriting frozen receipts.")
    arg=parser.parse_args()
    arg.cache.mkdir(parents=True, exist_ok=True)
    receipt=HERE/"inputs.json"
    if arg.freeze:
        if receipt.exists():
            raise FileExistsError(receipt)
        url=ENDPOINT+urlencode({"request":json.dumps(QUERY)})
        raw=urlopen(url,timeout=90).read()
        archive=json.loads(raw)
        if archive["status"] != "COMPLETE" or len(archive["data"]) >= QUERY["pagesize"]:
            raise ValueError("Incomplete MAST result")
        (HERE/"mast_response.json").write_bytes(raw)
        rows=[]
        for r in archive["data"]:
            target=r["target_name"].upper().replace("-","").replace(" ","")
            # GO 4857 target names use WD and occasionally different punctuation.
            if not ("1202" in target or "2105" in target or target in ("LAWD83","LP8527")):
                continue
            filename=r["dataURL"].rsplit("/",1)[1]
            if not filename.endswith("_i2d.fits"):
                raise ValueError(r["dataURL"])
            rows.append(dict(r,filename=filename,url=product_url(r["obs_id"],"i2d.fits")))
        if len(rows)!=4:
            raise ValueError([(r["target_name"],r["obs_id"]) for r in archive["data"]])
        rows.sort(key=lambda r:(r["target_name"],r["t_min"]))
        meta=dict(query_url=url,response_sha256=hashlib.sha256(raw).hexdigest())
    else:
        meta=json.loads(receipt.read_text())
        rows=meta["inputs"]
    with ThreadPoolExecutor(4) as pool:
        completed=list(pool.map(lambda r:download(r,arg.cache),rows))
    for r in completed:
        print(r["target_name"], r["filename"], r["bytes"], r["sha256"])
    if arg.freeze:
        receipt.write_text(json.dumps(dict(meta,inputs=completed),indent=2)+"\n")

if __name__=="__main__":
    main()
