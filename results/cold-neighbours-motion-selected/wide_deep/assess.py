"""Expand a six-pass association to each 182.625-day W2 epoch bin."""
import argparse,csv,os,sys,json
from pathlib import Path
from collections import defaultdict
import numpy as np
from astropy.io import fits
import gzip,io
BASE=Path(__file__).resolve().parent
sys.path.insert(0,str(BASE.parent/"widecheck"))
import check
def main():
    p=argparse.ArgumentParser();p.add_argument("sources",nargs="+")
    p.add_argument("--catalog",type=Path,default=BASE.parent/"wide"/"selected.csv")
    p.add_argument("--cache",type=Path,default=os.environ.get("UNTIMELY_CACHE"))
    a=p.parse_args();idx=check.fetch(check.INDEX,a.cache)
    table=fits.getdata(io.BytesIO(gzip.decompress(idx)),1);tiles=table["COADD_ID"].astype(str)
    sources=[s for s in csv.DictReader(a.catalog.open()) if s["source_name"] in a.sources]
    sortedgroups=defaultdict(list)
    for s in sources:sortedgroups[s["source_id"][:8]].append(s)
    def run(pair):
        tile,group=pair
        target=table[(tiles==tile)&(table["BAND"]==2)]
        bins=np.floor((target["MJDMEAN"]-55197)/182.625);epochs=[]
        for k in np.unique(bins):
            t=target[bins==k];epochs.append(t[np.argmax(t["N_DET"])])
        out_dir=BASE/tile;out_dir.mkdir(exist_ok=True)
        out=check.tile_job(tile,group,epochs,a.cache,out_dir)
        result=[]
        for s in group:
            r=[d for d in out["detections"] if d["source_name"]==s["source_name"]]
            z=check.score(s,r,len(epochs));result.append(z)
            print(s["source_name"],z["moving"],z["moving_advantage"],z["baseline_yr"],z["post_catwise"],flush=True)
        check.write_csv(out_dir/"summary.csv",result)
        check.write_csv(out_dir/"detections.csv",out["detections"])
        (out_dir/"manifest.json").write_text(json.dumps(out["epochs"],indent=1)+"\n")
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(3) as x:list(x.map(run,sortedgroups.items()))
if __name__=="__main__":main()
