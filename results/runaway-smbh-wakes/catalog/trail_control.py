"""Compare a detached VIS line across the individual Q1 dithers covering it."""
import argparse
import csv
import json
import sys
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'euclid'))
from dither import sky_files,measure_frame
from concurrent.futures import ThreadPoolExecutor

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cache',type=Path,required=True)
    p.add_argument('--measurements',type=Path,default=HERE/'targeted_measurements.jsonl')
    p.add_argument('--id',default='2675970977684150669')
    p.add_argument('--out',type=Path,default=HERE/'dithers')
    a=p.parse_args();a.cache.mkdir(exist_ok=True,parents=True);a.out.mkdir(exist_ok=True,parents=True)
    with a.measurements.open() as f: rows=[json.loads(l) for l in f]
    r=next(r for r in rows if str(r['object_id'])==a.id)
    ra,dec=float(r['ra']),float(r['dec']);pa=np.radians(float(r['mer_position_angle']))
    half=.5*r['forced_length_arcsec']/3600
    off=half*np.array([np.sin(pa)/np.cos(np.radians(dec)),np.cos(pa)])
    ends=np.array([[ra,dec]-off,[ra,dec]+off])
    keys,sia=sky_files(ra,dec,a.cache)
    with ThreadPoolExecutor(2) as pool:results=list(pool.map(lambda k:measure_frame(k,ends,a.cache),keys))
    out=dict(id=a.id,ends_deg=ends.tolist(),sia=sia,frames=results)
    (a.out/(a.id+'.json')).write_text(json.dumps(out,indent=1)+'\n')
    print([(v['status'],v.get('central_ridge_contrast_sigma')) for v in results],flush=True)
if __name__=='__main__':main()
