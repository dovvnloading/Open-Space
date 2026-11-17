"""Prefilter DR1 main/backup thresholds for cross-program auditing. No spectra used."""
import argparse
import csv
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import tempfile
import fitsio
import numpy as np
from desi_shard import BASE, product_groups, obtain


def roster(row):
    s, p, hp = row['survey'], row['program'], int(row['healpix'])
    _, (folder,names,checks) = product_groups(s,p,hp)
    with tempfile.TemporaryDirectory() as tmp:
        template = Path(tmp) / names[0]
        obtain(BASE + folder + '/' + names[0], template)
        rv = fitsio.read(str(template),ext='RVTAB')
        fib = fitsio.read(str(template),ext='FIBERMAP')
    groups={}
    for i,t in enumerate(rv['TARGETID']):
        groups.setdefault(int(t),[]).append(i)
    rows=[]
    for t,ii in groups.items():
        if len(ii)<2 or not np.all(rv['SUCCESS'][ii] & (rv['RVS_WARN'][ii]==0) & (fib['FIBERSTATUS'][ii]==0)):
            continue
        k=fib[ii[0]]
        rows.append([s,p,hp,t,k['TARGET_RA'],k['TARGET_DEC'],k['REF_ID'],k['REF_CAT'],k['PARALLAX'],
                     '|'.join(map(str,sorted(rv['EXPID'][ii])))])
    return hp,rows


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('selection',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--limit',type=int,default=60)
    args=parser.parse_args()
    fields=[r for r in csv.DictReader(args.selection.open()) if (r['survey'],r['program'])==('main','backup')][:args.limit]
    args.output.mkdir(exist_ok=True,parents=True)
    with ThreadPoolExecutor(max_workers=6) as ex:
        for hp,stars in ex.map(roster,fields):
            with (args.output/f'roster-{hp}.csv').open('w',newline='') as f:
                w=csv.writer(f)
                w.writerow(['survey','program','healpix','targetid','ra','dec','ref_id','ref_cat','parallax','expids'])
                w.writerows(stars)
            print(hp,len(stars),flush=True)
