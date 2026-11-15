"""Descriptive phase specificity and photometry-level injection recovery."""
import argparse
import csv
from pathlib import Path
import sys

import numpy as np

HERE=Path(__file__).resolve().parent



def main(out,locked):
    sys.path.insert(0,str(locked.resolve()))
    from holdout import quantify
    rows=[]
    with (locked/'phase_lock.csv').open() as f:
        locks=list(csv.DictReader(f))
    for lock in locks:
        number=lock['mpc'];p=float(lock['period_h']);duration=float(lock['duration_h'])
        center=float(lock['center_phase']);origin=float(lock['training_jd_start'])
        data=np.loadtxt(HERE/f'photometry_{number}.csv.gz',delimiter=',',skiprows=1)
        original,events=quantify(data,lock,8)
        # Alternating cyclic shifts preserve the timestamps and observed
        # residuals. Adjacent phases are correlated, so the rank is descriptive.
        phases=(center+np.arange(64)/64)%1
        ds=[]
        for phase in phases:
            tmp=dict(lock,center_phase=str(phase))
            ds.append(quantify(data,tmp,8)[0][9])
        row=dict(mpc=number,frozen_phase=center,unrestricted_best_phase=phases[np.argmax(ds)],
                 equally_deep_or_deeper_shifts=int(np.sum(np.asarray(ds)>=ds[0])),
                 n_offsets=len(ds))
        for amplitude in [.05,.10]:
            injected=data.copy()
            inside=abs(((data[:,0]-origin)*24/p-center+.5)%1-.5)<duration/2/p
            injected[inside,8] *= 10**(-.4*amplitude)
            s,_=quantify(injected,lock,8)
            row[f'increment_for_{int(100*amplitude)}_centimag']=s[9]-original[9]
        rows.append(row)
    out.mkdir(parents=True,exist_ok=True)
    with (out/'sensitivity.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,default=HERE)
    p.add_argument("--locked-code",type=Path,default=HERE/"locked")
    args=p.parse_args()
    main(args.out,args.locked_code)
