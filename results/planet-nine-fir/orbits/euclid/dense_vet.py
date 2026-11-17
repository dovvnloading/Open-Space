"""Apply independent stationary-source mask and single-exposure removals to every peak."""
import csv
import json
from pathlib import Path
import numpy as np
from dense import BASE,HERE,design,score_tile
from offgrid import catalog


def main():
    rows=list(csv.DictReader((HERE/'dense_triggers.csv').open()))
    names,t,bank,m=design()
    data=[np.load(BASE/(e+'_data.npy'),mmap_mode='r') for e in names]
    weights=[np.load(BASE/(e+'_weight.npy'),mmap_mode='r') for e in names]
    tree=catalog()
    out=[]
    for row in rows:
        x=int(row['x']);y=int(row['y'])
        dist,nn=tree.query([x,y])
        rng=np.s_[y-12:y+13,x-12:x+13]
        d=np.array([q[rng] for q in data]);w=np.array([q[rng] for q in weights])
        leave=[]
        for k in range(len(data)):
            wt=w.copy();wt[k]=0
            leave.append(float(score_tile(d,wt,bank)[:,12,12].max()))
        flags=[]
        if dist<=15:flags.append('DES stationary-source mask')
        if min(leave)<3:flags.append('fails single-exposure removals')
        out.append([row['ra_deg'],row['dec_deg'],row['score'],f'{dist*.1:.4f}',f'{min(leave):.4f}',names[np.argmin(leave)],' | '.join(flags) or 'requires further vetting'])
    with (HERE/'dense_vetting.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['ra_deg','dec_deg','score','distance_to_DES_arcsec','min_loo','removed_exposure','classification'])
        w.writerows(out)
    print('DES uncontaminated',sum(float(row[3])>1.5 for row in out))
    print('survivors',[(row[0],row[1],row[2]) for row in out if row[-1]=='requires further vetting'])

if __name__=='__main__': main()
