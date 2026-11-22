"""Summarize empirical injection recovery and control selection denominators."""
from collections import defaultdict
import csv
import gzip
import json
from math import sqrt
from pathlib import Path


def wilson(k,n):
    if not n:
        return None,None
    z=1.95996398454
    p=k/n
    den=1+z*z/n
    mid=(p+z*z/(2*n))/den
    half=z*sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return max(0,mid-half),min(1,mid+half)


def write_table(path, rows, fields):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def summarize(out):
    table=defaultdict(lambda:[0,0])
    with gzip.open(out/'injections.jsonl.gz','rt') as f:
        for line in f:
            r=json.loads(line)
            if r['status']!='injected' or r['was_blended']:
                continue
            selections=[('peak',r['recovered'])]
            if r['family']=='lya':
                selections.append(('EW_lower_gt_240',r['selected_ew']))
            else:
                selections.append(('HeII_Lya_pair_CIV_veto',r.get('selected_heii_pair',False)))
            center=r['wave']
            band=f'{int(center//250)*250}-{int(center//250)*250+250}'
            flux=r['input_flux']
            flux_bin=next(label for upper,label in [
                (1e-16,'<1e-16'),(3e-16,'1e-16..3e-16'),(1e-15,'3e-16..1e-15'),
                (float('inf'),'>=1e-15')] if flux<upper)
            groups=[('all','all'),('nominal_snr',r['target_snr']),('sigma_A',r['sigma_A']),
                    ('wave_air_A',band),('input_flux_cgs',flux_bin)]
            if r['family']=='lya' and r.get('input_ew0') is not None:
                groups.append(('input_EW0_A',r['input_ew0']))
            for sel,yes in selections:
                for dim,level in groups:
                    pair=table[(r['family'],sel,dim,str(level))]
                    pair[0]+=int(bool(yes))
                    pair[1]+=1
    rows=[]
    for (family,sel,dim,level),(k,n) in sorted(table.items()):
        low,hi=wilson(k,n)
        rows.append(dict(family=family,selection=sel,dimension=dim,level=level,
                         recovered=k,trials=n,fraction=k/n,wilson95_low=low,wilson95_high=hi))
    write_table(out/'completeness.csv',rows,['family','selection','dimension','level','recovered',
                                            'trials','fraction','wilson95_low','wilson95_high'])
    return rows
