"""Propagate catalogue-redshift intervals to geometric associations of ridges."""
import argparse
import ast
import csv
import json
from pathlib import Path
from astropy.cosmology import Planck18

from acquire import write_csv
from host_audit import ESA,tap


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('candidates','metrics','cache','out'):p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args()
    byid={r['id']:r for r in csv.DictReader(a.candidates.open())}
    targets=[byid[r['id']] for r in csv.DictReader(a.metrics.open()) if r['band']=='VIS' and r['sign']=='positive' and
        r['usable']=='1' and int(r['positive_bins'])>=4 and float(r['central_ridge_contrast_sigma'])>=3.]
    ids=sorted({s for r in targets for s in [r['host_object_id'],r['axis_host_id'],*r['interior_mer_ids'].split()] if s},key=int)
    phz={};receipts=[]
    for i in range(0,len(ids),50):
        batch=ids[i:i+50]
        q="SELECT object_id,phz_median,phz_95_int,phz_classification,phz_flags FROM q1.phz_photo_z WHERE object_id IN ("+','.join(batch)+") ORDER BY object_id"
        rows,receipt=tap(ESA,q,a.cache);receipts.append(dict(targets=batch,**receipt))
        phz.update({r['object_id']:r for r in rows})
    if len(phz)>len(ids) or set(phz)-set(ids):raise ValueError('unexpected PHZ IDs')
    write_csv(a.out/'component_phz.csv',[phz[k] for k in ids if k in phz])
    (a.out/'phz_queries.json').write_text(json.dumps(receipts,indent=1)+'\n')
    profiles=[]
    for r in targets:
        for host in dict.fromkeys(k for k in (r['host_object_id'],r['axis_host_id']) if k):
            interior={k for k in r['interior_mer_ids'].split() if k!=host}
            objs=[phz[k] for k in [host,*sorted(interior)] if k in phz]
            intervals=[]
            for o in objs:
                if not o['phz_95_int']:continue
                lo,hi=ast.literal_eval(o['phz_95_int'])
                if 0<=lo<=hi and hi>0:intervals.append((lo,hi))
            equal_z=''
            if len(intervals)>=2 and len(intervals)==len(interior)+1:
                equal_z=int(max(q[0] for q in intervals)<=min(q[1] for q in intervals))
            z=phz.get(host,{}).get('phz_median','')
            kpc=''
            if z and float(z)>0:
                kpc=float(Planck18.kpc_proper_per_arcmin(float(z)).value*float(r['length_arcsec'])/60)
            profiles.append(dict(id=r['id'],host_id=host,endpoint_match=int(host==r['host_object_id']),
                axial_match=int(host==r['axis_host_id']),interior_other_ids=' '.join(sorted(interior)),
                host_phz=z,number_valid_intervals=len(intervals),shared_95_interval_nonempty=equal_z,
                conditional_projected_length_kpc=kpc))
    write_csv(a.out/'host_redshift_geometry.csv',profiles)


if __name__=='__main__':main()
