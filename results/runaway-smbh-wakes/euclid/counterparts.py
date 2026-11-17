"""Catalogue checks of Euclid ridges passing the local longitudinal contrast test.

Search DESI DR1 at finite apertures and return measured redshifts for any
spectra found. Retrieve MER photometric intervals for distinct counterparts.
A blank radio catalogue cone is a coverage-dependent positional check.
"""
import argparse
import csv
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
from xml.etree import ElementTree as ET

import numpy as np

from acquire import sha
from host_audit import ESA,tap

DESI='https://datalab.noirlab.edu/tap/sync'
VIZ='https://vizier.cfa.harvard.edu/viz-bin/asu-xml'
RADIO={'VLASS_QL1':'J/ApJS/255/30/comp','LoTSS_DR2':'J/A+A/659/A1/catalog','RACS_low_DR1':'J/other/PASA/38.58/gausscut'}


def radio(task,cache):
    r,survey=task
    radius=20+float(r['length_arcsec'])/2
    p={'-source':RADIO[survey],'-c':f"{float(r['ra_deg']):.8f} {float(r['dec_deg']):.8f}",
       '-c.rs':f'{radius:.4f}','-out.max':100,'-out':'_r,_RAJ2000,_DEJ2000'}
    u=VIZ+'?'+urlencode(p)
    d=cache/(sha(u.encode())+'.xml')
    try:
        if not d.exists():
            with urlopen(u,timeout=30) as f:b=f.read()
            t=ET.fromstring(b)
            if t.find('.//{*}TABLE') is None:raise ValueError('missing VizieR table')
            d.write_bytes(b)
        raw=d.read_bytes();tree=ET.fromstring(raw)
        fs=[x.get('name') for x in tree.findall('.//{*}FIELD')]
        if fs!=['_r','_RAJ2000','_DEJ2000']:raise ValueError('unexpected catalogue fields')
        rows=[dict(zip(fs,[t.text for t in tr.findall('{*}TD')])) for tr in tree.findall('.//{*}TR')]
        return dict(id=r['id'],survey=survey,url=u,radius_arcsec=radius,sha256=sha(raw),bytes=len(raw),matches=rows)
    except (OSError,ValueError,ET.ParseError) as e:
        return dict(id=r['id'],survey=survey,url=u,error=repr(e))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['candidates','metrics','cache','out']:p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--jobs',type=int,default=5);a=p.parse_args()
    a.cache.mkdir(exist_ok=True,parents=True);a.out.mkdir(exist_ok=True,parents=True)
    f=[r for r in csv.DictReader(a.metrics.open()) if r['band']=='VIS' and r['sign']=='positive' and
       r['usable']=='1' and int(r['positive_bins'])>=4 and float(r['central_ridge_contrast_sigma'])>=3.]
    byid={r['id']:r for r in csv.DictReader(a.candidates.open())};targets=[byid[r['id']] for r in f]
    objids=sorted({s for r in targets for s in [r['host_object_id'],r['axis_host_id'],*r.get('interior_mer_ids','').split()] if s},key=int)
    queries=[]
    for i in range(0,len(objids),50):
        ids=objids[i:i+50]
        q="SELECT object_id,phz_median,phz_95_int,phz_classification,phz_flags FROM q1.phz_photo_z WHERE object_id IN ("+','.join(ids)+") ORDER BY object_id"
        queries.append(('Q1_PHZ',ids,ESA,q))
    # One small sky box per target; retain positions so on-axis vs field
    # spectra can be distinguished. Do not infer common z from one component.
    for i in range(0,len(targets),6):
        group=targets[i:i+6];regs=[]
        for r in group:
            x,y=float(r['ra_deg']),float(r['dec_deg'])
            rad=(float(r['length_arcsec'])/2+5)/3600;d=rad/np.cos(np.deg2rad(y))
            regs.append(f'(mean_fiber_ra BETWEEN {x-d:.8f} AND {x+d:.8f} AND mean_fiber_dec BETWEEN {y-rad:.8f} AND {y+rad:.8f})')
        q='SELECT targetid,mean_fiber_ra,mean_fiber_dec,z,zerr,zwarn,spectype,deltachi2,zcat_primary FROM desi_dr1.zpix WHERE '+' OR '.join(regs)+' ORDER BY targetid'
        queries.append(('DESI_DR1',[r['id'] for r in group],DESI,q))
    def request(task):
        name,ids,ep,q=task
        try:
            rows,receipt=tap(ep,q,a.cache)
            return dict(survey=name,targets=ids,**receipt,matches=rows)
        except (OSError,ValueError) as exc:
            return dict(survey=name,targets=ids,query=q,endpoint=ep,error=repr(exc))
    radios=[(r,s) for r in targets for s in RADIO if (s=='RACS_low_DR1' and float(r['dec_deg'])<41)
        or (s=='VLASS_QL1' and float(r['dec_deg'])>=-40) or (s=='LoTSS_DR2' and float(r['dec_deg'])>0)]
    with ThreadPoolExecutor(a.jobs) as pool:
        rf=[pool.submit(radio,t,a.cache) for t in radios]
        zf=[pool.submit(request,t) for t in queries]
        zr=[f.result() for f in zf];rr=[f.result() for f in rf]
    (a.out/'redshift_checks.json').write_text(json.dumps(zr,indent=1)+'\n')
    (a.out/'radio_checks.json').write_text(json.dumps(rr,indent=1)+'\n')


if __name__=='__main__':main()
