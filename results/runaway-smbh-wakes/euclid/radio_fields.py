"""Search VizieR radio source positions over complete fields, then cone-match all ridges.

Field extraction amortizes requests and preserves exact lookup coverage even
when individual target cones fail. Zero catalogue matches are compatible with
incomplete survey coverage or a radio-faint source.
"""
import argparse
import csv
import json
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np

from acquire import sha
RADIO = {'VLASS_QL1':'J/ApJS/255/30/comp', 'LoTSS_DR2':'J/A+A/659/A1/catalog',
         'RACS_low_DR1':'J/other/PASA/38.58/gausscut'}


def unit(ra,dec):
    a=np.deg2rad(ra);d=np.deg2rad(dec)
    return np.column_stack([np.cos(d)*np.cos(a),np.cos(d)*np.sin(a),np.sin(d)])


def load(table,targets,cache):
    points=unit([float(r['ra_deg']) for r in targets],[float(r['dec_deg']) for r in targets])
    mean=points.mean(axis=0);mean/=np.linalg.norm(mean)
    ra=np.rad2deg(np.arctan2(mean[1],mean[0]))%360;dec=np.rad2deg(np.arcsin(mean[2]))
    radius=np.rad2deg(np.max(np.arccos(np.clip(points@mean,-1,1))))+.02
    # Three distinct field clusters keep these queries compact.
    params={'-source':table,'-c':f'{ra:.8f} {dec:.8f}','-c.rd':f'{radius:.7f}',
            '-out.max':130000,'-out':'_RAJ2000,_DEJ2000'}
    base='https://vizier.cfa.harvard.edu/viz-bin/asu-tsv'
    url=base+'?'+urlencode(params);dest=cache/(sha(url.encode())+'.tsv')
    if not dest.exists():
        with urlopen(url,timeout=90) as f: raw=f.read()
        if b'#RESOURCE=' not in raw or ('#Name: '+table).encode() not in raw:raise ValueError('unexpected VizieR TSV')
        dest.write_bytes(raw)
    raw=dest.read_bytes();lines=raw.decode().splitlines()
    if any(re.match(r'#(INFO|ERROR).*?(truncat|overflow|limit reached|error=)',l,re.I) for l in lines):
        raise ValueError('VizieR error/row limit')
    coords=[]
    for l in lines:
        if l.startswith('#') or not l.strip():continue
        row=l.split('\t')
        if len(row)==2:
            try:coords.append([float(i) for i in row])
            except ValueError:continue
    if len(coords)>=130000:raise ValueError('VizieR row truncation')
    coords=np.array(coords).reshape((-1,2));sky=unit(*coords.T)
    result=[]
    for r,v in zip(targets,points):
        lim=20+float(r['length_arcsec'])/2
        da=np.rad2deg(np.arctan2(np.linalg.norm(np.cross(sky,v),axis=1),sky@v))*3600
        hit=np.flatnonzero(da<=lim)
        result.append(dict(id=r['id'],table=table,radius_arcsec=lim,matches=[
            dict(ra_deg=coords[i,0],dec_deg=coords[i,1],sep_arcsec=da[i]) for i in hit]))
    rec=dict(table=table,url=url,sha256=sha(raw),bytes=len(raw),rows=len(coords),radius_deg=radius,targets=len(targets))
    return rec,result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['candidates','metrics','cache','out']:p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args();a.cache.mkdir(exist_ok=True,parents=True)
    r={v['id']:v for v in csv.DictReader(a.candidates.open())}
    targets=[r[m['id']] for m in csv.DictReader(a.metrics.open()) if m['band']=='VIS' and m['sign']=='positive' and m['usable']=='1' and
      int(m['positive_bins'])>=4 and float(m['central_ridge_contrast_sigma'])>=3.]
    groups={}
    for v in targets:
        d=float(v['dec_deg']);field='N' if d>0 else 'S' if d< -40 else 'F'
        names={'N':['LoTSS_DR2','VLASS_QL1'],'S':['RACS_low_DR1'],'F':['RACS_low_DR1','VLASS_QL1']}[field]
        for s in names:groups.setdefault((field,s),[]).append(v)
    def go(k):
        field,survey=k
        try:
            rec,r=load(RADIO[survey],groups[k],a.cache);return dict(field=field,survey=survey,**rec),r
        except (OSError,ValueError) as e:
            return dict(field=field,survey=survey,table=RADIO[survey],error=repr(e)),[]
    with ThreadPoolExecutor(3) as p:all_groups=list(p.map(go,groups))
    (a.out/'radio_field_queries.json').write_text(json.dumps([v for v,_ in all_groups],indent=1)+'\n')
    (a.out/'radio_field_matches.json').write_text(json.dumps([r for _,records in all_groups for r in records],indent=1)+'\n')


if __name__=='__main__':main()
