"""Cross-check the Gaia mass screen against spectral identifications of double-lined systems."""
import argparse
import csv
import hashlib
import json
import math
import sys
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'epochs'))
from catalog_io import read_header, table_dtype

LAMOST=('https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?'+urllib.parse.urlencode(
        {'-source':'J/MNRAS/527/521/tableb1','-out.all':'1','-out.max':'15000'}))
APOGEE='https://data.sdss.org/sas/dr17/apogee/vac/apogee-sb2/apogee_sb2s-v1_0.fits'


def run(cache):
    cache.mkdir(parents=True,exist_ok=True)
    inputs={}
    for name,url in [('vizier.cfa.harvard.edu.tsv',LAMOST),('apogee_sb2.fits',APOGEE)]:
        dst=cache/name
        if not dst.exists():dst.write_bytes(urllib.request.urlopen(url,timeout=120).read())
        inputs[name]=dict(url=url,sha256=hashlib.sha256(dst.read_bytes()).hexdigest())
    with (HERE.parent/'gaia/source_diagnostics.csv').open() as f: stars=list(csv.DictReader(f))
    index={}
    for line in (cache/'vizier.cfa.harvard.edu.tsv').read_text().splitlines():
        row=line.split('\t')
        if len(row)==5 and row[0].strip().isdecimal():
            index[row[2].strip()]=row[1]
    if len(index)!=12426:raise ValueError('Incomplete Kovalev et al. table')
    with (cache/'apogee_sb2.fits').open('rb') as f:
        prim=read_header(f)
        n=int(prim['NAXIS1'])
        f.seek((n+2879)//2880*2880,1)
        h=read_header(f)
        d=np.fromfile(f,dtype=table_dtype(h),count=int(h['NAXIS2']))
    if len(d)!=31313:raise ValueError('Incomplete DR17 SB2 VAC')
    ra,dec=np.radians(d['RA']),np.radians(d['DEC'])
    hits=[]
    for star in stars:
        sid=star['source_id']
        if sid in index:hits.append(dict(source_id=sid,survey='LAMOST_MRS_Kovalev',name=index[sid],separation_arcsec=''))
        # The VAC uses 2MASS-like identifiers. Propagate DR3 to approximate epoch 2000.
        delta=float(star['dec'])
        dt=2000-float(star['ref_epoch'])
        x,y=np.radians([float(star['ra'])+dt*float(star['pmra'] or 0)/3600000/math.cos(math.radians(delta)),
                        delta+dt*float(star['pmdec'] or 0)/3600000])
        chord=np.sin((y-dec)/2)**2+np.cos(y)*np.cos(dec)*np.sin((x-ra)/2)**2
        sel=chord < np.sin(np.radians(2/3600)/2)**2
        for name in np.unique(d['OBJID'][sel]):
            ix=sel & (d['OBJID']==name)
            separation=np.degrees(2*np.arcsin(np.sqrt(chord[ix].min())))*3600
            hits.append(dict(source_id=sid,survey='APOGEE_Kounkel_DR17',name=name.decode().strip(),
                             separation_arcsec=round(float(separation),4)))
    with (HERE/'sb2_hits.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,['source_id','survey','name','separation_arcsec'],lineterminator='\n')
        w.writeheader();w.writerows(hits)
    return dict(inputs=inputs,screen_stars=len(stars),lamost_sb2_ids=len(index),
                apogee_sb2_ids=len(np.unique(d['OBJID'])),apogee_vac_rows=len(d),n_hits=len(hits))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cache',required=True,type=Path)
    (HERE/'sb2_depth.json').write_text(json.dumps(run(p.parse_args().cache),indent=2)+'\n')
