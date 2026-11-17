"""Query NED redshifts and SIMBAD bibliography for catalogue-linked screened positions."""
import argparse
import csv
import hashlib
import io
import json
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent.parent
NED='https://ned.ipac.caltech.edu/tap/sync'
SIMBAD='https://simbad.cfa.harvard.edu/simbad/sim-tap/sync'

def archive(path,name,base,q,cols):
    limit=50000
    url=base+'?'+urlencode(dict(REQUEST='doQuery',LANG='ADQL',QUERY=q,FORMAT='csv',MAXREC=limit))
    with urlopen(url,timeout=160) as r: raw=r.read(); status=r.status
    rows=pd.read_csv(io.BytesIO(raw))
    if not set(cols).issubset(rows.columns) or len(rows)>=limit:
        raise ValueError(f'Malformed or capped {name} response')
    (path/(name+'.csv')).write_bytes(raw)
    receipt=dict(url=url,query=q,status=status,utc=datetime.now(timezone.utc).isoformat(),rows=len(rows),
                 sha256=hashlib.sha256(raw).hexdigest(),bytes=len(raw))
    (path/(name+'_query.json')).write_text(json.dumps(receipt,indent=2)+'\n')
    return rows

def literature(field,ranks):
    path=field/'catalogs'
    m=pd.read_csv(path/'matches.csv')
    r=pd.read_csv(ranks/'ranks.csv')
    shortlist=r[(r.score>=5)|(r.peer_rank<=3)|(r.median_abs_residual.rank(ascending=False)<=3)]
    oids=list(m.loc[m.source_id.isin(shortlist.source_id),'oid'].unique())
    if not oids:return
    ids=','.join(str(v) for v in sorted(oids))
    q=f"""SELECT h.oidref,r.bibcode,r.title,r.doi FROM has_ref AS h
JOIN ref AS r ON h.oidbibref=r.oidbib WHERE h.oidref IN ({ids})
ORDER BY oidref,bibcode"""
    archive(path,'references',SIMBAD,q,['oidref','bibcode','title'])

def ned(field,ra,dec):
    path=field/'catalogs'
    q=f"""SELECT objid,prefname,ra,dec,prefphytype,z,zflag,z_bibcode FROM NEDTAP.objdir
WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{dec},0.52))
AND z IS NOT NULL ORDER BY objid"""
    cat=archive(path,'ned_redshifts',NED,q,['ra','dec','z'])
    ra0,dec0=np.deg2rad(cat.ra),np.deg2rad(cat.dec)
    rows=[]
    for s in pd.read_csv(field/'sources.csv').itertuples():
        a,d=np.deg2rad([s.ra,s.dec])
        sep=np.rad2deg(2*np.arcsin(np.sqrt(np.clip(np.sin((dec0-d)/2)**2+
                     np.cos(d)*np.cos(dec0)*np.sin((ra0-a)/2)**2,0,1))))*3600
        for idx in np.flatnonzero(sep<=3):
            rows.append(dict(source_id=s.source_id,offset_arcsec=float(sep[idx]),**cat.iloc[idx].to_dict()))
    pd.DataFrame(rows,columns=['source_id','offset_arcsec',*cat.columns]).to_csv(path/'ned_matches.csv',index=False)
if __name__=='__main__':
    jobs=[(literature,(HERE,HERE/'nep_two_band')),(literature,(HERE/'wide',HERE/'wide/full')),
          (ned,(HERE,270,66.56)),(ned,(HERE/'wide',209.53125,30))]
    with ThreadPoolExecutor(4) as ex:
        futures=[ex.submit(fn,*args) for fn,args in jobs]
        for fut in futures:fut.result()
