"""Check three degeneracy-limited PRIME peaks against earlier VVV and local KMT stars."""
import csv
import json
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode

import numpy as np

from crossmatch import coords
from fetch import OUT,get,write_csv,SHA
from independent import retrieve,KMT

CACHE=Path('.survey-cache/prime_followups')
NAMES={'PRIME-2024-BLG-245','PRIME-2025-BLG-0239','PRIME-2025-BLG-0360'}


def inspect(row):
    n=row['event'];ra,dec=coords(row['ra'],row['dec']);out=[];receipts=[]
    query=KMT+urlencode(dict(search_type='radec',ra=ra,dec=dec,ra_type='decimal_degree',arcsec=2))
    src=retrieve(CACHE/(n+'.html'),query)
    hits=json.loads(re.search(r'const kmtnetSearchMarkers = (.*?);',src.decode())[1])
    receipts.append(dict(event=n,url=query,sha256=SHA(src),bytes=len(src),product='field lookup',rows=len(hits)))
    seeds=json.loads((Path('staging/flag-refits')/(n+'.json')).read_text())
    t0,te,rho,u=[seeds[k] for k in ('fs_t0','fs_te_or_rise','rho','fs_u0_or_fall')]
    import VBMicrolensing as vbm
    eng=vbm.VBMicrolensing();eng.LoadESPLTable(str(Path(vbm.__file__).parent/'data/ESPL.tbl'))
    for hit in hits:
        name=hit['lcid']
        url=KMT+urlencode(dict(action='download_lc',lcid=name));raw=retrieve(CACHE/(name+'.csv'),url)
        rr=list(csv.DictReader(raw.decode().splitlines()))
        receipts.append(dict(event=n,url=url,sha256=SHA(raw),bytes=len(raw),product='KMT field photometry',rows=len(rr)))
        year=n.split('-')[1]
        for site in ('BLG','SAO','SSO'):
            subset=[r for r in rr if r['year']==year and r['site']==site]
            if not subset:continue
            data=np.array([[float(r[k]) for k in ('date','flux','err','chi2')] for r in subset])
            t,f,e,ch=data.T;t=t+2400000-t0
            ok=(e>0)&(ch>0)&(ch<20)&np.isfinite(t+f+e)
            t,f,e=t[ok],f[ok],e[ok];off=abs(t)>max(10,2*te*max(rho,1))
            if len(t)<5 or off.sum()<4:continue
            e=np.maximum(e,.15*np.median(e))
            z=(f[off]-np.median(f[off]))/e[off];scale=max(1.,np.median(abs(z))/.67448975)
            w=1/(e*scale)**2
            fs=np.array(eng.ESPLLightCurve([np.log(u),np.log(te),0.,np.log(rho)],t.tolist())[0])-1
            fs-=np.average(fs,weights=w)
            y=f-np.average(f,weights=w)
            flux=max(0,np.dot(w*fs,y)/np.dot(w,fs*fs))
            gain=np.dot(w,y*y)-np.dot(w,(y-flux*fs)**2)
            out.append(dict(event=n,id=name,separation_arcsec=hit['dist'],site=site,source_year=year,usable=len(t),
                            n_inside_source_radius_crossing=int(sum(abs(t)<rho*te)),before_closest_d=round(-max(t[t<0]),4) if (t<0).any() else '',
                            after_closest_d=round(min(t[t>0]),4) if (t>0).any() else '',fs_shape_gain=round(float(gain),2),
                            positive_template_flux=round(float(flux),3),empirical_error_factor=round(scale,4)))
    vurl='https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?'+urlencode({'-source':'II/387/virac2','-c':f'{ra:.9f} {dec:.9f}','-c.rs':2,
       '-out.max':12,'-out.add':'_r','-sort':'_r','-out':'srcid,RAJ2000,DEJ2000,Jmag,Hmag,Ksmag,o_Ksmag,MedKsmag,MADKs,etaKs,KsEp1,KsEpl'})
    v=retrieve(CACHE/(n+'.tsv'),vurl)
    lines=[line for line in v.decode().splitlines() if line and line[0]!='#']
    vrows=list(csv.DictReader([*lines[:1],*lines[3:]],delimiter='\t')) if lines else []
    receipts.append(dict(event=n,url=vurl,sha256=SHA(v),bytes=len(v),product='VIRAC2',rows=len(vrows),
                         status='error' if b'#ERROR' in v or not lines else 'ok'))
    return receipts,out,[dict(event=n,**r) for r in vrows]


if __name__=='__main__':
    CACHE.mkdir(exist_ok=True)
    rows=[r for r in csv.DictReader((OUT/'index.csv').open()) if r['event'] in NAMES]
    with ThreadPoolExecutor(3) as pool:records=list(pool.map(inspect,rows))
    for j,name in enumerate(('prime_candidate_receipts','prime_candidate_kmt','prime_candidate_virac2')):
        write_csv(OUT/(name+'.csv'),[r for result in records for r in result[j]])
