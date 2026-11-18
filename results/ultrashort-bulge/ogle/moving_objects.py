"""Known-body SkyBoT cones at three observed bright epochs per OGLE short alert."""
import argparse
import concurrent.futures as cf
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

import erfa
import numpy as np
from counterparts import coords, HERE


def times_for(cache,row):
    y,_,n=row["event"].split("-")
    t,m,e=np.loadtxt(cache/y/f"blg-{n}.dat",usecols=(0,1,2)).T
    t0=float(row["t0_hjd"])
    floor=np.hypot(e,.004)
    good=np.isfinite(m)&(e>0)&(e<1)&(m>5)&(m<25)
    baseline=np.median(m[good&(abs(t-t0)>30)])
    excess=(baseline-m)/floor
    ind=np.flatnonzero(good&(abs(t-t0)<1.5)&(excess>4))
    selected=[ind[0],ind[np.argmax(excess[ind])],ind[-1]]
    return t[np.unique(selected)]


def query(job):
    row,hjd=job
    ra,dec=coords(row["ra_hms"],row["dec_dms"])
    a,d=np.radians([ra,dec]); direction=np.array([np.cos(a)*np.cos(d),np.sin(a)*np.cos(d),np.sin(d)])
    pvh,_=erfa.epv00(hjd,0.)
    # HJD_UTC ~= JD_UTC + r_helio(Earth).n/c. Seconds of timing uncertainty are
    # negligible relative to this generous 120-arcsecond cone.
    jd=hjd-(pvh['p']@direction)*499.004783836/86400
    args={'-ep':f'{jd:.7f}','-ra':f'{ra:.8f}','-dec':f'{dec:.8f}','-rs':120,'-observer':'304',
          '-filter':0,'-output':'basic','-from':'OpenSpaceResearch','-mime':'votable'}
    url='https://ssp.imcce.fr/webservices/skybot/api/conesearch.php?'+urllib.parse.urlencode(args)
    out=dict(event=row['event'],hjd_observed=float(hjd),jd_utc_query=float(jd),url=url,
             queried_utc=datetime.now(timezone.utc).isoformat(timespec='seconds'))
    try:
        with urllib.request.urlopen(url,timeout=60) as resp: data=resp.read()
        doc=ET.fromstring(data)
        infos={e.get('name'):e.get('value') for e in doc.findall('.//{*}INFO')}
        errors=[e.text or '' for e in doc.findall('.//{*}INFO') if e.get('name')=='QUERY_STATUS' and e.get('value')=='ERROR']
        empty=bool(errors) and all(m.startswith('No solar system object was found in the requested FOV') for m in errors)
        if infos.get('QUERY_STATUS')!='OK' and not empty: raise ValueError(str(infos))
        if empty: out['service_detail']='explicit empty field' 
        out.update(sha256=hashlib.sha256(data).hexdigest(),result_n=len(doc.findall('.//{*}TR')),
                   columns=[e.get('name') for e in doc.findall('.//{*}FIELD')],
                   matches=[[e.text for e in tr.findall('{*}TD')] for tr in doc.findall('.//{*}TR')],
                   status='ok')
    except Exception as ex:
        out['status']=str(ex)
    return out


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('cache',type=Path);cache=ap.parse_args().cache
    rows=list(csv.DictReader((HERE/'catalog.csv').open()))
    jobs=[(r,hjd) for r in rows if 0<float(r['tau_days'])<.5 for hjd in times_for(cache,r)]
    with cf.ThreadPoolExecutor(3) as pool: records=list(pool.map(query,jobs))
    (HERE/'moving_objects.json').write_text(json.dumps(records,indent=2)+'\n')

if __name__=='__main__':
    main()
