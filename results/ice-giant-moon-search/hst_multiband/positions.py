"""Batch HST reception epochs through the official Horizons file API."""
import argparse
import csv
import gzip
import hashlib
import json
import time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import requests
from astropy.time import Time
ROOT=Path(__file__).resolve().parent
URL='https://ssd.jpl.nasa.gov/api/horizons_file.api'
SPECIAL={75051:'2023U1',85051:'S2002_N5',85052:'S2021_N1'}


def run(cache):
    cache.mkdir(exist_ok=True,parents=True)
    inventory=list(csv.DictReader((ROOT.parent/'manifests/observations.csv').open()))
    utc_stamps={}; coords=[];receipts=[]
    for planet,codes in [('Neptune',[899,*range(801,815),85051,85052]),('Uranus',[799,*range(701,728),75051])]:
        rows=sorted([r for r in inventory if r['obs_collection']=='HST' and r['planet']==planet],key=lambda r:(float(r['t_min'])+float(r['t_max']))/2)
        times=[(float(r['t_min'])+float(r['t_max']))/2+2400000.5 for r in rows]
        for code in codes:
            for start in range(0,len(rows),1000):
                subset=rows[start:start+1000];jd=times[start:start+1000]
                payload=(f"!$$SOF\nCOMMAND='{code}'\nOBJ_DATA='NO'\nMAKE_EPHEM='YES'\nEPHEM_TYPE='OBSERVER'\nCENTER='500@-48'\n"
                         "TLIST_TYPE='JD'\nQUANTITIES='1,19,20,21,24'\nCSV_FORMAT='YES'\nANG_FORMAT='DEG'\nEXTRA_PREC='YES'\nTIME_DIGITS='FRACSEC'\nTLIST="
                         +',\n'.join("'%.9f'"%t for t in jd)+'\n')
                key=hashlib.sha256(payload.encode()).hexdigest();out=cache/(key+'.json')
                for attempt in range(5):
                    try:
                        if out.exists():raw=out.read_bytes()
                        else:
                            response=requests.post(URL,data={'format':'json'},files={'input':('times.txt',payload)},timeout=70)
                            response.raise_for_status();raw=response.content
                        result=json.loads(raw)['result']; ident=next(v for v in result.splitlines() if v.startswith('Target body name:'))
                        if f'({code})' not in ident or SPECIAL.get(code,'') not in ident or '_merged' not in ident:
                            raise ValueError(f'Wrong target: {ident}')
                        table=list(csv.reader(result.split('$$SOE')[1].split('$$EOE')[0].strip().splitlines()))
                        if len(table)!=len(subset):raise ValueError((code,len(table),len(subset),result[-160:]))
                        if not out.exists():out.write_bytes(raw)
                        break
                    except (requests.RequestException,KeyError,ValueError,StopIteration,IndexError):
                        if attempt==4:raise
                        time.sleep(2**attempt)
                receipts.append(dict(body=code,first=subset[0]['obs_id'],last=subset[-1]['obs_id'],epochs=len(jd),
                                     request_sha256=key,response_sha256=hashlib.sha256(raw).hexdigest(),identity=ident.strip()))
                for row,t,val in zip(subset,jd,table,strict=True):
                    stamp=val[0].strip()
                    if stamp not in utc_stamps:utc_stamps[stamp]=Time(datetime.strptime(stamp,'%Y-%b-%d %H:%M:%S.%f').replace(tzinfo=timezone.utc)).jd
                    if abs(utc_stamps[stamp]-t)*86400>.003:raise ValueError('wrong returned epoch '+row['obs_id'])
                    coords.append(dict(rootname=row['obs_id'],body=code,ra=val[3].strip(),dec=val[4].strip(),r_au=val[5].strip(),
                                       delta_au=val[7].strip(),light_min=val[9].strip(),phase_deg=val[10].strip()))
            print(planet,code,len(coords),flush=True)
        with gzip.open(ROOT/(planet.lower()+'_positions.csv.gz'),'wt',newline='') as fp:
            w=csv.DictWriter(fp,fieldnames=list(coords[0]));w.writeheader();w.writerows(coords)
        coords=[]
        with (ROOT/'position_queries.csv').open('w',newline='') as fp:
            w=csv.DictWriter(fp,fieldnames=list(receipts[0]));w.writeheader();w.writerows(receipts)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('cache',type=Path);run(ap.parse_args().cache)
