"""Fetch all candidate-neighborhood matches. Raw response hash is before gzip."""
import argparse
import concurrent.futures
import csv
import datetime
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import time
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from urllib.error import HTTPError

ROOT = Path(__file__).resolve().parent
TAP = 'https://irsa.ipac.caltech.edu/TAP/sync'
TABLES = {
    'neo': ('neowiser_p1bs_psd', 'source_id,scan_id,frame_num,ra,dec,mjd,qual_frame,qi_fact,saa_sep,cc_flags,moon_masked,ph_qual,nb,w1sat,w2sat,w1mpro,w1sigmpro,w1snr,w1rchi2,w2mpro,w2sigmpro,w2snr,w2rchi2', 2024.25),
    'catwise': ('catwise_2020', 'source_name,source_id,unwise_objid,ra,dec,meanobsmjd,cc_flags,ab_flags,w1mpro,w1sigmpro,w2mpro,w2sigmpro,w1mpro_pm,w1sigmpro_pm,w2mpro_pm,w2sigmpro_pm,w1rchi2,w2rchi2', 2015.4),
    'unwise': ('unwise_2019', 'unwise_objid,ra,dec,flux_1,dflux_1,flux_2,dflux_2,fracflux_1,fracflux_2,flags_unwise_1,flags_unwise_2,primary_1,primary_2,qf_1,qf_2', 2014.5),
    'allwise': ('allwise_p3as_psd', 'designation,ra,dec,cc_flags,ext_flg,nb,ph_qual,var_flg,w1mpro,w1sigmpro,w1snr,w2mpro,w2sigmpro,w2snr,w3mpro,w3sigmpro,w3snr,w4mpro,w4sigmpro,w4snr', 2010.5)
}

def position(target, epoch):
    ra, dec, delta = float(target['ra']), float(target['dec']), epoch-float(target['ref_epoch'])
    return ra+delta*float(target['pmra'])/(3.6e6*math.cos(math.radians(dec))), dec+delta*float(target['pmdec'])/3.6e6

def separation(ra, dec, ra0, dec0):
    return 3600*math.hypot(((ra-ra0+180)%360-180)*math.cos(math.radians(dec0)), dec-dec0)

def targets():
    with (ROOT.parent/'literature'/'candidates.csv').open() as f:
        return list(csv.DictReader(f))

def request(kind,target):
    table, columns, epoch = TABLES[kind]
    ra,dec=position(target,epoch)
    q = f"SELECT {columns} FROM {table} WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra:.8f},{dec:.8f},0.0033333333))=1"
    if kind=='neo':
        q+=' AND mjd>=60291 AND mjd<60524'
    url=TAP+'?'+urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',MAXREC=20000,QUERY=q))
    return q,url

def fetch(kind,target,dest):
    dest = Path(dest); dest.mkdir(parents=True,exist_ok=True)
    key = target['label']+'_'+kind
    out, receipt = dest/(key+'.csv.gz'), dest/(key+'.json')
    q,url=request(kind,target)
    columns=TABLES[kind][1]
    error = dest/(key+'.error.json')
    error_body = dest/(key+'.error.gz')
    if error.exists() or error_body.exists():
        raise ValueError(f'Previous error receipt exists for {key}; use a new campaign')
    if out.exists() != receipt.exists():
        raise ValueError(f'Incomplete frozen pair for {key}; use a new campaign')
    if out.exists():
        try:
            r=json.loads(receipt.read_text(encoding='utf-8'))
            data=gzip.decompress(out.read_bytes())
            valid=(r['url']==url and r['query']==q and r['bytes']==len(data)
                   and r['sha256']==hashlib.sha256(data).hexdigest())
        except (OSError,ValueError,KeyError,EOFError) as exc:
            raise ValueError(f'Invalid frozen pair for {key}') from exc
        if not valid: raise ValueError(f'Frozen pair mismatch for {key}; use a new campaign')
        return
    payload=None; status=None
    try:
        with urlopen(Request(url,headers={'User-Agent':'Open-Space research, SPHEREx candidate follow-up'}), timeout=120) as resp:
            status=resp.status; payload=resp.read()
        reader=csv.DictReader(io.StringIO(payload.decode('utf-8-sig')))
        if reader.fieldnames!=columns.split(','): raise ValueError('Unexpected response columns')
        n=sum(1 for _ in reader)
        if n>=20000: raise ValueError('Reached MAXREC')
        record=dict(url=url,query=q,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),status=status,bytes=len(payload),rows=n,sha256=hashlib.sha256(payload).hexdigest())
        with out.open('xb') as f:f.write(gzip.compress(payload,mtime=0))
        with receipt.open('x',encoding='utf-8') as f:json.dump(record,f,sort_keys=True,indent=2);f.write('\n')
        print(key,n,flush=True)
    except Exception as exc:
        if isinstance(exc,HTTPError): status=exc.code;payload=exc.read()
        record=dict(url=url,query=q,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),status=status,error=repr(exc))
        if payload is not None:
            record.update(bytes=len(payload),sha256=hashlib.sha256(payload).hexdigest())
            with error_body.open('xb') as f:f.write(gzip.compress(payload,mtime=0))
        with error.open('x',encoding='utf-8') as f:json.dump(record,f,sort_keys=True,indent=2);f.write('\n')
        raise

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--workers',type=int,default=4)
    ap.add_argument('--output',type=Path,required=True,help='Campaign directory for immutable request pairs')
    args=ap.parse_args()
    with concurrent.futures.ThreadPoolExecutor(args.workers) as pool:
        jobs=[pool.submit(fetch,k,t,args.output) for t in targets() for k in TABLES]
        for job in concurrent.futures.as_completed(jobs):
            job.result()
