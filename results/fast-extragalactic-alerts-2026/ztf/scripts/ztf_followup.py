"""Collect broker catalog checks, raw photometry/limits and candidate image triplets."""
from concurrent.futures import ThreadPoolExecutor
import csv
import datetime as dt
import gzip
import hashlib
import json
import os
from pathlib import Path
import urllib.request
from ztf_tap import ROOT, RAW, META, get

def fink(oid):
    if not os.environ.get('FAST_ALERTS_NEW_CAMPAIGN') or Path(os.environ['FAST_ALERTS_NEW_CAMPAIGN']).resolve() != ROOT.resolve():
        raise RuntimeError('Use a separate new campaign for live retrieval')
    if ROOT.is_relative_to(Path(__file__).resolve().parents[2]):
        raise ValueError('Live responses must be outside the frozen result')
    url='https://api.ztf.fink-portal.org/api/v1/objects'
    for tag,params in [('fink_history',{'objectId':oid,'withupperlim':True,'output-format':'json'}),('fink_stamps',{'objectId':oid,'withcutouts':True,'output-format':'json','columns':'i:objectId,i:candid,i:jd,i:fid,i:magpsf,i:ra,i:dec'})]:
        if (RAW / f'{oid}_{tag}.json.gz').exists() or (META / f'{oid}_{tag}.json').exists():
            raise FileExistsError(oid + '_' + tag)
        req=urllib.request.Request(url,data=json.dumps(params).encode(),headers={'Content-Type':'application/json'})
        with urllib.request.urlopen(req,timeout=65) as h:
            raw=h.read()
            info=dict(method='POST',url=url,params=params,retrieved_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
                      http_status=h.status,bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest())
        path=RAW/f'{oid}_{tag}.json.gz'
        info['file']=str(path.relative_to(ROOT))
        path.write_bytes(gzip.compress(raw,mtime=0))
        (META/f'{oid}_{tag}.json').write_text(json.dumps(info,indent=2)+'\n')
        print(oid,tag,len(json.loads(raw)),len(raw),flush=True)

def main():
    with (ROOT/'derived'/'ztf_pairs.csv').open() as f: ids=sorted({x['oid'] for x in csv.DictReader(f)})
    where=' WHERE oid IN ('+','.join("'"+oid+"'" for oid in ids)+')'
    def tables(table):
        rows=get(table+'_followup','SELECT * FROM ztf.'+table+where)
        print(table,len(rows),flush=True)
    with ThreadPoolExecutor(max_workers=3) as ex:
        list(ex.map(tables,['ss_ztf','ps1_ztf','gaia_ztf','xmatch','probability','non_detection','forced_photometry','detection','dataquality']))
    with ThreadPoolExecutor(max_workers=3) as ex: list(ex.map(fink,ids))
if __name__=='__main__':
    main()
