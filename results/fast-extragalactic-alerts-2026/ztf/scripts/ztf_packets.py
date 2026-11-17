"""Fetch original ZTF candidate metadata and science/reference/difference stamps.

An unstamped historical detection is traced to its stamped parent. This preserves
the distinction between images of the selected epoch and an alert's history.
"""
from concurrent.futures import ThreadPoolExecutor
import csv
import datetime as dt
import gzip
import hashlib
import json
from urllib.parse import urlencode
from urllib.request import urlopen
from ztf_tap import ROOT, RAW, META
from ztf_screen import write_csv

def collect(target):
    oid,cid=target
    name=oid+'_'+cid
    for kind in ['info','science','template','difference']:
        params={'oid':oid,'candid':cid}
        if kind=='info':
            endpoint='get_avro_info'
        else:
            endpoint='get_stamp'
            params.update(type=kind,format='png')
        url='https://avro.alerce.online/'+endpoint+'?'+urlencode(params)
        with urlopen(url,timeout=45) as h:
            raw=h.read()
            rec=dict(url=url,method='GET',http_status=h.status,retrieved_utc=dt.datetime.now(dt.timezone.utc).isoformat(),sha256=hashlib.sha256(raw).hexdigest(),bytes=len(raw),mime=h.headers.get('Content-Type'))
        if kind=='info':
            if str(json.loads(raw)['candid'])!=cid: raise ValueError('Wrong candidate returned')
            dest=RAW/'packets'/(name+'.json.gz')
            dest.write_bytes(gzip.compress(raw,mtime=0))
        else:
            if not raw.startswith(b'\x89PNG'): raise ValueError('Non-PNG returned')
            dest=RAW/'packets'/(name+'_'+kind+'.png')
            dest.write_bytes(raw)
        rec['file']=str(dest.relative_to(ROOT))
        (META/'packets'/(name+'_'+kind+'.json')).write_text(json.dumps(rec,indent=2)+'\n')
    print(oid,cid,flush=True)

def main():
    for path in (RAW/'packets',META/'packets'): path.mkdir(exist_ok=True)
    with (ROOT/'derived'/'ztf_pairs.csv').open() as f:
        pairs={r['oid']:r for r in csv.DictReader(f)}
    epochs={}
    with gzip.open(RAW/'detection_followup.csv.gz','rt') as f:
        for row in csv.DictReader(f):
            oid,cid=row['oid'],row['candid']
            if cid in (pairs[oid]['candid1'],pairs[oid]['candid2']):
                source=cid if row['has_stamp']=='True' else row['parent_candid']
                epochs[(oid,cid)]=source
    with ThreadPoolExecutor(max_workers=4) as ex: list(ex.map(collect,sorted({(oid,src) for (oid,cid),src in epochs.items()})))
    audit=[]
    for (oid,cid),src in sorted(epochs.items()):
        with gzip.open(RAW/'packets'/f'{oid}_{src}.json.gz','rt') as f: packet=json.load(f)
        epoch=next((e for e in [packet['candidate'],*(packet['prv_candidates'] or [])] if str(e.get('candid'))==cid),None)
        if epoch is None: raise ValueError(f'{oid}/{cid} missing from {src}')
        audit.append(dict(oid=oid,candid=cid,packet_candid=src,at_image_epoch=cid==src,programid=epoch.get('programid'),fid=epoch.get('fid'),jd=epoch.get('jd'),positive=epoch.get('isdiffpos')))
    write_csv(ROOT/'derived'/'ztf_public_program.csv',audit,audit[0].keys())
if __name__=='__main__':
    main()
