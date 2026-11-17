"""Identify program IDs of all eligible ZTF exposure times from public packet metadata."""
import collections
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
import hashlib
import datetime as dt
from urllib.parse import urlencode
from urllib.request import urlopen
import itertools
import json
import math
from ztf_tap import get, RAW, ROOT
from ztf_screen import num,write_csv
def eligible():
    objects=collections.defaultdict(list)
    for f in sorted(RAW.glob('*_detections.csv.gz')):
        for row in csv.DictReader(gzip.open(f,'rt')):
            if num(row,'isdiffpos')==1 and 0<num(row,'sigmapsf')<=2.5/math.log(10)/5 and math.isfinite(num(row,'magpsf')):
                objects[row['oid']].append(row)
    result=[]
    for oid, rows in sorted(objects.items()):
        for a,b in itertools.combinations(sorted(rows,key=lambda r:num(r,'mjd')),2):
            if a['fid']==b['fid'] and 0.02<=num(b,'mjd')-num(a,'mjd')<=1 and a['pid']!=b['pid']:
                result.append((oid,a,b))
    return result
def key(r): return (str(r['fid']),round(float(r['mjd']),7))
def main():
    pairs=eligible(); ids=sorted({oid for oid,_,_ in pairs})
    batches=[ids[i:i+40] for i in range(0,len(ids),40)]
    def run(inp):
        batch,oids=inp
        q=("SELECT DISTINCT mjd,fid,programid FROM ztf.forced_photometry WHERE oid IN ("+
           ','.join("'"+oid+"'" for oid in oids)+") AND mjd>=61345 AND mjd<61361")
        tag='public_visit_batch_'+str(batch)
        path=RAW/(tag+'.csv.gz')
        if path.exists():
            return list(csv.DictReader(gzip.open(path,'rt')))
        return get(tag,q)
    programs=collections.defaultdict(set)
    with ThreadPoolExecutor(max_workers=3) as ex:
        for rows in ex.map(run,enumerate(batches)):
            for r in rows: programs[key(r)].add(r['programid'])
    needed={key(r) for oid,a,b in pairs for r in (a,b)}
    print('Needed',len(needed),'matched',len(needed.intersection(programs)),flush=True)
    unknown=needed.difference(programs)
    # Choose one available direct packet per visit independent of photometry.
    # Historical-only visits use the packet's parent to read the original program.
    lookup=collections.defaultdict(list)
    for f in sorted(RAW.glob('*_detections.csv.gz')):
        for r in csv.DictReader(gzip.open(f,'rt')):
            if key(r) in unknown: lookup[key(r)].append(r)
    targets={}
    for k in sorted(unknown):
        r=min(lookup[k],key=lambda row:(row['has_stamp']!='True',row['oid'],row['candid']))
        targets[k]=r
    hist=[r for r in targets.values() if r['has_stamp']!='True']
    parents={}
    if hist:
        q="SELECT oid,candid,parent_candid FROM ztf.detection WHERE oid IN ("+','.join("'"+r['oid']+"'" for r in hist)+')'
        path=RAW/'public_visit_parents.csv.gz'
        if path.exists():
            rows=list(csv.DictReader(gzip.open(path,'rt')))
        else:
            rows=get('public_visit_parents',q)
        for r in rows: parents[r['candid']]=r['parent_candid']
    folder=RAW/'program_packets'
    folder.mkdir(exist_ok=True)
    sources = {}
    for k, row in targets.items():
        cid = row['candid'] if row['has_stamp']=='True' else parents[row['candid']]
        sources[k] = (row['oid'], cid)

    def fetch(source):
        oid, cid = source
        url = 'https://avro.alerce.online/get_avro_info?' + urlencode(dict(oid=oid, candid=cid))
        path = folder / (cid + '.json.gz')
        receipt = folder / (cid + '.receipt.json')
        if not path.exists() or not receipt.exists():
            with urlopen(url, timeout=50) as reply:
                body = reply.read()
                status = reply.status
            path.write_bytes(gzip.compress(body, mtime=0))
            metadata = dict(url=url, method='GET', http_status=status,
                            sha256=hashlib.sha256(body).hexdigest(), bytes=len(body),
                            file=str(path.relative_to(ROOT)),
                            retrieved_utc=dt.datetime.now(dt.timezone.utc).isoformat())
            receipt.write_text(json.dumps(metadata) + '\n')
        with gzip.open(path, 'rt') as stream:
            return source, json.load(stream)

    with ThreadPoolExecutor(max_workers=4) as ex:
        packets = dict(ex.map(fetch, sorted(set(sources.values()))))
    for k, row in targets.items():
        packet = packets[sources[k]]
        photometry = [packet['candidate'], *(packet['prv_candidates'] or [])]
        epoch = next(e for e in photometry if str(e.get('candid')) == row['candid'])
        if (str(epoch['fid']), round(float(epoch['jd'])-2400000.5, 7)) != k:
            raise ValueError('Packet exposure does not match')
        programs[k].add(str(epoch['programid']))
    rows=[dict(fid=fid,mjd=mjd,programids=','.join(sorted(programs.get((fid,mjd),set())))) for fid,mjd in sorted(needed)]
    write_csv(ROOT/'derived'/'ztf_visits.csv',rows,rows[0].keys())
    accepted=[(oid,a,b) for oid,a,b in pairs if programs[key(a)]=={'1'} and programs[key(b)]=={'1'}]
    counts=dict(eligible_visits=len(needed),program1_visits=sum(programs[k]=={'1'} for k in needed),unverified_visits=sum(not programs[k] for k in needed),mixed_program_visits=sum(len(programs[k])>1 for k in needed),program1_pairs=len(accepted),program1_objects=len({oid for oid,a,b in accepted}))
    (ROOT/'derived'/'ztf_visit_counts.json').write_text(json.dumps(counts,indent=2)+'\n')
    print(counts)
if __name__=='__main__':main()
