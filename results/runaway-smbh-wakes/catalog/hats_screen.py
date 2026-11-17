"""Apply the frozen Q1 shape cuts to every partition of IRSA's merged HATS catalogue."""
import argparse
import csv
import hashlib
import io
import json
import os
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from xml.etree import ElementTree as ET

import numpy as np
import pyarrow as pa
import pyarrow.fs as afs
import pyarrow.parquet as pq

HERE = Path(__file__).resolve().parent
S3 = 'https://nasa-irsa-euclid-q1.s3.us-east-1.amazonaws.com/'
PREFIX = 'contributed/q1/merged_objects/hats/euclid_q1_merged_objects-hats/'
CUTS = ['mer_vis_det', 'mer_ellipticity', 'mer_semimajor_axis']
COLUMNS = ['object_id','tileid','ra','dec',*CUTS,'mer_position_angle','mer_kron_radius',
           'mer_segmentation_area','mer_det_quality_flag','phz_phz_median']
SHA = lambda b: hashlib.sha256(b).hexdigest()

def download(url, headers=None):
    for attempt in range(4):
        try:
            with urlopen(Request(url, headers=headers or {}), timeout=75) as response:
                return response.status, response.headers, response.read()
        except OSError:
            if attempt == 3: raise
            time.sleep(2**attempt)

class RangeFile(io.RawIOBase):
    def __init__(self,url,size,etag):
        self.url,self.size,self.etag,self.pos=url,size,etag,0
        self.receipts=[]
    def readable(self): return True
    def seekable(self): return True
    def tell(self): return self.pos
    def seek(self,offset,whence=0):
        self.pos = [0,self.pos,self.size][whence]+offset
        return self.pos
    def read(self,n=-1):
        if n<0: n=self.size-self.pos
        n=min(n,self.size-self.pos)
        if n<=0: return b''
        a,b=self.pos,self.pos+n-1
        # Unique URLs prevent intermediary caches that ignore Range from mixing responses.
        url=self.url+f'?range_start={a}&range_end={b}'
        status,h,raw=download(url,{'Range':f'bytes={a}-{b}','If-Match':self.etag,'Accept-Encoding':'identity'})
        if status!=206 or h.get('Content-Range')!=f'bytes {a}-{b}/{self.size}' or len(raw)!=n:
            raise IOError(f'Unexpected range response for {url}')
        self.receipts.append([a,b,SHA(raw)])
        self.pos+=n
        return raw

def inspect(item,cache,filesystem=None):
    key,size,etag=item
    name=key.split('/Npix=')[1].replace('/','_')+'.json'
    part=cache/name
    if part.exists():
        prev=json.loads(part.read_text())
        if prev['etag']==etag:return prev
    if filesystem is None:
        f=RangeFile(S3+key,size,etag)
    else:
        f=filesystem.open_input_file('nasa-irsa-euclid-q1/'+key)
        if f.size()!=size: raise IOError('Parent size changed: '+key)
    pf=pq.ParquetFile(f,pre_buffer=filesystem is not None)
    projected=[]
    if filesystem is not None:
        allcols=pf.read(columns=COLUMNS,use_threads=False)
        # Hash decoded Arrow columns, separate from the publisher object ETag.
        for c in COLUMNS:
            arr=allcols[c].combine_chunks()
            sink=pa.BufferOutputStream()
            batch=pa.record_batch([arr],names=[c])
            with pa.ipc.new_stream(sink,batch.schema) as writer: writer.write_batch(batch)
            projected.append([c,hashlib.sha256(sink.getvalue()).hexdigest()])
    rows=[];nvis=0
    for i in range(1 if filesystem is not None else pf.metadata.num_row_groups):
        t=allcols if filesystem is not None else pf.read_row_group(i, columns=CUTS,use_threads=False)
        vis,ell,major=[np.asarray(t[c].to_numpy()) for c in CUTS]
        nvis+=int((vis==1).sum())
        selected=np.flatnonzero((vis==1)&(ell>0.95)&(ell<1)&(major>=7.5))
        if len(selected):
            if filesystem is not None:
                extra=t
            else:
                extra=pf.read_row_group(i,columns=[c for c in COLUMNS if c not in CUTS],use_threads=False)
                for c in CUTS: extra=extra.append_column(c,t[c])
            rows.extend(extra.select(COLUMNS).take(selected).to_pylist())
    record=dict(key=key,bytes=size,etag=etag,groups=pf.metadata.num_row_groups,n=pf.metadata.num_rows,
                n_vis=nvis,n_selected=len(rows),ranges=getattr(f,"receipts",[]),projected=projected,selected=rows)
    f.close()
    part.write_text(json.dumps(record,allow_nan=False)+'\n')
    print(name,record['n'],len(rows),flush=True)
    return record

def table(path,rows,fields):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cache',type=Path,required=True)
    p.add_argument('--workers',type=int,default=6)
    p.add_argument('--transport',choices=['native','ranges'],default='native')
    a=p.parse_args();a.cache.mkdir(exist_ok=True,parents=True)
    listurl=S3+'?'+urlencode({'list-type':2,'prefix':PREFIX+'dataset/Norder'})
    _,_,raw=download(listurl);root=ET.fromstring(raw)
    assert root.find('{*}IsTruncated').text=='false'
    items=[(x.find('{*}Key').text,int(x.find('{*}Size').text),x.find('{*}ETag').text) for x in root.findall('{*}Contents')]
    _,_,meta=download(S3+PREFIX+'partition_info.csv')
    expected={f"Norder={r['Norder']}/Dir={r['Dir']}/Npix={r['Npix']}":int(r['num_rows']) for r in csv.DictReader(io.StringIO(meta.decode()))}
    assert len(items)==len(expected)
    native=None
    if a.transport=='native':
        proxy=os.environ.get('HTTPS_PROXY',os.environ.get('https_proxy',''))
        native=afs.S3FileSystem(anonymous=True,region='us-east-1',proxy_options=proxy or None,
                                request_timeout=90,connect_timeout=20)
    def attempt(item):
        for n in range(4):
            try:return inspect(item,a.cache,native if n<2 else None)
            except OSError:
                if n==3:raise
                time.sleep(7*2**n)
    with ThreadPoolExecutor(a.workers) as pool: records=list(pool.map(attempt,items))
    records.sort(key=lambda r:r['key'])
    for r in records:assert expected[r['key'].split('dataset/')[1].rsplit('/',1)[0]]==r['n']
    rows=sorted([s for r in records for s in r['selected']],key=lambda s:s['object_id'])
    table(HERE/'hats_selected.csv',rows,COLUMNS)
    table(HERE/'hats_partitions.csv',[{k:r[k] for k in ('key','bytes','etag','n','n_vis','n_selected')} for r in records],
          ['key','bytes','etag','n','n_vis','n_selected'])
    import gzip
    with (HERE/'hats_ranges.jsonl.gz').open('wb') as dest, gzip.GzipFile(fileobj=dest,mode='wb',mtime=0) as gz:
        for r in records:
            gz.write((json.dumps({'key':r['key'],'ranges':r['ranges'],'projected_ipc_sha256':r.get('projected',[])},separators=(',',':'))+'\n').encode())
    report=dict(list_url=listurl,list_sha256=SHA(raw),partition_info_url=S3+PREFIX+'partition_info.csv',partition_info_sha256=SHA(meta),
                expected_partitions=len(expected),partitions=len(records),n=sum(r['n'] for r in records),
                n_vis=sum(r['n_vis'] for r in records),n_selected=len(rows),unique_selected_ids=len({s['object_id'] for s in rows}))
    (HERE/'hats_summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(report,flush=True)

if __name__=='__main__':main()
