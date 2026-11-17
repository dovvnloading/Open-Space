"""Batch 1.5-arcsecond crossmatches for the complete released KMT alert inventory, via VizieR ASU."""
import argparse
import concurrent.futures
import csv
import json
import re
import time
import urllib.request
from pathlib import Path
from freeze import HERE,digest
from neighbours import angle

ENDPOINT='https://vizier.cfa.harvard.edu/viz-bin/asu-tsv'
CATALOGS={
 'VIRAC2':('II/387/virac2','srcid,RAJ2000,DEJ2000,pmRA,pmDE,Jmag,Hmag,Ksmag,e_Ksmag,o_Ksmag,MedKsmag,etaKs,KsEp1,KsEpl'),
 'GaiaDR3':('I/355/gaiadr3','Source,RA_ICRS,DE_ICRS,pmRA,pmDE,Gmag,BP-RP,VarFlag,RUWE')}
BATCH=80
# The frozen request batches each included one source from the other year.
# Replay its actual coordinates for hash reproducibility and omit the duplicate
# from the year-specific scientific match table.
REPEATED_TARGET={2024:'KMT-2025-BLG-2024',2025:'KMT-2024-BLG-2025'}


def query(job):
    db,rows,year,i,cache,cached_only=job
    name=f'{db}-{year}-{i:04d}'
    cat,cols=CATALOGS[db]
    coords=[f"{angle(r['ra'],True):.9f} {angle(r['dec']):+.9f}" for r in rows]
    targets='\n'.join(coords)
    payload=(f'-source={cat}\n-out.max=10000\n-out.add=_q,_r\n-out={cols}\n-c.rs=1.5\n'
             f'-c=<<======EndList\n{targets}\n======EndList\n').encode()
    rawfile=cache/(name+'.raw'); receipt=cache/(name+'.json')
    if rawfile.exists() and receipt.exists():
        raw=rawfile.read_bytes(); rec=json.loads(receipt.read_text())
        if rec.get('request_sha256')!=digest(payload) or rec.get('sha256')!=digest(raw):
            raise ValueError(f'Cached VizieR request or response mismatch in {name}')
        if rec.get('status')==200:rec.pop('error',None)
    else:
        rec=dict(batch=name,source=cat,url=ENDPOINT,method='POST',request_sha256=digest(payload),
                 first_event=rows[0]['event'],last_event=rows[-1]['event'],targets=len(rows))
        if cached_only:
            if receipt.exists():rec=json.loads(receipt.read_text())
            else:rec['status']='pending'
            return db,year,[],rec
        request=urllib.request.Request(ENDPOINT,data=payload,headers={'Content-Type':'application/x-www-form-urlencoded',
                                                                     'User-Agent':'Open-Space-research/1.0'})
        for attempt in range(3):
            try:
                with urllib.request.urlopen(request,timeout=140) as response:
                    raw=response.read();rec.update(status=response.status,utc=response.headers.get('date'),
                                                  bytes=len(raw),sha256=digest(raw))
                if b'#Column\t_q' not in raw and b'No object found' not in raw:
                    raise ValueError('no VizieR result table or empty indicator')
                rec.pop('error',None)
                rawfile.write_bytes(raw);break
            except (OSError,ValueError) as exc:
                rec.update(status='failed',error=str(exc));raw=b''
                if attempt<2:time.sleep(2+attempt)
        receipt.write_text(json.dumps(rec)+'\n')
    lines=raw.decode(errors='replace').splitlines();fields=None;matches=[]
    for line in lines:
        if line.startswith('_q\t'):
            fields=line.split('\t')
        elif re.match(r'^\s*\d+\t',line) and fields:
            values=[s.strip() for s in line.split('\t')]
            r=dict(zip(fields,values));q=int(r.pop('_q'));r['sep_arcsec']=r.pop('_r')
            if q<1 or q>len(rows):raise ValueError((name,q))
            matches.append(dict(event=rows[q-1]['event'],**r))
    # The VizieR server echoes every target evaluated, including empty cones.
    rec['echoed_cones']=sum(l.startswith('#INFO\t-c=') for l in lines)
    rec['matches']=len(matches)
    if rec.get('status')==200 and rec['echoed_cones']!=len(rows):
        rec['coverage_warning']='VizieR did not echo every target; unmatched IDs unverified'
    return db,year,matches,rec


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cache',type=Path,required=True)
    p.add_argument('--workers',type=int,default=3)
    p.add_argument('--cached-only',action='store_true',help='Export verified results acquired so far without network access')
    args=p.parse_args()
    data=[json.loads(x) for x in (args.cache/'index.jsonl').read_text().splitlines()]
    out=args.cache/'batch';out.mkdir(exist_ok=True)
    jobs=[]
    for year in (2024,2025):
        rows=sorted([r for r in data if r['event'][4:8]==str(year) or r['event']==REPEATED_TARGET[year]],
                    key=lambda r:r['event'])
        for i in range(0,len(rows),BATCH):
            for db in CATALOGS:jobs.append((db,rows[i:i+BATCH],year,i,out,args.cached_only))
    allrows={};receipts=[]
    with concurrent.futures.ThreadPoolExecutor(args.workers) as pool:
        for db,year,rows,rec in pool.map(query,jobs):
            receipts.append(rec)
            allrows.setdefault((db,year),[]).extend(r for r in rows if r['event'][4:8]==str(year))
    for (db,year),rows in allrows.items():
        if not rows:continue
        cols=list(dict.fromkeys(k for r in rows for k in r))
        with (HERE/f'{db}-{year}.tsv').open('w') as h:
            w=csv.DictWriter(h,cols,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    (HERE/'batch-neighbour-receipts.json').write_text(json.dumps(receipts,indent=2)+'\n')


if __name__=='__main__':
    main()
