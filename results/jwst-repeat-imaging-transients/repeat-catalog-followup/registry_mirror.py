"""Freeze regional rows of Kyoto's public TNS mirror and expose each page's coverage."""
import hashlib
from urllib.error import HTTPError
import re
import time
from concurrent.futures import ThreadPoolExecutor
from urllib.request import Request,urlopen

import numpy as np
import pandas as pd
from lxml import html
from counterparts import FIELDS
from campaign import HERE

BASE='https://www.o.kwasan.kyoto-u.ac.jp/~kawabata/tns/'
CENTERS=list(FIELDS.values())+[(268.53,65.2),(3.595,-30.4),(64.038,-24.068),(101.959,70.248)]
POSITION=re.compile(r'R\.A\.\s*=\s*(\d+):(\d+):([\d.]+),\s*Decl\.\s*=\s*([+\-])(\d+):(\d+):([\d.]+)')


def fetch(page):
    for n in range(4):
        try:
            with urlopen(Request(BASE+page,headers={'User-Agent':'Open-Space astronomical research'}),timeout=45) as f:
                data=f.read()
            return data
        except HTTPError as error:
            if error.code==404 or n==3: raise
            time.sleep(1+n)
        except OSError:
            if n==3:raise
            time.sleep(1+n)


def scan(page):
    data=fetch(page)
    tree=html.fromstring(data)
    rows=[];seen=parsed=0;names=[]
    for tr in tree.xpath('//tr[count(*[@rowspan])>=3]'):
        name=tr.xpath('string(*[1])').strip()
        if not re.fullmatch(r'20\d\d[A-Za-z]+',name):continue
        seen+=1; names.append(name)
        text=tr.xpath('string(*[@rowspan][3])')
        pos=POSITION.search(text)
        if not pos:continue
        h,m,s,sign,d,dm,ds=pos.groups()
        ra=15*(int(h)+int(m)/60+float(s)/3600)
        dec=(1 if sign=='+' else -1)*(int(d)+int(dm)/60+float(ds)/3600)
        parsed+=1
        if any(abs(dec-dd)<.6 and abs(ra-rr)*np.cos(np.deg2rad(dd))<.6 for rr,dd in CENTERS):
            rows.append((name,ra,dec,page))
    return rows,dict(file=page,url=BASE+page,sha256=hashlib.sha256(data).hexdigest(),bytes=len(data),
                     objects_seen=seen,positions_parsed=parsed,region_rows=len(rows),
                     first=names[0] if names else '',last=names[-1] if names else '')


def main():
    pages=[]
    for year in range(2026,2021,-1):
        page=f'{year}.html'
        try:
            doc=html.fromstring(fetch(page))
        except HTTPError as error:
            if error.code==404:
                continue
            raise
        pages+=sorted({x for x in doc.xpath('//div//a/@href') if re.fullmatch(fr'{year}[a-z]*\.html',x)})
    if not pages:raise ValueError('empty registry index')
    with ThreadPoolExecutor(8) as pool:
        output=list(pool.map(scan,pages))
    records=pd.DataFrame([r for rows,_ in output for r in rows],columns=['name','ra','dec','page'])
    receipts=pd.DataFrame([rec for _,rec in output])
    if (receipts.objects_seen==0).any():raise ValueError(receipts.loc[receipts.objects_seen==0,'file'].to_list())
    records.to_csv(HERE/'kyoto_field_positions.csv',index=False)
    receipts.to_csv(HERE/'kyoto_pages.csv',index=False)
    print(len(receipts),'pages,',receipts.positions_parsed.sum(),'position records,',len(records),'regional rows')


if __name__=='__main__':main()
