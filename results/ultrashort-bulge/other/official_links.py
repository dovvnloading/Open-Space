"""Freeze aliases/positional neighbours explicitly supplied by the MOA/PRIME team.

The two numerical columns are retained without imposing a cut. The list can include
chance neighbours from other years and requires per-event position/epoch validation.
"""
from concurrent.futures import ThreadPoolExecutor
import csv
from pathlib import Path
from fetch import get, SHA, write_csv, OUT, BASE

CACHE=Path('.survey-cache')
def one(pair):
    url=BASE+'stuff/datafile/xref/'+pair+'.txt'
    path=CACHE/pair
    if not path.exists():path.write_bytes(get(url))
    raw=path.read_bytes()
    index={r['event'] for r in csv.DictReader((OUT/'index.csv').open())}
    rows=[]
    for line in raw.decode().splitlines():
        bits=line.split()
        if len(bits)==4 and (bits[0] in index or bits[1] in index):
            rows.append(dict(pair=pair,id_1=bits[0],id_2=bits[1],survey_column_3=bits[2],survey_column_4=bits[3]))
    return dict(pair=pair,url=url,bytes=len(raw),sha256=SHA(raw),rows_in_scope=len(rows)),rows

if __name__=='__main__':
    with ThreadPoolExecutor(5) as ex:data=list(ex.map(one,['moa-kmt','prime-kmt','moa-ogle','prime-ogle','prime-moa']))
    write_csv(OUT/'official_link_receipts.csv',[r for r,_ in data])
    write_csv(OUT/'official_links.csv',[r for _,rows in data for r in rows])
