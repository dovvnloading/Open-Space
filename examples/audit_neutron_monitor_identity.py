"""Audit station identity corruption in public neutron monitor Parquet.

Audits labels only; it does not assert that plausible numeric rate series
are invented. An HTML DOCTYPE is never a neutron-monitor station ID.
"""
from __future__ import annotations
from collections import Counter
from pathlib import Path
import argparse
import hashlib
import json
from urllib.request import Request, urlopen

URL="https://huggingface.co/datasets/juliensimon/neutron-monitor-cosmic-rays/resolve/refs%2Fconvert%2Fparquet/default/train/0000.parquet"
# Exact fragments of a W3C HTML 4.01 DTD declaration, not station IDs.
HTML_HEADER_TOKENS=frozenset(('" -//W3C//DTD'.strip().replace(' ',''),'4.01//EN"','<!DOCTYPE','HTML','PUBLIC','"-//W3C//DTD'))

def looks_like_html_header(value: object) -> bool:
    return value in HTML_HEADER_TOKENS if isinstance(value,str) else False

def summarize(rows: list[dict]) -> dict:
    ids=Counter(str(r['station']) for r in rows)
    names=Counter(str(r['station_name']) for r in rows)
    header_rows=sum(n for code,n in ids.items() if looks_like_html_header(code))
    missing_rates=sum(r.get('count_rate') is None for r in rows)
    return {
        'source':URL, 'row_count':len(rows),
        'distinct_station_ids':len(ids),
        'station_id_counts':dict(sorted(ids.items())),
        'station_name_counts':dict(sorted(names.items())),
        'html_doctype_fragment_rows':header_rows,
        'html_doctype_fragment_fraction':round(header_rows/len(rows),8) if rows else None,
        'identical_station_name_to_station_id_rows':sum(r.get('station')==r.get('station_name') for r in rows),
        'missing_count_rate_rows':missing_rates,
        'data_quality_conclusion':('station_id_provenance_unusable' if header_rows else 'no_html_header_signature_detected'),
        'caution':'Station labels are invalid; this does NOT prove the count-rate values are fraudulent or synthetic.'
    }

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fetch',action='store_true')
    p.add_argument('--input',type=Path,default=Path('/tmp/neutron-monitor.parquet'))
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.fetch:
        req=Request(URL,headers={'User-Agent':'Open-Space-research/0.2'})
        with urlopen(req,timeout=180) as src,a.input.open('wb') as dst:
            while chunk:=src.read(2*1024*1024):dst.write(chunk)
    import pyarrow.parquet as pq
    rows=pq.read_table(a.input,columns=['station','station_name','count_rate']).to_pylist()
    result=summarize(rows)
    result['input_sha256']=hashlib.sha256(a.input.read_bytes()).hexdigest()
    result['input_bytes']=a.input.stat().st_size
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=='__main__':
    main()
