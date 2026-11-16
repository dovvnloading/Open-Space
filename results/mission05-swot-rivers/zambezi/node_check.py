"""Diagnose the Chavuma 2025-04-18 release discrepancy at five nodes along the reach.

python node_check.py retrieve NODE_OUTPUT_DIRECTORY
python node_check.py compare NODE_OUTPUT_DIRECTORY
"""
import csv
import io
import json
import math
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from acquire import HYDROCRON, fetch, sha
from coverage import value
from reprocessing import write

FIELDS = ('node_id,time_str,lat,lon,p_wse,wse,wse_u,width,width_u,node_dist,node_q,node_q_b,'
          'n_good_pix,pass_id,cycle_id,granuleUR')
NODES = ('12297000160051','12297000160201','12297000160361','12297000160521','12297000160681')


def retrieve(out):
    out.mkdir(exist_ok=True, parents=True)
    def single(task):
        node, version = task
        query = dict(feature='Node',feature_id=node,start_time='2025-04-14T00:00:00Z',
                     end_time='2025-04-19T23:59:59Z',output='csv',
                     collection_name='SWOT_L2_HR_RiverSP_'+version,fields=FIELDS)
        url = HYDROCRON+'?'+urlencode(query)
        raw = fetch(url)
        content = json.loads(raw)['results']['csv'].encode()
        return dict(node_id=node, version=version, url=url, response_sha256=sha(raw), csv_sha256=sha(content)), content.decode()
    tasks = [(node,v) for node in NODES for v in ('2.0','D')]
    manifest, rows = [], []
    with ThreadPoolExecutor(max_workers=4) as pool:
        for record, content in pool.map(single,tasks):
            manifest.append(record)
            rows.extend(dict(version=record['version'],**r) for r in csv.DictReader(io.StringIO(content)))
    write(out/'node_observations.csv',rows)
    (out/'node_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')


def compare(out):
    with (out/'node_observations.csv').open() as f:
        rows = {(r['node_id'],r['pass_id'],r['version']): r for r in csv.DictReader(f)}
    result = []
    for node in NODES:
        for pid in ('208','305'):
            a,b = rows[node,pid,'2.0'],rows[node,pid,'D']
            valid = all(math.isfinite(value(r,'wse')) for r in (a,b))
            d = math.nan
            if valid:
                mean_lat = math.radians((value(a,'lat') + value(b,'lat'))/2)
                d = 111320*math.hypot((value(a,'lon') - value(b,'lon'))*math.cos(mean_lat),value(a,'lat')-value(b,'lat'))
            result.append(dict(node_id=node,pass_id=pid,date_v2=a['time_str'],date_d=b['time_str'],
                               q_v2=a['node_q'],q_d=b['node_q'],prior_height_v2=a['p_wse'],prior_height_d=b['p_wse'],
                               displacement_m=round(d,2) if valid else '',
                               d_minus_v2_wse_m=round(value(b,'wse')-value(a,'wse'),5) if valid else ''))
    write(out/'node_comparison.csv',result)


if __name__ == '__main__':
    mode,out=sys.argv[1:]
    if mode=='retrieve':
        retrieve(Path(out))
    elif mode!='compare':
        raise ValueError(mode)
    compare(Path(out))
