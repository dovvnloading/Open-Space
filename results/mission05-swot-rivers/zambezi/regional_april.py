"""Replay April 2025 Chavuma neighborhood release diagnostic from its query manifest.

python regional_april.py DATA_DIRECTORY D_DIRECTORY
"""
import csv
import io
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from acquire import fetch, sha
from reprocessing import write


def main(data, out):
    manifest = json.loads((out/'regional_april_manifest.json').read_text())
    def get(record):
        reach = record['translation']['v16_reach_id']
        raw = fetch(record['url'])
        content = json.loads(raw)['results']['csv']
        with (data/f'reach_{reach}.csv').open() as f:
            v2 = {r['pass_id']:r for r in csv.DictReader(f) if r['cycle_id']=='031'}
        results = []
        for row in csv.DictReader(io.StringIO(content)):
            if row['pass_id'] not in v2:
                continue
            a = v2[row['pass_id']]
            fields = ['reach_id','pass_id','time_str','granuleUR','wse','width','reach_q','node_dist','loc_offset','p_dist_out']
            results.append(dict(**{k:row[k] for k in fields}, v2_wse=a['wse'],v2_width=a['width'],
                                v2_reach_q=a['reach_q'],v2_granule=a['granuleUR']))
        return dict(**record,csv_payload_sha256=sha(content.encode())),results
    vals, urls = [], []
    with ThreadPoolExecutor(max_workers=4) as pool:
        for record, rows in pool.map(get,manifest):
            urls.append(record)
            vals.extend(rows)
    write(out/'regional_april_comparison.csv',vals)
    (out/'regional_april_payloads.json').write_text(json.dumps(urls,indent=2)+'\n')


if __name__=='__main__':
    main(*map(Path,sys.argv[1:]))
