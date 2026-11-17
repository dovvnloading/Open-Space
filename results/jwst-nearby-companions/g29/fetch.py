"""Download just the two GO 9077 mosaics and four date-bracketing SUB256 standards."""
import argparse
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen
import sys

HERE = Path(__file__).resolve().parent
def main(cache, freeze):
    cache.mkdir(parents=True, exist_ok=True)
    if freeze:
        # Re-run the exact archived CAOM requests (saved as plain URLs) and freeze response digests.
        from datetime import datetime, timezone
        queries = json.loads((HERE / 'query_urls.json').read_text())
        inputs, receipts = [], {}
        for name, url in queries.items():
            raw = urlopen(url, timeout=90).read()
            result = json.loads(raw)
            if result['status'] != 'COMPLETE' or len(result['data']) >= 3000:
                raise ValueError('CAOM result incomplete')
            rows = result['data']
            receipts[name] = dict(url=url, sha256=hashlib.sha256(raw).hexdigest(),
                                  rows=len(rows), retrieved_utc=datetime.now(timezone.utc).isoformat())
            if name == 'science':
                inputs.extend(rows)
            else:
                # Closest observation preceding and following the target in each band.
                tm = inputs[0]['t_min']
                for band in ('F2100W', 'F2550W'):
                    rr = [r for r in rows if r['filters']==band and r['obs_id'].endswith('sub256')]
                    inputs += [max((r for r in rr if r['t_min']<tm),key=lambda r:r['t_min']),
                               min((r for r in rr if r['t_min']>tm),key=lambda r:r['t_min'])]
        for row in inputs:
            obs = row['obs_id']
            program,visit = obs.split('_')[0].split('-')
            row['filename'] = obs+'_i2d.fits'
            row['url'] = f'https://stpubdata.s3.amazonaws.com/jwst/public/{program}/L3/t/{visit}/{row["filename"]}'
    else:
        inputs = json.loads((HERE/'inputs.json').read_text())
    def fetch(row):
        p = cache / row['filename']
        if not p.exists():
            raw = urlopen(row['url'],timeout=120).read()
            p.write_bytes(raw)
        dig = hashlib.sha256(p.read_bytes()).hexdigest()
        if row.get('sha256',dig)!=dig:
            raise ValueError(f'Hash mismatch: {p}')
        return dict(row,sha256=dig,bytes=p.stat().st_size)
    with ThreadPoolExecutor(4) as ex:
        checked=list(ex.map(fetch,inputs))
    if freeze:
        for name,obj in [('inputs',checked),('mast_receipts',receipts)]:
            (HERE/(name+'.json')).write_text(json.dumps(obj,indent=2)+'\n')
if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cache',type=Path); parser.add_argument('--freeze',action='store_true')
    a=parser.parse_args(); main(a.cache,a.freeze)
