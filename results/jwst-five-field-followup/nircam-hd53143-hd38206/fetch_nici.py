"""Get independent Gemini/NICI H-channel images. Check the Gemini archive's byte hashes;
keep the raw-index receipt and a manifest with SHA-256 for both stored and decoded files."""
import argparse,bz2,csv,hashlib,json,time
from pathlib import Path
from urllib.request import urlopen
from datetime import datetime,timezone
from concurrent.futures import ThreadPoolExecutor
def main(cache,out):
    url='https://archive.gemini.edu/jsonsummary/20081218/NICI/'
    with urlopen(url,timeout=90) as f:raw=f.read()
    rows=[r for r in json.loads(raw) if r['object']=='HD_38206' and r['observation_id']=='GS-2008B-Q-500-32' and r['observation_class']=='science' and r['exposure_time']>=25]
    (out/'gemini_query.json').write_text(json.dumps(dict(url=url,sha256=hashlib.sha256(raw).hexdigest(),target='HD_38206',observation_id='GS-2008B-Q-500-32',min_exp_seconds=25,rows=len(rows),retrieved_utc=datetime.now(timezone.utc).isoformat()),indent=2)+'\n')
    def get(r):
        name=r['filename'];u='https://archive.gemini.edu/file/'+name;path=cache/name
        for i in range(3):
            try:
                if not path.exists():
                    with urlopen(u,timeout=65) as f:b=f.read()
                    assert hashlib.md5(b).hexdigest()==r['file_md5'];path.write_bytes(b)
                break
            except OSError:
                if i==2:raise
                time.sleep(2**i)
        b=path.read_bytes();d=bz2.decompress(b);assert hashlib.md5(d).hexdigest()==r['data_md5']
        return dict(file=name,url=u,compressed_SHA256=hashlib.sha256(b).hexdigest(),decompressed_SHA256=hashlib.sha256(d).hexdigest(),archive_md5=r['file_md5'],bytes=len(b),date=r['ut_datetime'],filter=r['filter_name'],exposure_time=r['exposure_time'])
    with ThreadPoolExecutor(6) as pool:
        results=list(pool.map(get,rows))
    with (out/'gemini_inputs.csv').open('w') as f:w=csv.DictWriter(f,list(results[0]));w.writeheader();w.writerows(results)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('cache',type=Path);p.add_argument('out',type=Path);a=p.parse_args();main(a.cache,a.out)
