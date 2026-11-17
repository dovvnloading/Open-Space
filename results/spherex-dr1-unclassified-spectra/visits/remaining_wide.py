"""Complete native-image follow-up of every wide-field raw score >=5 winner."""
import argparse
import csv
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace

from measure import run

HERE=Path(__file__).resolve().parent
OUT=HERE/'wide_remaining'

def read(path):
    with path.open() as stream: return list(csv.DictReader(stream))

def save(path,rows):
    with path.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=rows[0],lineterminator='\n')
        writer.writeheader();writer.writerows(rows)

def select():
    wide=HERE.parent/'wide'
    done={(r['source_id'],r['channel']) for r in read(HERE/'wide/targets.csv')}
    spectra={(r['source_id'],r['channel']):r for r in read(wide/'spectra.csv')}
    selected=[]
    for r in read(wide/'full/ranks.csv'):
        k=(r['source_id'],r['peak_channel'])
        if float(r['score'])>=5 and k not in done:
            y=spectra[k]
            selected.append(dict(source_id=k[0],detector=y['detector'],wave_um=y['wave_um'],
                             scale_um=round(float(y['wave_hi_um'])-float(y['wave_lo_um']),6),
                             channel=k[1],selection='remaining_raw_score_ge5',label=''))
    return selected

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('index_cache',type=Path);p.add_argument('science_cache',type=Path)
    p.add_argument('--per-bin',type=int,default=3)
    p.add_argument('--workers-each',type=int,default=2,help='Concurrent cutouts in each of two target groups')
    p.add_argument('--retry-errors',action='store_true')
    p.add_argument('--replay',action='store_true')
    args=p.parse_args()
    OUT.mkdir(exist_ok=True)
    targets=select();save(OUT/'targets.csv',targets)
    if args.retry_errors:
        targets=[r for r in targets if not (OUT/(r['source_id']+'.json')).is_file() or
                 any(k.startswith('read_error') for k in json.loads((OUT/(r['source_id']+'.json')).read_text())['status_counts'])]
    # Preserve the initial / quality-fill labels if a transient read failure
    # occurs after the additional-exposure pass.
    previous={}
    if args.replay or args.retry_errors:
        for t in targets:
            path=OUT/(t['source_id']+'.csv')
            if path.is_file():
                saved={r['url']:r.get('selection_reason','') for r in read(path)}
                if any(saved.values()): previous[path]=saved
    with TemporaryDirectory() as tmpdir:
        jobs=[]
        for replay in (False,True):
            group=[t for t in targets if bool(args.replay or (args.retry_errors and
                   (OUT/(t['source_id']+'.csv')).is_file()))==replay]
            for i in range(2):
                part=group[i::2]
                if not part: continue
                file=Path(tmpdir)/f'targets_{replay}_{i}.csv';save(file,part)
                jobs.append(SimpleNamespace(targets=file,sources=HERE.parent/'wide/sources.csv',
                    index_cache=args.index_cache,calibration=HERE/'wide',out=OUT,science_cache=args.science_cache,
                    replay_dir=OUT if replay else None,per_bin=args.per_bin,workers=args.workers_each))
        with ThreadPoolExecutor(2) as pool: list(pool.map(run,jobs))
    for path,labels in previous.items():
        data=read(path)
        for row in data: row['selection_reason']=labels[row['url']]
        save(path,data)
    print('requested_pairs',len(targets),flush=True)

if __name__=='__main__':main()
