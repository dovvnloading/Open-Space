"""Freeze KMTNet annual inventory and download released pySIS event archives."""
import argparse
import concurrent.futures
import csv
import shutil
import hashlib
import gzip
import io
import json
import re
import tarfile
import time
import threading
from email.utils import parsedate_to_datetime
import urllib.error
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from lxml import html

BASE = 'https://kmtnet.kasi.re.kr/ulens/event/'
CONTROLS = {'KMT-2024-BLG-0816', 'KMT-2024-BLG-0792', 'KMT-2024-BLG-3237'}
HERE = Path(__file__).resolve().parent
_RATE_GATE=threading.Lock()
_RESUME_AT=0.0


def digest(data):
    return hashlib.sha256(data).hexdigest()


def get(url):
    global _RESUME_AT
    while True:
        with _RATE_GATE:
            wait=_RESUME_AT-time.monotonic()
        if wait<=0: break
        time.sleep(min(60,wait))
    request = urllib.request.Request(url, headers={'User-Agent': 'Open-Space-research/1.0 (public microlensing inventory)'})
    try:
        opened=urllib.request.urlopen(request,timeout=65)
    except urllib.error.HTTPError as exc:
        if exc.code==429:
            header=exc.headers.get('Retry-After','60')
            try:
                pause=float(header)
            except ValueError:
                try:pause=(parsedate_to_datetime(header)-datetime.now(timezone.utc)).total_seconds()
                except (ValueError,TypeError):pause=60.
            with _RATE_GATE:
                _RESUME_AT=max(_RESUME_AT,time.monotonic()+max(30.,pause))
        raise
    with opened as response:
        raw = response.read()
        return raw, dict(url=url, final_url=response.geturl(), status=response.status,
                         utc=response.headers.get('Date'), size=len(raw), sha256=digest(raw))


def index(cache):
    receipts, summary, seeds = [], [], []
    for year in (2024, 2025, 2026):
        url = f'{BASE}{year}/'
        path = cache / f'{year}.html'
        frozen=HERE/f'index-{year}.html.gz'
        # Reproduction can use the frozen public 2026 roster and aggregate receipt
        # without re-reading an alert service that continues to change.
        if year==2026 and not path.exists() and (HERE/'event-pages-2026.tsv').exists() and (HERE/'index-receipts.json').exists():
            prior=json.loads((HERE/'index-receipts.json').read_text())
            receipts.append(next(r for r in prior['queries'] if r['url']==url))
            summary.append(next(r for r in prior['seasons'] if r['year']==year))
            continue
        if not path.exists() and frozen.exists() and (HERE/'index-receipts.json').exists():
            original=next(v for v in json.loads((HERE/'index-receipts.json').read_text())['queries'] if v['url']==url)
            raw=gzip.decompress(frozen.read_bytes())
            if digest(raw)!=original['sha256']:raise ValueError(f'Frozen index mismatch: {year}')
            path.write_bytes(raw)
            (cache/f'{year}.headers').write_text('date: '+original['utc']+'\n')
        if path.exists():
            raw = path.read_bytes()
            header = (cache / f'{year}.headers').read_text()
            dates = re.findall(r'^date: (.+)$', header, re.M | re.I)
            rec = dict(url=url, final_url=url, status=200, utc=dates[-1].strip(), size=len(raw), sha256=digest(raw))
        else:
            raw, rec = get(url)
            path.write_bytes(raw)
            (cache / f'{year}.headers').write_text('date: '+rec['utc']+'\n')
        doc = html.fromstring(raw)
        trs = doc.xpath('//tr')
        columns = [x.text_content().strip() for x in trs[0].xpath('./td')]
        rows = []
        for tr in trs[1:]:
            cells = [x.text_content().strip() for x in tr.xpath('./td')]
            if len(cells) != len(columns):
                raise ValueError((year, cells))
            rows.append(dict(zip(columns, cells)))
        reported = doc.xpath('//h1')[0].text_content()
        public = '(public)' in reported
        if public != (year in (2024, 2025)):
            raise ValueError(f'Check new season release status: {reported}')
        vals = [float(r['t_E']) for r in rows if r['t_E'] not in ('','-') and float(r['t_E'])>0]
        counts = Counter(re.search(r'\d{2}',r['Field']).group() for r in rows)
        summary.append(dict(year=year, release_label=reported, listed=len(rows), positive_tE=len(vals),
                            smallest_positive_tE_d=min(vals), below_0p50=sum(x<.50 for x in vals),
                            at_0p50=sum(x==.50 for x in vals), at_or_below_2=sum(x<=2 for x in vals),
                            fields=dict(sorted(counts.items()))))
        receipts.append(rec)
        if not public:
            with (HERE/f'event-pages-{year}.tsv').open('w') as h:
                w=csv.writer(h,delimiter='\t',lineterminator='\n')
                w.writerow(('event','event_page'))
                for r in rows:
                    w.writerow((r['Event'],f"{url}view.php?event={r['Event']}"))
        if public:
            if not frozen.exists():frozen.write_bytes(gzip.compress(raw,compresslevel=9,mtime=0))
            dest = HERE / f'catalog-{year}.tsv'
            with dest.open('w') as out:
                w = csv.DictWriter(out, fieldnames=columns, delimiter='\t',lineterminator='\n')
                w.writeheader()
                w.writerows(rows)
            for r in rows:
                name = r['Event']
                code = f'KB{year%100:02d}{name[-4:]}'
                forkey = lambda k: float(r[k]) if r[k] not in ('','-') else None
                seeds.append(dict(event=name, field=r['Field'], ra=r['RA'], dec=r['Dec'],
                                  related=r['Related event'],t0=forkey('t_0'), te=forkey('t_E'),
                                  u0=forkey('u_0'), url=f'{url}data/{code}/pysis/pysis.tar.gz',
                                  page=f'{url}view.php?event={name}'))
    (HERE/'index-receipts.json').write_text(json.dumps(dict(queries=receipts,seasons=summary), indent=2)+'\n')
    with (cache/'index.jsonl').open('w') as out:
        for r in seeds:
            out.write(json.dumps(r)+'\n')
    return seeds


def fetch(row, cache):
    event = row['event']
    url = row['url']
    path = cache / f'{event}.tar.gz'
    receipt = cache / f'{event}.receipt.json'
    if receipt.is_file() and (path.is_file() or (cache / f'{event}.evicted').is_file()):
        return
    while shutil.disk_usage(cache).free < 120*1024**2:
        time.sleep(10)
    for attempt in range(3):
        try:
            if path.exists():
                data = path.read_bytes()
                dates = re.findall(r'^date: (.+)$', (cache/f'{event}.headers').read_text(), re.M|re.I)
                rec = dict(url=url,final_url=url,status=200,utc=dates[-1].strip(),size=len(data),sha256=digest(data))
            else:
                data, rec = get(url)
            with tarfile.open(fileobj=io.BytesIO(data),mode='r:gz') as f:
                curves = [m.name for m in f if m.isfile() and m.name.endswith('.pysis')]
                if not curves:
                    raise ValueError('no pySIS light curves')
            rec.update(event=event,curves=curves)
            temp = path.with_suffix('.part')
            temp.write_bytes(data)
            temp.replace(path)
            receipt.write_text(json.dumps(rec)+'\n')
            (cache/f'{event}.error.json').unlink(missing_ok=True)
            break
        except (OSError,tarfile.TarError,ValueError) as exc:
            if attempt == 2 or (isinstance(exc,urllib.error.HTTPError) and exc.code == 404):
                failure = dict(event=event,url=url,status=getattr(exc,'code',None),error=str(exc),
                               utc=datetime.now(timezone.utc).isoformat())
                (cache/f'{event}.error.json').write_text(json.dumps(failure)+'\n')
                break
            time.sleep(1+attempt)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cache',type=Path,required=True)
    p.add_argument('--download',action='store_true')
    p.add_argument('--jobs',type=int,default=4)
    p.add_argument('--limit',type=int)
    a=p.parse_args()
    a.cache.mkdir(exist_ok=True,parents=True)
    rows=index(a.cache)
    if a.download:
        ordered=sorted(rows,key=lambda r:(r['event'] not in CONTROLS,
                       r['te'] if r['te'] and r['te']>0 else float('inf'),r['event']))[:a.limit]
        with concurrent.futures.ThreadPoolExecutor(a.jobs) as ex:
            list(ex.map(lambda r:fetch(r,a.cache),ordered))


if __name__=='__main__':
    main()
