"""Freeze KMTNet annual inventory and download released pySIS event archives."""
import argparse
import concurrent.futures
import csv
import gzip
import hashlib
import io
import json
import re
import shutil
import tarfile
import time
import urllib.error
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from lxml import html

BASE = 'https://kmtnet.kasi.re.kr/ulens/event/'
CONTROLS = {'KMT-2024-BLG-0816', 'KMT-2024-BLG-0792', 'KMT-2024-BLG-3237'}
HERE = Path(__file__).resolve().parent


def digest(data):
    return hashlib.sha256(data).hexdigest()


def get(url):
    request = urllib.request.Request(url, headers={'User-Agent': 'Open-Space-research/1.0 (public microlensing inventory)'})
    with urllib.request.urlopen(request, timeout=65) as response:
        raw = response.read()
        return raw, dict(url=url, final_url=response.geturl(), status=response.status,
                         utc=response.headers.get('Date'), size=len(raw), sha256=digest(raw))


def index(cache):
    outputs, receipts, summary, seeds = [], [], [], []
    for year in (2024, 2025, 2026):
        url = f'{BASE}{year}/'
        path = cache / f'{year}.html'
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
        # 2026 is labeled proprietary; keep an aggregate and index receipt only.
        if public:
            dest = HERE / f'catalog-{year}.tsv'
            with dest.open('w') as out:
                w = csv.DictWriter(out, fieldnames=columns, delimiter='\t',lineterminator='\n')
                w.writeheader()
                w.writerows(rows)
            outputs.append(dest)
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
    if path.is_file() and receipt.is_file():
        return
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
