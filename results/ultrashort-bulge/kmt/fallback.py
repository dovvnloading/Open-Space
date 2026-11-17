"""Acquire advertised DIA data when the pySIS archive is unavailable."""
import argparse
import concurrent.futures
import io
import json
import tarfile
from pathlib import Path
from urllib.parse import urljoin
from lxml import html

from freeze import get, HERE


def acquire(row,cache):
    name=row['event']
    path=cache/(name+'.tar.gz')
    receipt=cache/(name+'.json')
    if path.exists() and receipt.exists():
        return
    rec={'event':name}
    try:
        page,rec['page']=get(row['page'])
        links=html.fromstring(page).xpath('//a/@href')
        dia=next((s for s in links if s.endswith('/diapl.tar.gz')),None)
        if dia is None:
            rec['status']='no_DIA_archive_link_on_page'
        else:
            raw,rec['curve']=get(urljoin(row['page'],dia))
            with tarfile.open(fileobj=io.BytesIO(raw),mode='r:gz') as f:
                members=[m.name for m in f if m.isfile() and m.name.endswith('.diapl')]
            if not members:
                raise ValueError('no KMT DIA curve members')
            rec.update(status='ok',members=members)
            temp=path.with_suffix('.part')
            temp.write_bytes(raw);temp.replace(path)
    except (OSError,ValueError,tarfile.TarError) as exc:
        rec.update(status='failed',error=str(exc))
    receipt.write_text(json.dumps(rec)+'\n')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cache',type=Path,required=True)
    a=p.parse_args()
    target=a.cache/'fallback';target.mkdir(exist_ok=True)
    rows={r['event']:r for r in map(json.loads,(a.cache/'index.jsonl').read_text().splitlines())}
    missing=[p.name.split('.')[0] for p in a.cache.glob('*.error.json') if json.loads(p.read_text()).get('status')==404 and not (a.cache/(p.name.split('.')[0]+'.tar.gz')).exists()]
    with concurrent.futures.ThreadPoolExecutor(3) as ex:
        list(ex.map(lambda n:acquire(rows[n],target),missing))
    entries=[json.loads(p.read_text()) for p in sorted(target.glob('*.json'))]
    (HERE/'fallback-receipts.json').write_text(json.dumps(entries,indent=2)+'\n')


if __name__=='__main__':
    main()
