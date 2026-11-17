"""Run the held-out Sector 105 corrected-light-curve search with a signed dip gate."""
import argparse
import csv
import gzip
import hashlib
import io
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from acquire import fetch, product_url, select
from scan import AUX_COLUMNS, Curve, rates, scan_one, write_csv


def header(raw, offset):
    """Parse only scalar keywords in standard header cards used by these QLP files."""
    out = {}
    for i in range(offset, len(raw), 80):
        card = raw[i:i+80].decode('ascii')
        if card.strip() == 'END':
            return out, ((i + 80 + 2879)//2880)*2880
        if card[8:10] == '= ':
            out[card[:8].strip()] = card[10:].split('/')[0].strip().strip("'").strip()
    raise ValueError('missing FITS END card')


def unpack(raw):
    h, offset = header(raw,0)
    if h.get('SIMPLE') != 'T' or h.get('NAXIS') != '0':
        raise ValueError('unexpected primary HDU')
    ext, offset = header(raw,offset)
    if ext.get('XTENSION') != 'BINTABLE':
        raise ValueError('expected QLP binary light curve table')
    columns = []
    dtypes = {'E': '>f4', 'D': '>f8', 'J': '>i4'}
    for i in range(1, int(ext['TFIELDS'])+1):
        if float(ext.get(f'TSCAL{i}',1)) != 1 or float(ext.get(f'TZERO{i}',0)) != 0:
            raise ValueError('scaled FITS column requires a general FITS reader')
        columns.append((ext[f'TTYPE{i}'], dtypes[ext[f'TFORM{i}']]))
    dtype = np.dtype(columns)
    if int(ext['NAXIS1']) != dtype.itemsize or len(raw) % 2880:
        raise ValueError('invalid QLP product length')
    return h, np.frombuffer(raw, dtype=dtype, count=int(ext['NAXIS2']), offset=offset)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pilot-sample', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--workers', type=int, default=6)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    cache = args.out/'trigger-inputs'
    cache.mkdir(exist_ok=True)
    pilot = json.loads(args.pilot_sample.read_text())
    raw = fetch(pilot['url'])
    digest = hashlib.sha256(raw).hexdigest()
    if digest != pilot['sha256'] or pilot['sector'] != 105:
        raise ValueError('target list changed')
    _, _, ranked = select(raw, 384)
    if ranked[:128] != pilot['targets']:
        raise ValueError('first 128 differ from frozen pilot')
    targets = ranked[128:384]
    fields = ['rank', 'tic', 'sector', 'ra', 'dec', 'url', 'bytes', 'sha256', 'status', 'detail']
    def load(target):
        i, star = target
        tic = star['tic']
        url = product_url(105,tic)
        receipt = dict(rank=i,tic=tic,sector=105,ra=star['ra'],dec=star['dec'],url=url)
        try:
            raw = fetch(url)
            receipt.update(bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest())
            hd, table = unpack(raw)
            if int(hd['TICID']) != tic or int(hd['SECTOR']) != 105:
                raise ValueError('FITS identifiers differ from sample')
            metadata = dict(camera=int(hd['CAMERA']),ccd=int(hd['CCD']),tessmag=float(hd['TESSMAG']),
                            ra_obj=star['ra'],dec_obj=star['dec'])
            aux = {k:table[k].astype(float) for k in [*AUX_COLUMNS,'SAP_FLUX'] if k in table.dtype.names}
            c = Curve(table['TIME'].astype(float), table['DET_FLUX'].astype(float),
                      table['QUALITY']==0, aux, metadata)
            return receipt,c,table
        except (OSError, ValueError,KeyError) as exc:
            receipt.update(status='failed',detail=str(exc))
            return receipt,None,None

    audit, events, controls, receipts = [], [], [], []
    # Bounded batches keep large in-memory FITS from accumulating behind the scanner.
    with ThreadPoolExecutor(args.workers) as pool:
        for start in range(0,len(targets),args.workers):
            for receipt,curve,table in pool.map(load, enumerate(targets[start:start+args.workers],129+start)):
                if curve is not None:
                    tic = receipt['tic']
                    a,e,c = scan_one(curve,tic,105,8462852,8,require_shoulders=True)
                    a.update(tic=tic,sector=105,sha256=receipt['sha256'],**curve.meta)
                    audit.append(a)
                    events.extend(e)
                    controls.extend(c)
                    receipt.update(status='ok',detail='')
                    if e:
                        cols = {k:table[k] for k in table.dtype.names}
                        meta = {k.upper():v for k,v in curve.meta.items()}
                        np.savez_compressed(cache/f'{tic}.npz',**cols,**meta,TICID=tic,SECTOR=105)
                receipts.append(receipt)
            write_csv(args.out/'products.csv',receipts,fields)
            print(f'{len(receipts)} listed targets, {len(audit)} analyzed, {len(events)} triggers',flush=True)
    write_csv(args.out/'audit.csv',audit,list(audit[0]))
    write_csv(args.out/'triggers.csv',events,['tic','sector','camera','ccd','tessmag','ra_obj',
              'dec_obj','peak_btjd','start','stop','duration_h','ingress_fraction','depth','score','diagnostics'])
    # Fixed timestamp for byte-reproducible compression.
    with (args.out/'injections.csv.gz').open('wb') as out:
        with gzip.GzipFile(fileobj=out,mode='wb',mtime=0,filename='') as compressed:
            with io.TextIOWrapper(compressed, encoding='utf-8',newline='') as stream:
                writer=csv.DictWriter(stream,fieldnames=list(controls[0]))
                writer.writeheader()
                writer.writerows(controls)
    report = dict(target_list_sha256=digest,rank_start=129,rank_end=384,registered=256,
                  n_analyzed=len(audit),failures=len(receipts)-len(audit),all_trigger_count=len(events),
                  stars_triggered=len({e['tic'] for e in events}),seed=8462852,trials_per_cell=8,
                  wall_star_days=sum(a['span_d'] for a in audit),limits=rates(audit,controls,len(events)))
    (args.out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
