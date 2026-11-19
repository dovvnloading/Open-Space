#!/usr/bin/env python3
"""Fetch only selected DJA v4.4 individual .spec.fits (1D + possible 2D HDUs).

Must be preceded by publisher-MD5-verified select_targets.py output.
Never infers redshift, 2D validation, NIRCam match, or discovery from download.
"""
import argparse
import csv
import datetime as dt
import hashlib
import json
from pathlib import Path
import re
import subprocess
import urllib.request

HOST = 's3.amazonaws.com'
BASE = 'https://s3.amazonaws.com/msaexp-nirspec/extractions'
SAFE = re.compile(r'^[A-Za-z0-9_.-]+$')


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as fh:
        for chunk in iter(lambda: fh.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def safe(row):
    root, name = row['root'], row['file']
    if not root or not name or not SAFE.fullmatch(root) or not SAFE.fullmatch(name) or '..' in root or '..' in name or not name.endswith('.spec.fits'):
        raise ValueError(f'Unsafe path: {root!r}/{name!r}')
    return root, name, f'{BASE}/{root}/{name}'


def read_plan(selected, receipt_path):
    receipt = json.loads(receipt_path.read_text())
    observed = sha256(selected)
    if observed != receipt['selection_sha256']:
        raise ValueError('Selection file SHA-256 does not match selection receipt')
    if receipt['catalog_rows'] != 80367:
        raise ValueError('Full 80,367-row catalog freeze not evidenced')
    records = []
    seen = set()
    with selected.open(newline='') as fh:
        for row in csv.DictReader(fh):
            if row['targeted_selected'] != 'true':
                raise ValueError('Unexpected non-selected row in target list')
            root, name, url = safe(row)
            key = (root,name)
            if key not in seen:
                records.append({'root':root,'file':name,'url':url})
                seen.add(key)
    if receipt['selected_spectrum_rows'] < len(records) or receipt['selected_spectrum_rows'] <= 0:
        raise ValueError('Selected source count inconsistent with unique FITS paths')
    return receipt, records


def write_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj,indent=2)+'\n')


def preflight():
    try:
        proc = subprocess.run(['getent','ahostsv4',HOST],capture_output=True,text=True,timeout=4)
        if proc.returncode != 0:
            raise RuntimeError(f'DNS getent {HOST}: exit={proc.returncode}, stderr={proc.stderr.strip()!r}')
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError('DNS getent s3.amazonaws.com: timeout after 4s') from exc


def download_one(item, output, cap_mb=1024):
    target = output/item['root']/item['file']
    target.parent.mkdir(parents=True,exist_ok=True)
    if target.exists():
        return {'status':'EXISTING_BYTES_NOT_PUBLISHER_VERIFIED', 'sha256':sha256(target),'bytes':target.stat().st_size}
    temp=target.with_suffix('.fits.partial')
    total=0
    try:
        req=urllib.request.Request(item['url'],headers={'User-Agent':'Open-Space-DJA-targeted-research/1.0'})
        with urllib.request.urlopen(req,timeout=30) as response, temp.open('wb') as fh:
            if response.status != 200:
                raise ValueError(f'Unexpected HTTP {response.status}')
            while True:
                chunk=response.read(1024*1024)
                if not chunk: break
                total+=len(chunk)
                if total > cap_mb*1024*1024:
                    raise ValueError(f'Over per-product cap of {cap_mb} MiB')
                fh.write(chunk)
        if not total:
            raise ValueError('Empty body')
        temp.replace(target)
    finally:
        temp.unlink(missing_ok=True)
    return {'status':'RECEIVED_UNINDEPENDENTLY_VALIDATED_FITS_BYTES','sha256':sha256(target),'bytes':total}


def run(selected, receipt_file, out, plan_only=False):
    receipt, entries = read_plan(selected, receipt_file)
    out.mkdir(parents=True,exist_ok=True)
    plan={'made_at_utc':utc(), 'status':'PLAN_ONLY_NO_BYTES' if plan_only else 'PENDING_DOWNLOAD',
          'selection_sha256':receipt['selection_sha256'], 'selected_spectrum_rows':receipt['selected_spectrum_rows'],
          'distinct_fits_urls':len(entries), 'inputs':entries}
    write_json(out/'targeted_plan.json',plan)
    if plan_only:
        return {'status':'PLAN_ONLY_NO_BYTES','distinct_fits_urls':len(entries)}
    ledger=out/'download_receipts.jsonl'
    try:
        preflight()
    except Exception as exc:
        with ledger.open('a') as fh:
            fh.write(json.dumps({'at_utc':utc(),'stage':'preflight','status':'BLOCKED_DNS',
                                 'host':HOST, 'error':f'{type(exc).__name__}: {exc}',
                                 'downloaded':0})+'\n')
        return {'status':'BLOCKED_DNS','downloaded':0,'error':str(exc)}
    successes=0; failures=0
    with ledger.open('a') as fh:
        for item in entries:
            row={'at_utc':utc(),'url':item['url'],'root':item['root'],'file':item['file']}
            try:
                row.update(download_one(item,out/'fits'));successes+=1
            except Exception as exc:
                row.update({'status':'DOWNLOAD_FAILED','error':f'{type(exc).__name__}: {exc}'})
                failures+=1
            fh.write(json.dumps(row)+'\n'); fh.flush()
    return {'status':'DOWNLOADS_COMPLETED_NOT_VALIDATED' if not failures else 'PARTIAL_FAILURE',
            'downloaded_or_cached':successes,'failed':failures}


def self_test():
    assert safe({'root':'jades-v4','file':'jades-v4_prism-clear_1180_1.spec.fits'})[2].endswith('.spec.fits')
    for root in ('../evil','evil/foo','..'):
        try: safe({'root':root,'file':'x.spec.fits'})
        except ValueError: pass
        else: raise AssertionError('path escaped allowlist')
    print('PASS scoped per-source URL and traversal rejection (synthetic only)')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--selected',type=Path);p.add_argument('--selection-receipt',type=Path)
    p.add_argument('--out',type=Path);p.add_argument('--plan-only',action='store_true');p.add_argument('--self-test',action='store_true')
    args=p.parse_args()
    if args.self_test:self_test()
    else:
        if not args.selected or not args.selection_receipt or not args.out:
            p.error('--selected --selection-receipt and --out are required')
        print(json.dumps(run(args.selected,args.selection_receipt,args.out,args.plan_only),indent=2))
