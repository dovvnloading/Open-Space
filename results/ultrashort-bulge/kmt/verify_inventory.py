"""Reconcile released event IDs, raw or safely evicted bundles, fits and detector checks offline."""
import argparse
import csv
import gzip
import hashlib
import json
import tarfile
from collections import Counter
from pathlib import Path

from freeze import HERE


def sha256(path):
    with path.open('rb') as h:
        return hashlib.file_digest(h,'sha256').hexdigest()


def verify(cache,stage):
    rows=[json.loads(line) for line in (cache/'index.jsonl').read_text().splitlines()]
    catalog={r['event'] for r in rows}
    if len(catalog)!=len(rows):raise ValueError('Duplicate released event IDs')
    if len(catalog)!=6789 or any(n[4:8] not in {'2024','2025'} for n in catalog):
        raise ValueError('Unexpected frozen event list')
    receipts=json.loads((HERE/'index-receipts.json').read_text())
    for year in (2024,2025):
        rec=next(r for r in receipts['queries'] if r['url'].endswith(f'/{year}/'))
        raw=gzip.decompress((HERE/f'index-{year}.html.gz').read_bytes())
        if len(raw)!=rec['size'] or hashlib.sha256(raw).hexdigest()!=rec['sha256']:
            raise ValueError(f'Frozen {year} annual catalog checksum mismatch')
    seasons={y:Counter() for y in ('2024','2025')};problems=[];fitted={};expected={}
    for r in rows:
        n=r['event'];c=seasons[n[4:8]];c['listed']+=1
        receipt=cache/(n+'.receipt.json');err=cache/(n+'.error.json');alt=cache/'fallback'/(n+'.json')
        path=cache/(n+'.tar.gz');marker=cache/(n+'.evicted');mode='pySIS'
        if receipt.is_file():
            rec=json.loads(receipt.read_text());members=rec['curves']
        elif alt.exists() and err.exists() and json.loads(err.read_text())['status']==404 and json.loads(alt.read_text())['status']=='ok':
            result=json.loads(alt.read_text());rec=result['curve'];members=result['members'];mode='DIAPL';path=alt.with_name(n+'.tar.gz')
        else:
            problems.append([n,'missing successful acquisition']);continue
        if path.exists():
            if path.stat().st_size!=rec['size'] or sha256(path)!=rec['sha256']:
                problems.append([n,'stored bytes mismatch']);continue
            with tarfile.open(path) as h:
                suffix='.pysis' if mode=='pySIS' else '.diapl'
                exact=[m.name for m in h if m.isfile() and m.name.endswith(suffix)]
            if members!=exact or not exact:
                problems.append([n,'archive member mismatch']);continue
            c['retained_rehashed']+=1
        elif marker.exists() and mode=='pySIS' and json.loads(marker.read_text())['sha256']==rec['sha256']:
            c['eviction_hash_receipt_verified']+=1
        else:
            problems.append([n,'missing stored bytes and matching eviction receipt']);continue
        c[mode+'_bundles']+=1;c[mode+'_I']+=sum('_I.' in s for s in members);c[mode+'_V']+=sum('_V.' in s for s in members)
        c['downloaded_bytes']+=rec['size']
        for key in ('screen','bound'):
            g=stage/key/(n+'.json')
            if not g.is_file():problems.append([n,key+' grid pending']);continue
            model=json.loads(g.read_text())
            if model.get('status')!='ok':problems.append([n,key+' grid '+str(model.get('status'))]);continue
            c[key+'_grid']+=1
            if key=='bound':
                if model['ps_te']<.5:
                    expected[(n,'ps')]=(model['ps_te'],0.,model['ps_t0']-2400000.)
                if model['fs_te']<.5 or model['fs_rho']>1:
                    expected[(n,'fs')]=(model['fs_te'],model['fs_rho'],model['fs_t0']-2400000.)
                fitted[n]=model['n']
    with (stage/'all-cadence.csv').open() as h:
        timed={r['event'] for r in csv.DictReader(h) if not r.get('error')}
    for n in sorted(catalog):
        if n in timed:seasons[n[4:8]]['with_sampling']+=1
        else:problems.append([n,'sampling pending'])
    with (HERE/'short-lead-audit.tsv').open() as h:
        actual={(r['event'],r['model']):(float(r['te_d']),float(r['rho']),float(r['t0']))
                for r in csv.DictReader(h,delimiter='\t') if r['status']=='ok'}
    for n,m in sorted(expected.keys() ^ actual.keys()):
        problems.append([n,m+' detector audit differs from capped grid'])
    for key in expected.keys() & actual.keys():
        n,m=key
        if not all(abs(a-b)<1e-7 for a,b in zip(actual[key],expected[key])):
            problems.append([n,m+' detector audit has stale geometry'])
        else:seasons[n[4:8]]['audited_'+m]+=1
    return dict(frozen_public_ids=len(catalog),by_year=seasons,distinct_detector_audit_events=len({n for n,_ in actual}),
                fitted_I_band_observations=sum(fitted.values()),problems=problems)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cache',type=Path,required=True);p.add_argument('--stage',type=Path,default=Path('staging'))
    p.add_argument('--out',type=Path);args=p.parse_args();result=verify(args.cache,args.stage)
    message=json.dumps(result,indent=2)+'\n'
    if args.out:args.out.write_text(message)
    else:print(message,end='')
    if result['problems']:raise SystemExit(1)


if __name__=='__main__':main()
