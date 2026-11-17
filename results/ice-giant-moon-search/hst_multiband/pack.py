"""Collect replayable HST epoch results without retaining copies of the public FITS."""
import csv
import gzip
import hashlib
import io
from pathlib import Path
import tarfile
ROOT=Path(__file__).resolve().parent


def rows(path):
    return list(csv.DictReader(path.open())) if path.exists() else []


def write(path,records):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)


def main():
    epochs=[];measurements=[];seen=set()
    for planet in ['uranus','neptune']:
        dense=rows(ROOT/(planet+'_coverage.csv'))
        sparse={r['epoch']:r for r in rows(ROOT/(planet+'_sparse.csv'))}
        vetted={r['epoch']:r for r in rows(ROOT/(planet+'_vet.csv'))}
        for r in dense:
            tag=r['epoch'];d=ROOT/'epochs'/tag;s=sparse.get(tag)
            if r['status']!='screened' and s is None:raise ValueError('incomplete: '+tag)
            v=vetted.get(tag);branch='multi-exposure' if s is None else 'sparse'
            if s is None and int(r['dropout_hits']) and v is None:raise ValueError('unvetted '+tag)
            rec=dict(planet=planet,epoch=tag,filter=r['filter'],selection=branch,status=(s or r)['status'],
                     archive_exposures=r['archive_exposures'],registered=(s or r)['used_exposures'],
                     visits=(s or r)['used_visits'],grid=(s or r)['grid'],supported=(s or r)['supported'],
                     initial=(s or r)['screen_hits'],
                     after_dropout=r['dropout_hits'] if s is None else '',
                     remaining_after_static=(s or v or {'after_static_subtraction':0})['after_static_subtraction'])
            epochs.append(rec)
            reg={x['rootname']:x for x in rows(d/'registration.csv')}
            for x in rows(d/'inputs.csv'):
                name=x['file'];n=name.split('_')[0]
                if n in seen:raise ValueError('duplicate '+n)
                seen.add(n);rr=reg.get(n,{})
                measurements.append(dict(rootname=n,planet=planet,filter=r['filter'],epoch=tag,visit=n[4:6],status=rr.get('status',''),
                                        selection_status=rec['status'],bytes=x['bytes'],sha256=x['sha256'],uri='mast:HST/product/'+name))
    obs=rows(ROOT.parent/'manifests/observations.csv')
    expected={r['obs_id'] for r in obs if r['obs_collection']=='HST'}
    if seen!=expected:raise ValueError(f'{len(seen^expected)} HST files do not reconcile')
    write(ROOT/'epoch_index.csv',epochs);write(ROOT/'exposure_hashes.csv',measurements)
    followup=[]
    for path in sorted((ROOT/'epochs').glob('*/sparse_vet_filters.csv')):
        followup.extend(rows(path))
    if followup:write(ROOT/'sparse_followup.csv',followup)
    # Deterministic archive of the detector- and epoch-level diagnostic tables.
    with (ROOT/'epoch_details.tar.gz').open('wb') as raw:
        with gzip.GzipFile(filename='',fileobj=raw,mode='wb',mtime=0) as compressed:
            with tarfile.open(fileobj=compressed,mode='w') as bundle:
                for f in sorted((ROOT/'epochs').rglob('*')):
                    if not f.is_file():continue
                    content=f.read_bytes();info=tarfile.TarInfo(f.relative_to(ROOT).as_posix())
                    info.size=len(content);info.mtime=0;info.mode=0o644
                    bundle.addfile(info,io.BytesIO(content))


if __name__=='__main__':main()
