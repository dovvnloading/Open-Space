"""Join per-alert fits and immutable input receipts; keep failures in the denominator."""
import argparse
import collections
import csv
import gzip
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent


def rows(path,delimiter=','):
    with path.open() as f:
        yield from csv.DictReader(f,delimiter=delimiter)


def inputs():
    items={}
    ogle={r['event']:r for r in rows(HERE/'ogle/catalog.csv')}
    for r in rows(HERE/'ogle/photometry.csv'):
        n=r['event'];c=ogle[n]
        items['OGLE-'+n]=dict(survey='OGLE',year=n[:4],field=c['field'],input_sha256=r['sha256'],
                url=r['url'],catalog_te_d=c['tau_days'],time_system='HJD')
    for p in (HERE/'kmt').glob('catalog-*.tsv'):
        for r in rows(p,'\t'):
            n=r['Event'];items[n]=dict(survey='KMT',year=n[4:8],field=r['Field'],input_sha256='',
                url='',catalog_te_d=r['t_E'],time_system='HJD')
    for p in (HERE/'kmt').glob('pySIS-*.tsv'):
        for r in rows(p,'\t'):
            n=r['event']
            if r['sha256']:items[n].update(input_sha256=r['sha256'],url=r['url'])
    # Native receipts take precedence if an export was made while downloads ran.
    for r in rows(HERE/'other/photometry.csv'):
        n=r['event'];ss,y=n.split('-')[:2]
        items[n]=dict(survey=ss,year=y,field=r['field'],input_sha256=r['input_sha256'],
                url=r['phot_url'],catalog_te_d=r['page_te_d'],time_system='survey JD')
    for r in rows(HERE/'other/gaia.csv'):
        n=r['event'];items[n]=dict(survey='Gaia',year='20'+n[4:6],field='bulge rectangle',
                input_sha256=r['sha256'],url=r['url'],catalog_te_d='',time_system='JD TCB')
    return items


def write_gz(path,table):
    cols=list(dict.fromkeys(col for row in table for col in row))
    import io
    with path.open('wb') as raw, gzip.GzipFile(fileobj=raw,mode='wb',mtime=0,filename='') as zz, io.TextIOWrapper(zz,newline='') as f:
        w=csv.DictWriter(f,cols,lineterminator='\n');w.writeheader();w.writerows(table)


def read_results(folder):
    return {p.stem:json.loads(p.read_text()) for p in folder.glob('*.json')}


def run(cache):
    spec=inputs()
    for p in (cache/'staging/kmt').glob('*.receipt.json'):
        v=json.loads(p.read_text());n=p.name.split('.')[0]
        if n in spec:spec[n].update(input_sha256=v['sha256'],url=v['url'])
    for p in (cache/'staging/kmt/fallback').glob('*.json'):
        v=json.loads(p.read_text());n=p.stem
        if n in spec and v.get('sha256'):spec[n].update(input_sha256=v['sha256'],url=v['url'])
    broad=read_results(cache/'staging/screen');bound=read_results(cache/'staging/bound')
    cc=list(rows(cache/'staging/all-cadence.csv'))
    cad={r['event']:r for r in cc}
    gg=[];kk=[];flags=[]
    totals=collections.defaultdict(collections.Counter)
    for n,meta in sorted(spec.items()):
        g=broad.get(n)
        status=('ok' if g and g.get('status')=='ok' else 'fit failed' if g else
                'insufficient epochs' if meta['survey']=='Gaia' else 'pending')
        entry={'event':n,**meta,'fit_status':status}
        cnt=totals[(meta['survey'],meta['year'])]
        cnt['inventoried']+=1;cnt['hashed']+=bool(meta['input_sha256']);cnt['fitted']+=(status=='ok')
        if g:entry.update({k:v for k,v in g.items() if k not in ('event','status')})
        if status=='ok':
            cnt['usable_epochs']+=g['n'];cnt['ps_short']+=g['ps_te']<.5
            cnt['fs_short']+=g['fs_te']<.5;cnt['fs_large']+=g['fs_rho']>1
            for label,r in [('primary',g),('baseline-capped',bound.get(n))]:
                if not r or r.get('status')!='ok':continue
                reasons=[]
                if r['ps_te']<.5:reasons.append('PS-short')
                if r['fs_te']<.5:reasons.append('FS-short')
                if r['fs_rho']>1:reasons.append('large-source')
                if not reasons:continue
                flags.append(dict(event=n,pass_name=label,flags='|'.join(reasons),
                    ps_te=r['ps_te'],fs_te=r['fs_te'],fs_rho=r['fs_rho'],
                    fs_minus_ps_gain=round(r['fs_gain']-r['ps_gain'],2),
                    fs_minus_bin_gain=round(r['fs_gain']-r['bin_gain'],2),
                    n_before_best_0p5d=r['n_before_0.5d'],n_after_best_0p5d=r['n_after_0.5d']))
        gg.append(entry)
        if n in bound:kk.append(bound[n])
        if n in cad and not cad[n].get('error'):cnt['cadence_rows']+=1
    out=HERE/'tables';out.mkdir(exist_ok=True)
    write_gz(out/'grid-fits.csv.gz',gg)
    write_gz(out/'kmt-reference-flux-grid.csv.gz',kk)
    write_gz(out/'all-screening-flags.csv.gz',flags)
    write_gz(out/'cadence.csv.gz',sorted(cad.values(),key=lambda r:r['event']))
    summary=[dict(survey=k[0],year=k[1],**v) for k,v in sorted(totals.items())]
    (out/'coverage.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--cache',type=Path,default=Path.cwd())
    args=p.parse_args();print(json.dumps(run(args.cache)))
