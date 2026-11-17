"""Join per-alert fits and immutable input receipts; keep failures in the denominator."""
import argparse
import collections
import csv
import gzip
import json
import hashlib
from datetime import date
from pathlib import Path
from replay_guard import external,fresh_for_inputs,sha,verify_receipt,frozen_parameters,parameters_digest

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
                url=r['url'],catalog_te_d=c['tau_days'],time_system='HJD',native_usable=r['n'])
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
                url=r['phot_url'],catalog_te_d=r['page_te_d'],time_system='survey JD',native_usable=r['usable'])
    for r in rows(HERE/'other/gaia.csv'):
        n=r['event'];items[n]=dict(survey='Gaia',year='20'+n[4:6],field='bulge rectangle',
                input_sha256=r['sha256'],url=r['url'],catalog_te_d='',time_system='JD TCB',native_usable=r['usable'])
    return items


def write_gz(path,table):
    cols=list(dict.fromkeys(col for row in table for col in row))
    import io
    with path.open('wb') as raw, gzip.GzipFile(fileobj=raw,mode='wb',mtime=0,filename='') as zz, io.TextIOWrapper(zz,newline='') as f:
        w=csv.DictWriter(f,cols,lineterminator='\n');w.writeheader();w.writerows(table)


def read_results(folder,spec,cap):
    folder=external(folder);meta=json.loads((folder/'campaign.json').read_bytes())
    if meta['cap']!=cap:raise ValueError('wrong checkpoint cap')
    timing=frozen_parameters(cap)
    if meta.get('parameters_sha256')!=parameters_digest(timing):raise ValueError('campaign timing map mismatch')
    for name,digest in meta['code_sha256'].items():
        if sha(name)!=digest:raise ValueError('checkpoint code changed')
    result={}
    for path in folder.glob('*.json'):
        if path.name=='campaign.json':continue
        row=json.loads(path.read_bytes());ident=row.pop('_identity',None);name=path.stem
        if not ident or ident['event']!=name or ident['cap']!=cap or ident['code_sha256']!=meta['code_sha256'] or ident['runtime']!=meta['runtime'] or name not in spec or ident['input_sha256']!=spec[name]['input_sha256'] or ident.get('parameters')!=timing.get(name):raise ValueError('checkpoint identity mismatch: '+name)
        result[name]=row
    expected={n for n in spec if not cap or n.startswith('KMT-')}
    if set(result)!=expected:raise ValueError('incomplete checkpoint membership')
    return result


def run(cache,screen,bound,cadence,out):
    cache=external(cache);screen=external(screen);bound=external(bound);cadence=external(cadence)
    verify_receipt('tables/grid-fits.csv.gz')
    for name in ('ogle/catalog.csv','ogle/photometry.csv','kmt/catalog-2024.tsv','kmt/catalog-2025.tsv','kmt/pySIS-2024.tsv','kmt/pySIS-2025.tsv','other/photometry.csv','other/gaia.csv'):verify_receipt(name)
    spec=inputs()
    broad=read_results(screen,spec,False);bound=read_results(bound,spec,True)
    verify_receipt('tables/cadence.csv.gz')
    frozen_cadence=gzip.decompress((HERE/'tables/cadence.csv.gz').read_bytes())
    if cadence.read_bytes()!=frozen_cadence:raise ValueError('cadence bytes differ from frozen table')
    cc=list(rows(cadence))
    cad={r['event']:r for r in cc}
    # Two Gaia CSVs have only two numeric transits and never enter the flux
    # reducer. Their two-exposures-per-wing exposure is nonetheless exactly zero.
    for n,meta in spec.items():
        if meta['survey']!='Gaia' or int(meta['native_usable'])>=3 or n in cad:continue
        raw=(cache/'.survey-cache/gaiaraw'/(n+'.csv')).read_bytes()
        assert hashlib.sha256(raw).hexdigest()==meta['input_sha256']
        times=[]
        for record in csv.reader(raw.decode().splitlines()[2:]):
            if record and record[2].strip().replace('.','',1).isdigit():times.append(float(record[1]))
        assert len(times)==int(meta['native_usable'])
        s=2451544.5+(date(int(meta['year']),3,1)-date(2000,1,1)).days
        cad[n]=dict(event=n,season_d=245.,distinct_epochs_in_season=len({t for t in times if s<=t<s+245}),
                    **{'days_'+w:0. for w in ('0.03','0.1','0.25','0.5','1.5')},error='')
    gg=[];kk=[];flags=[]
    totals=collections.defaultdict(collections.Counter)
    for n,meta in sorted(spec.items()):
        g=broad.get(n)
        status=('ok' if g and g.get('status')=='ok' else
                'insufficient epochs' if meta['survey']=='Gaia' and int(meta['native_usable'])<6 else
                'fit failed' if g else 'pending')
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
    out=fresh_for_inputs(out,cache,screen,bound,cadence)
    write_gz(out/'grid-fits.csv.gz',gg)
    write_gz(out/'kmt-reference-flux-grid.csv.gz',kk)
    write_gz(out/'all-screening-flags.csv.gz',flags)
    write_gz(out/'cadence.csv.gz',sorted(cad.values(),key=lambda r:r['event']))
    summary=[dict(survey=k[0],year=k[1],**v) for k,v in sorted(totals.items())]
    (out/'coverage.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--cache',type=Path,required=True);p.add_argument('--screen',type=Path,required=True);p.add_argument('--bound',type=Path,required=True);p.add_argument('--cadence',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    args=p.parse_args();print(json.dumps(run(args.cache,args.screen,args.bound,args.cadence,args.out)))
