from common import D,np
import csv,json
rows=list(csv.DictReader((D/'swift-photometry.csv').open()))
rows=sorted((r for r in rows if float(r['target_3.5_exposure'])>20 and float(r['ref1_3.5_exposure'])>20 and 10<float(r['ref1_3.5_rate'])<70),key=lambda r:r['start_utc'])
out={}
for year in sorted({r['start_utc'][:4] for r in rows}):
    sel=[r for r in rows if r['start_utc'][:4]==year]
    out[year]=dict(exposures=len(sel),first=sel[0]['start_utc'],last=sel[-1]['date_end_utc'],
        within_half_depth=[int(sum(abs(float(x['delta_'+which+'_cycles']))<np.sqrt(2*np.log(2))*w for x in sel)) for which,w in [('primary',.024834876740958547),('secondary',.024377945718026622)]])
total=list((D/'swift_stamps').glob('*.json'))
receipts=[]
for p in sorted(total):
    r=json.loads(p.read_text()); ex=p.with_suffix('.expmeta')
    receipts.append({'obsid':r['obsid'],'url_sky':r['url'],'sha256_sky':r['sha256'],'n_whole_stamps':len(r['planes']),'url_exposure':r['url'].replace('_sk.img','_ex.img') if ex.exists() else None,'sha256_exposure':ex.read_text() if ex.exists() else None})
(D/'swift-receipts-replayed.json').write_text(json.dumps(receipts,indent=2))
(D/'coverage.json').write_text(json.dumps(dict(first=rows[0]['start_utc'],last=rows[-1]['date_end_utc'],n=len(rows),years=out,receipts=len(receipts)),indent=2))
