"""Vet priority continuous results using independent alerts, previous epochs and colors."""
import csv
import json
import re
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode

import numpy as np
from lxml import html

from fetch import OUT,SHA,write_csv,get
from crossmatch import coords,extern
from known import frozen
from flag_alternatives import candidate_ids
from flag_refits import WORK
from run_grid import inputs,read_curve


def catalog():
    result=frozen()
    results={x[0]:x for x in result}
    source_path=Path('staging/kmt/2026.html')
    doc=html.fromstring(source_path.read_bytes())
    trs=doc.xpath('//tr');cols=[td.text_content().strip() for td in trs[0].xpath('./td')]
    for tr in trs[1:]:
        r=dict(zip(cols,[e.text_content().strip() for e in tr.xpath('./td')]))
        if r['RA']!='-' and r['Dec']!='-':
            results[r['Event']]=(r['Event'],*coords(r['RA'],r['Dec']),None)
    # Ancillary preceding season for 2024 MOA/PRIME; never included in discovery denominator.
    ogle=Path('.ogle-cache/raw/2024/lenses.par')
    if ogle.exists():
        for line in ogle.read_text().splitlines():
            c=line.split()
            if len(c)>5 and re.match(r'^2024-BLG-\d+$',c[0]):
                results['OGLE-'+c[0]]=('OGLE-'+c[0],*coords(c[3],c[4]),float(c[5]))
    return results


def add_periods(a):
    times,flux,err,g,diff=read_curve(a)
    use=np.isfinite(times+flux+err)&(times>2450000)&(times<2470000)&(err>0)
    t,f,e=times[use],flux[use],err[use];o=np.argsort(t);t,f,e=t[o],f[o],e[o]
    year=(1970+(t-2440587.5)/365.2425).astype(int)
    summary=[]
    for s in np.unique(year):
        mask=year==s
        if mask.sum()<3:continue
        tt,ff,ee=t[mask],f[mask],e[mask]
        mid=np.median(ff);mad=1.4826*np.median(abs(ff-mid))
        z=(ff-mid)/np.hypot(ee,mad)
        triple=np.minimum.reduce([z[:-2],z[1:-1],z[2:]])
        triple[tt[2:]-tt[:-2]>1]=-np.inf
        if len(triple) and np.max(triple)>5:summary.append(str(s))
    return '|'.join(summary)


if __name__=='__main__':
    ids=candidate_ids()
    cats=catalog(); jobs={j[0]:j for j in inputs(Path.cwd()) if j[0] in ids}
    phot={r['event']:r for r in csv.DictReader((OUT/'photometry.csv').open())}
    official=list(csv.DictReader((OUT/'official_links.csv').open()))
    alerts=list(cats.values())
    task=[(n,*cats[n][1:3],survey) for n in sorted(ids) for survey in ('Gaia DR3','VVV DR2')]
    with ThreadPoolExecutor(7) as work: results=list(work.map(extern,task))
    receipts=[v for v,rr in results];neighbors=[r for _,rr in results for r in rr]
    write_csv(OUT/'flag_catalog_queries.csv',receipts);write_csv(OUT/'flag_catalog_neighbours.csv',neighbors)
    ga=[r for r in neighbors if r['catalog']=='Gaia DR3' and r['phot_variable_flag']=='VARIABLE']
    gids=sorted({r['source_id'] for r in ga})
    if gids:
        q='SELECT source_id,best_class_name,best_class_score FROM gaiadr3.vari_classifier_result WHERE source_id IN ('+','.join(gids)+')'
        u='https://gea.esac.esa.int/tap-server/tap/sync?'+urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=q))
        cache=Path('.survey-cache/catalogs')/('priority-classes-'+SHA(u.encode())[:8])
        if not cache.exists():cache.write_bytes(get(u))
        raw=cache.read_bytes()
        classified={r['source_id']:r for r in csv.DictReader(raw.decode().splitlines())}
        write_csv(OUT/'flag_gaia_classes.csv',[dict(event=r['event'],**classified[r['source_id']]) for r in ga if r['source_id'] in classified])
        receipts.append(dict(event='classified variables',catalog='Gaia DR3 classifier',url=u,bytes=len(raw),sha256=SHA(raw),status='ok',neighbours=len(classified)))
        write_csv(OUT/'flag_catalog_queries.csv',receipts)
    else:classified={}
    records=[]
    for n in sorted(ids):
        _,ra,dec,epoch=cats[n];year=int(n.split('-')[1])
        aliases=[]
        for src,x,y,tt in alerts:
            if n==src or not np.isfinite(x+y):continue
            s=3600*np.hypot((x-ra)*np.cos(np.deg2rad(dec)),y-dec)
            if s>3:continue
            yr=2000+int(src[4:6]) if src.startswith('Gaia') else int(src.split('-')[1])
            dt=(float(tt)-float(epoch)) if epoch and tt else None
            aliases.append(f'{src} [{s:.2f} arcsec; '+(f'{dt:+.2f} d' if dt is not None else 'epoch omitted')+(']' if yr==year else '; other season]'))
        off=set(add_periods(jobs[n]).split('|'))-{str(year),''}
        rr=phot.get(n,{})
        cp=[r for r in official if r['id_1']==n or r['id_2']==n]
        linked=sorted(set(r['id_2'] if r['id_1']==n else r['id_1'] for r in cp))
        records.append(dict(event=n,ra_deg=round(ra,7),dec_deg=round(dec,7),assessment=rr.get('assessment',''),
            position_matches_3arcsec=' | '.join(aliases),project_associations='|'.join(linked),
            offseason_3consecutive_gt5_in_1d='|'.join(sorted(off)),
            gaia_variable_neighbours='|'.join(r['source_id']+':'+classified.get(r['source_id'],{}).get('best_class_name','unclassified') for r in ga if r['event']==n)))
    write_csv(OUT/'flag_vetting.csv',records)
