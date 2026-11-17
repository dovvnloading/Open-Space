"""Freeze comparison stars, three-arcsecond catalogs and 2016–25 field curves."""
import argparse
import csv
import io
import json
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode

from support import HERE, KMT, PORTAL, read
from neighbours import angle

EVENTS = ['KMT-2024-BLG-0824','KMT-2024-BLG-3278','KMT-2025-BLG-1495','KMT-2025-BLG-0220']
CFA = 'https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?'
GAIA = 'https://gea.esac.esa.int/tap-server/tap/sync?'
CATS = [('virac','II/387/virac2','srcid,RAJ2000,DEJ2000,Jmag,Hmag,Ksmag,e_Ksmag,o_Ksmag,MedKsmag,SJKsmag,etaKs,KsEp1,KsEpl'),
        ('vvv','II/348/vvv2','srcid,raJ2000,deJ2000,Jmag3,Hmag3,Ksmag3,Jperrbits,Ksperrbits'),
        ('vsx','B/vsx/vsx','Name,RAJ2000,DEJ2000,Type,Period')]


def catalog(event):
    path=KMT/('catalog-'+event.split('-')[1]+'.tsv')
    with path.open() as f:
        return next(r for r in csv.DictReader(f,delimiter='\t') if r['Event']==event)


def fetch_catalog(args):
    cache, name, url=args
    raw,receipt=read(cache,name,url)
    if name=='gaia':
        if not raw.startswith(b'source_id,'):
            raise ValueError(receipt)
        rows=list(csv.DictReader(io.StringIO(raw.decode())))
    else:
        lines=raw.decode().splitlines()
        if not any(s.startswith('#RESOURCE=yCat') for s in lines):
            (cache/(name+'.receipt')).unlink()
            raise ValueError(f'catalog service error for {name}')
        lines=[s for s in lines if s.strip() and not s.startswith('#')]
        rows=[]
        if lines:
            rows=[dict(zip(lines[0].split('\t'),[v.strip() for v in s.split('\t')])) for s in lines[3:]]
    return dict(name=name,receipt=receipt,rows=rows)


def freeze_event(event,cache):
    row=catalog(event);ra,dec=angle(row['RA'],True),angle(row['Dec'])
    c=cache/event;c.mkdir(parents=True,exist_ok=True)
    url=PORTAL+urlencode(dict(search_type='radec',ra=f'{ra:.9f}',dec=f'{dec:.9f}',ra_type='decimal_degree',arcsec=15))
    buf,rec=read(c,'cone',url)
    markers=json.loads(re.search(rb'const kmtnetSearchMarkers = (\[.*?\]);',buf)[1])
    field='BLG'+row['Field'][3:]
    matches=[m for m in markers if m['lcid'].startswith(field+'.')]
    target=next((m for m in matches if m['lcid'].endswith('.'+row['Star ID'])),None)
    if target is None:
        target=min(matches or markers,key=lambda x:x['dist'])
    # Compare the two field centers at the same sky position when available and
    # two separate stars in the primary field (nearest and widest listed).
    picks=[target]
    samepos=[m for m in markers if m['lcid']!=target['lcid'] and m['ra']==target['ra'] and m['dec']==target['dec']]
    others=[m for m in matches if m['lcid']!=target['lcid']]
    picks.extend(samepos[:1])
    if others:
        picks.extend([others[0],others[-1]])
    picks=list({m['lcid']:m for m in picks}.values())
    def get_star(m):
        lcid=m['lcid'];u=PORTAL+urlencode(dict(action='download_lc',lcid=lcid))
        _,r=read(c,lcid,u)
        return dict(**m,receipt=r)
    with ThreadPoolExecutor(2) as pool:
        meta=list(pool.map(get_star,picks))
    urls=[]
    qry=f"""SELECT source_id,ra,dec,phot_g_mean_mag,bp_rp,phot_variable_flag,ruwe,pmra,pmdec FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra:.9f},{dec:.9f},{3/3600:.10f}))"""
    urls.append((c,'gaia',GAIA+urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=qry))))
    for name,cat,cols in CATS:
        u=CFA+urlencode({'-source':cat,'-c':f'{ra:.9f} {dec:.9f}','-c.rs':3,'-out.add':'_r','-sort':'_r','-out.max':50,'-out':cols})
        urls.append((c,name,u))
    cats=[]
    for params in urls:
        try: cats.append(fetch_catalog(params))
        except (OSError,ValueError) as exc: cats.append(dict(name=params[1],url=params[2],error=str(exc)))
    return dict(event=event,ra_deg=ra,dec_deg=dec,cone=rec,reference_stars=meta,catalogs=cats)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cache',type=Path,required=True)
    args=p.parse_args()
    with ThreadPoolExecutor(2) as pool:
        events=list(pool.map(lambda event: freeze_event(event,args.cache),EVENTS))
    (HERE/'sightline-receipts.json').write_text(json.dumps(events,indent=2)+'\n')


if __name__=='__main__':
    main()
