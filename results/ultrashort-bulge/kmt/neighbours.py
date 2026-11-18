"""Query Gaia DR3 and VIRAC2 in 1.5-arcsec cones centered on short-list survey coordinates."""
import argparse
import concurrent.futures
import csv
import io
import json
import math
import urllib.parse
from pathlib import Path

from freeze import HERE, get

GAIA = 'https://gea.esac.esa.int/tap-server/tap/sync'
VIRAC = 'https://vizier.cfa.harvard.edu/viz-bin/asu-tsv'
RADIUS = 1.5
EXTRAS = {f'KMT-2025-BLG-{i}' for i in ('0220','1213','1495','1587','2047','2688','2874')} | {'KMT-2024-BLG-'+s for s in ('0251','0909','3326','3418','3429')}


def angle(s, ra=False):
    a,b,c=map(float,s.lstrip('-+').split(':'))
    v=a+b/60+c/3600
    return v*(15 if ra else -1 if s.startswith('-') else 1)


def query(event, database, cache):
    name,ra,dec = event['event'],angle(event['ra'],True),angle(event['dec'])
    outfile=cache/(name+'-'+database+'.json')
    if outfile.exists() and json.loads(outfile.read_text()).get('status') == 200:
        return
    if database=='Gaia':
        adql = (f'SELECT source_id, ra, dec, phot_g_mean_mag, bp_rp, phot_variable_flag, ruwe, pmra, pmdec '
                f"FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec),"
                f"CIRCLE('ICRS',{ra:.9f},{dec:.9f},{RADIUS/3600:.10f}))")
        url=GAIA+'?'+urllib.parse.urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=adql))
    else:
        cols='srcid,RAJ2000,DEJ2000,pmRA,pmDE,Jmag,Hmag,Ksmag,e_Ksmag,o_Ksmag,MedKsmag,SJKsmag,etaKs,KsEp1,KsEpl'
        url=VIRAC+'?'+urllib.parse.urlencode({'-source':'II/387/virac2','-c':f'{ra:.9f} {dec:.9f}',
                                               '-c.rs':str(RADIUS),'-out.max':'100','-out':cols})
    try:
        raw,rec=get(url)
        text=raw.decode()
        if database=='Gaia':
            if not text.startswith('source_id,'):
                raise ValueError('missing Gaia CSV header')
            rows=list(csv.DictReader(io.StringIO(text)))
        else:
            if '#INFO\tError=' in text or '#RESOURCE=yCat_2387' not in text:
                raise ValueError('unexpected VizieR response')
            lines=[s for s in text.splitlines() if s.strip() and not s.startswith('#')]
            if not lines:
                rows=[]
            elif lines[0].startswith('srcid\t'):
                cols=lines[0].split('\t')
                rows=[dict(zip(cols,[v.strip() for v in l.split('\t')])) for l in lines[3:]]
            else:
                raise ValueError('unexpected VizieR columns')
        for r in rows:
            rra=float(r.get('ra',r.get('RAJ2000')))
            rdec=float(r.get('dec',r.get('DEJ2000')))
            r['sep_arcsec']=round(math.hypot((rra-ra)*math.cos(math.radians(dec)),rdec-dec)*3600,4)
        rec.update(event=name,database=database,radius_arcsec=RADIUS,matches=sorted(rows,key=lambda x:x['sep_arcsec']))
        (cache/(name+'-'+database+'.raw')).write_bytes(raw)
    except (OSError,ValueError) as exc:
        rec=dict(event=name,database=database,url=url,status='failed',error=str(exc))
    outfile.write_text(json.dumps(rec)+'\n')


def summarize(cache):
    receipt_fields=['event','database','radius_arcsec','status','utc','size','sha256','url','error']
    records=[json.loads(p.read_text()) for p in sorted(cache.glob('*.json'))]
    with (HERE/'neighbour-receipts.tsv').open('w') as h:
        w=csv.DictWriter(h,fieldnames=receipt_fields,extrasaction='ignore',delimiter='\t',lineterminator='\n')
        w.writeheader()
        for r in records:
            w.writerow(r)
    for db in ['Gaia','VIRAC2']:
        keep=[]
        for rec in records:
            if rec['database']==db and rec['status']==200:
                for r in rec['matches']:
                    keep.append(dict(event=rec['event'],**r))
        if keep:
            with (HERE/(db+'-neighbours.tsv')).open('w') as h:
                w=csv.DictWriter(h,list(dict.fromkeys(k for r in keep for k in r)),delimiter='\t',lineterminator='\n')
                w.writeheader()
                w.writerows(keep)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cache',type=Path,required=True)
    p.add_argument('--event',action='append',default=[])
    p.add_argument('--summarize-only',action='store_true')
    args=p.parse_args()
    cache=args.cache/'neighbours'
    cache.mkdir(exist_ok=True,parents=True)
    if not args.summarize_only:
        rows=[r for r in map(json.loads,(args.cache/'index.jsonl').read_text().splitlines())
              if (r['te'] is not None and 0<r['te']<=2) or r['event'] in EXTRAS or r['event'] in args.event]
        rows.sort(key=lambda r:(r['event'] not in args.event,r['te'] or 1e6))
        with concurrent.futures.ThreadPoolExecutor(4) as ex:
            list(ex.map(lambda v:query(*v,cache),((r,d) for r in rows for d in ['VIRAC2','Gaia'])))
    summarize(cache)


if __name__=='__main__':
    main()
