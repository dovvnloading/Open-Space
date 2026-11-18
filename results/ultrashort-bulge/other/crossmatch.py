"""Spatial neighbours in the frozen alert lists and external catalogs for short leads.

A neighbour is a match in position only; check the epochs, light curve, and crowding.
"""
import argparse
import csv
import io
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from scipy.spatial import cKDTree
import numpy as np

from fetch import OUT, SHA, get, write_csv

CACHE=Path('.survey-cache')


def coords(ra,dec):
    def degree(v):
        bits=v.strip().replace('+','').split(':')
        return (-1 if v.strip().startswith('-') else 1)*(abs(float(bits[0]))+float(bits[1])/60+float(bits[2])/3600.)
    return 15*degree(ra),degree(dec)


def catalog():
    out=[]
    for r in csv.DictReader((OUT/'index.csv').open()):
        ra,dec=coords(r['ra'],r['dec'])
        meta=CACHE/'records'/(r['event']+'.json')
        epoch=json.loads(meta.read_text()).get('page_t0_jd') if meta.exists() else None
        out.append((r['event'],ra,dec,epoch))
    for p in sorted((OUT.parent/'kmt').glob('catalog-*.tsv')):
        for r in csv.DictReader(p.open(),delimiter='\t'):
            try: out.append((r['Event'],*coords(r['RA'],r['Dec']),float(r['t_0'])+2400000))
            except (ValueError,TypeError):continue
    for r in csv.DictReader((OUT.parent/'ogle/catalog.csv').open()):
        out.append(('OGLE-'+r['event'],*coords(r['ra_hms'],r['dec_dms']),r['t0_hjd']))
    for r in csv.DictReader((OUT/'gaia.csv').open()):
        meta=CACHE/'records'/(r['event']+'.json')
        epoch=json.loads(meta.read_text()).get('page_t0_jd') if meta.exists() else None
        out.append((r['event'],float(r['ra_deg']),float(r['dec_deg']),epoch))
    return out


def pairs(items):
    a,d=np.deg2rad([[r[1],r[2]] for r in items]).T
    xyz=np.column_stack([np.cos(d)*np.cos(a),np.cos(d)*np.sin(a),np.sin(d)])
    neighbours=cKDTree(xyz).query_pairs(np.deg2rad(2./3600.))
    out=[]
    for i,j in sorted(neighbours):
        r,s=items[i],items[j]
        if not r[0].startswith(('PRIME','MOA','Gaia')) and not s[0].startswith(('PRIME','MOA','Gaia')):
            continue
        theta=np.rad2deg(2*np.arcsin(np.linalg.norm(xyz[i]-xyz[j])/2))*3600
        dt=float(s[3])-float(r[3]) if r[3] and s[3] else None
        out.append(dict(id_a=r[0],id_b=s[0],separation_arcsec=round(theta,6),delta_peak_or_trigger_d=round(dt,6) if dt is not None else ''))
    write_csv(OUT/'matches.csv',out)


def extern(job):
    event,ra,dec,survey=job
    if survey=='Gaia DR3':
        q=(f"SELECT TOP 100 source_id,ra,dec,phot_g_mean_mag,bp_rp,phot_variable_flag,ruwe,pmra,pmdec "
           f"FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra:.8f},{dec:.8f},0.000555556))")
        url='https://gea.esac.esa.int/tap-server/tap/sync?'+urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=q))
    else:
        params={'-source':'II/348/vvv2','-c':f'{ra:.8f} {dec:.8f}','-c.rs':'2','-out.max':'10',
                '-out.add':'_r','-out':'srcid,Jmag3,Hmag3,Ksmag3','-sort':'_r'}
        url='https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?'+urlencode(params)
    receipt=dict(event=event,catalog=survey,url=url)
    try:
        path=CACHE/'catalogs'/(event+('-gaia.csv' if survey=='Gaia DR3' else '-vvvsmall.tsv'))
        if not path.exists():
            path.write_bytes(get(url))
        raw=path.read_bytes()
        if survey=='Gaia DR3':
            if not raw.startswith(b'source_id,'):
                raise ValueError('unexpected TAP return')
            rows=list(csv.DictReader(io.StringIO(raw.decode())))
            for r in rows:
                sep=3600*np.hypot((float(r['ra'])-ra)*np.cos(np.deg2rad(dec)),float(r['dec'])-dec)
                r['sep_arcsec']=round(sep,5)
        else:
            if b'#INFO\tservice_protocol=ASU' not in raw or b'#ERROR' in raw:
                raise ValueError('unexpected ASU return')
            lines=[line for line in raw.decode().splitlines() if line and not line.startswith('#')]
            if not lines and b'#Column' in raw:
                raise ValueError('response missing data header; completeness unverified')
            rows=list(csv.DictReader([*lines[:1],*lines[3:]],delimiter='\t')) if lines else []
            for r in rows:
                r.update(sep_arcsec=float(r.pop('_r')))
        receipt.update(status='ok',bytes=len(raw),sha256=SHA(raw),neighbours=len(rows))
        return receipt,[dict(event=event,catalog=survey,**r) for r in sorted(rows,key=lambda x:x['sep_arcsec'])]
    except Exception as e:
        receipt.update(status='error',error=str(e))
        return receipt,[]


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--add',nargs='*',default=[])
    args=p.parse_args()
    items=catalog()
    pairs(items)
    (CACHE/'catalogs').mkdir(exist_ok=True)
    names=set(args.add)
    for r in csv.DictReader((OUT/'photometry.csv').open()):
        try:
            if 0<float(r['page_te_d'])<.5:
                names.add(r['event'])
        except (ValueError,KeyError):pass
    for r in csv.DictReader((OUT.parent/'ogle/catalog.csv').open()):
        if 0<float(r['tau_days'])<.5:
            names.add('OGLE-'+r['event'])
    tasks=[(name,ra,dec,source) for name,ra,dec,_ in items if name in names for source in ['Gaia DR3','VVV DR2']]
    with ThreadPoolExecutor(8) as pool: data=list(pool.map(extern,tasks))
    write_csv(OUT/'catalog_query_receipts.csv',[r for r,_ in data])
    write_csv(OUT/'catalog_neighbours.csv',[r for _,rr in data for r in rr])
