"""Archive typed SIMBAD, 2MASS, AllWISE and GALEX responses for followed-up stars."""
import csv
import hashlib
import json
import math
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE=Path(__file__).parent
IDS={'3378588057203660160','3431326755205579264','1830677037678583680'}
SIMBAD='https://simbad.cfa.harvard.edu/simbad/sim-tap/sync'
IRSA='https://irsa.ipac.caltech.edu/TAP/sync'
MAST='https://mast.stsci.edu/api/v0/invoke'
MAM='https://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt'


def tap(endpoint,adql):
    return endpoint+'?'+urllib.parse.urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=adql))


def position(r,year):
    dt=year-float(r['ref_epoch'])
    dec=float(r['dec'])
    return (float(r['ra'])+dt*float(r['pmra'])/(3600000*math.cos(math.radians(dec))),
            dec+dt*float(r['pmdec'])/3600000)


def jobs(r):
    sid=r['source_id']
    sub=(f" FROM ident AS i JOIN basic AS b ON i.oidref=b.oid "
         f"JOIN allfluxes AS v ON b.oid=v.oidref WHERE i.id='Gaia DR3 {sid}'")
    yield sid+'_optical.csv',tap(SIMBAD,"SELECT b.main_id,b.otype,b.sp_type,b.sp_bibcode,v.B,v.V,v.G"+sub)
    for key,year,table,fields in [('2mass',2000,'fp_psc','j_m,j_cmsig,h_m,h_cmsig,k_m,k_cmsig,ph_qual,rd_flg,bl_flg,cc_flg'),
                    ('allwise',2010,'allwise_p3as_psd','w1mpro,w1sigmpro,w2mpro,w2sigmpro,w3mpro,w3sigmpro,w4mpro,w4sigmpro,ph_qual,cc_flags,ext_flg')]:
        ra,dec=position(r,year)
        sql=(f"SELECT designation,ra,dec,{fields} FROM {table} WHERE 1=CONTAINS("
             f"POINT('J2000',ra,dec),CIRCLE('J2000',{ra:.9f},{dec:.9f},0.0008333333))")
        yield sid+'_'+key+'.csv',tap(IRSA,sql)
    ra,dec=position(r,2007)
    req=dict(service='Mast.Galex.Catalog',params=dict(ra=ra,dec=dec,radius=3/3600),format='json')
    yield sid+'_galex.json',MAST+'?'+urllib.parse.urlencode({'request':json.dumps(req)})


def get(item):
    filename,url=item
    try:
        data=urllib.request.urlopen(url,timeout=35).read()
        if filename.endswith('.json'):
            payload=json.loads(data)
            if payload.get('status')!='COMPLETE':
                raise ValueError(str(payload.get('msg',payload.get('status'))))
        elif data.lstrip().startswith(b'<') or b',' not in data.splitlines()[0]:
            raise ValueError('unexpected response format')
        (HERE/'photometry'/filename).write_bytes(data)
        return dict(file=filename,url=url,sha256=hashlib.sha256(data).hexdigest(),size=len(data))
    except (OSError,ValueError) as err:
        return dict(file=filename,url=url,error=str(err))


if __name__=='__main__':
    (HERE/'photometry').mkdir(exist_ok=True)
    tasks=[]
    for r in csv.DictReader((HERE.parent/'gaia/source_diagnostics.csv').open()):
        if r['source_id'] in IDS: tasks.extend(jobs(r))
    with ThreadPoolExecutor(5) as ex: receipts=list(ex.map(get,tasks))
    (HERE/'photometry'/'inputs.json').write_text(json.dumps(receipts,indent=2)+'\n')
