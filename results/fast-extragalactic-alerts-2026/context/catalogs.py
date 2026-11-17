"""Retrieve astrometric, morphology and known minor-planet context for the rapid-pair objects."""
import argparse
import csv
import gzip
import hashlib
import io
import json
import math
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone, timedelta
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parent
STORE = ROOT / 'queries'
GAIA = 'https://gea.esac.esa.int/tap-server/tap/sync'
PS1 = 'https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/mean.csv'
SB = 'https://ssd-api.jpl.nasa.gov/sb_ident.api'
HORIZONS = 'https://ssd.jpl.nasa.gov/api/horizons.api'

def positions():
    sources = {}
    for filename, names in [('rubin/pairs.csv', ('ra', 'dec', 't1')), ('ztf/derived/ztf_pairs.csv', ('meanra', 'meandec', 'mjd1'))]:
        with (ROOT.parent / filename).open() as f:
            for r in csv.DictReader(f):
                oid = r['oid']
                ra, dec, mjd = (float(r[key]) for key in names)
                if not oid.startswith('ZTF'):
                    mjd -= 37 / 86400
                if oid not in sources or mjd < sources[oid][2]:
                    sources[oid] = (ra, dec, mjd)
    return sources

def query(name, url, params):
    STORE.mkdir(exist_ok=True)
    record = STORE / (name+'.receipt.json')
    bodyfile = STORE / (name+'.txt.gz')
    url += '?' + urlencode(params)
    if record.exists() and bodyfile.exists():
        row=json.loads(record.read_text())
        if row.get('url') == url and row.get('status') == 200:
            return gzip.decompress(bodyfile.read_bytes())
    row = dict(url=url, retrieved_utc=datetime.now(timezone.utc).isoformat())
    try:
        with urlopen(Request(url,headers={'User-Agent':'Open-Space astronomy research/1.0'}),timeout=75) as response:
            body = response.read()
            row.update(status=response.status, content_type=response.headers.get('Content-Type'), bytes=len(body),sha256=hashlib.sha256(body).hexdigest())
        bodyfile.write_bytes(gzip.compress(body,mtime=0))
        return body
    except Exception as exc:
        row.update(status=getattr(exc,'code',None),error=str(exc))
    finally:
        record.write_text(json.dumps(row,indent=2)+'\n')

def separation(ra, dec, ra2, dec2):
    a,b = math.radians((ra2-ra)/2),math.radians((dec2-dec)/2)
    return 7200*math.degrees(math.asin(min(1,math.sqrt(math.sin(b)**2+math.cos(math.radians(dec))*math.cos(math.radians(dec2))*math.sin(a)**2))))

def catalogs(item):
    oid,(ra,dec,mjd)=item
    sql=(f"SELECT source_id,ra,dec,pmra,pmra_error,pmdec,pmdec_error,parallax,parallax_error,phot_g_mean_mag,bp_rp,ruwe "
         f"FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra:.10f},{dec:.10f},0.002778))")
    query('gaia_'+oid,GAIA,dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=sql))
    if dec > -30.0:
        cols='objID,raMean,decMean,nDetections,gMeanPSFMag,rMeanPSFMag,iMeanPSFMag,zMeanPSFMag,iMeanKronMag'
        query('ps1_'+oid,PS1,dict(ra=ra,dec=dec,radius=10/3600,columns='['+cols+']'))

def sex(x,ra=False):
    seconds=round(abs(x)/(15 if ra else 1)*3600,3)
    minutes,s=divmod(seconds,60); h,m=divmod(int(minutes),60)
    return ('M' if x<0 else '')+f'{h:02d}-{m:02d}-{s:06.3f}'

def small_body(item):
    oid,(ra,dec,mjd)=item
    epoch=datetime(1858,11,17)+timedelta(days=mjd)
    params={'mpc-code':'I41' if oid.startswith('ZTF') else 'X05','obs-time':epoch.isoformat(timespec='seconds'),
            'fov-ra-center':sex(ra,True),'fov-dec-center':sex(dec),
            'fov-ra-hwidth': round(0.01/max(0.05,math.cos(math.radians(dec))),6),'fov-dec-hwidth':0.01,
            'two-pass':'true','suppress-first-pass':'true','req-elem':'false'}
    query('jpl_'+oid,SB,params)

def horizons():
    items = []
    for path in sorted((ROOT.parent/'ztf/data/ztf').glob('*_detections.csv.gz')):
        with gzip.open(path,'rt') as f:
            for row in csv.DictReader(f):
                if row['oid'] in ['ZTF26acgcmxl','ZTF26acgrjsd']:
                    items.append(row)
    items.sort(key=lambda row:float(row['mjd']))
    if not items: raise ValueError('Missing ZTF input detections')
    times=sorted({float(r['mjd'])+2400000.5 for r in items})
    raw=query('horizons_18150',HORIZONS,dict(format='json',COMMAND="'18150;'",CENTER="'I41@399'",EPHEM_TYPE="'OBSERVER'",
               TLIST="'"+' '.join(f'{t:.8f}' for t in times)+"'",ANG_FORMAT="'DEG'",QUANTITIES="'1'",EXTRA_PREC="'YES'",
               CAL_FORMAT="'JD'",CSV_FORMAT="'YES'",OBJ_DATA="'NO'"))
    report=json.loads(raw)['result']; lines=report.split('$$SOE')[1].split('$$EOE')[0].strip().splitlines()
    predictions=[(float(r[0]),float(r[-3]),float(r[-2])) for r in csv.reader(lines)]
    with (ROOT/'horizons_18150.csv').open('w') as f:
        writer=csv.writer(f);writer.writerow(['oid','candid','mjd','ra','dec','ephem_ra','ephem_dec','separation_arcsec'])
        for r in items:
            jd=float(r['mjd'])+2400000.5
            t,ra,dec=min(predictions,key=lambda p:abs(p[0]-jd))
            if abs(t-jd)>0.00001:raise ValueError('No matching calculated epoch')
            writer.writerow([r['oid'],r['candid'],r['mjd'],r['ra'],r['dec'],ra,dec,separation(float(r['ra']),float(r['dec']),ra,dec)])

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=['catalogs','small-body','horizons'])
    mode=parser.parse_args().mode
    if mode=='catalogs':
        with ThreadPoolExecutor(3) as pool:list(pool.map(catalogs,positions().items()))
    elif mode=='small-body':
        for item in sorted(positions().items(),key=lambda x:not x[0].startswith('ZTF')):small_body(item)
    else:horizons()
