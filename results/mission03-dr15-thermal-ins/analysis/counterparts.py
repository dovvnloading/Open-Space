"""Query nearby catalogue entries; an empty or failed query is never a flux limit."""
import argparse
import csv
import io
import json
import math
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

def separation(ra, dec, other_ra, other_dec):
    a, b = math.radians(dec), math.radians(other_dec)
    h = math.sin((a-b)/2)**2 + math.cos(a)*math.cos(b)*math.sin(math.radians(ra-other_ra)/2)**2
    return math.degrees(2*math.asin(min(1, math.sqrt(h))))*3600

def requests_for(ra, dec):
    def cone(r): return f"1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{dec},{r/3600:.8f}))"
    def tap(url, sql): return Request(url+'?'+urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=sql)))
    gaia = 'source_id,ra,dec,pmra,pmdec,parallax,parallax_error,phot_g_mean_mag,bp_rp,ruwe'
    wise = 'designation,ra,dec,sigra,sigdec,w1mpro,w1sigmpro,w2mpro,w2sigmpro,cc_flags,ph_qual'
    optical = 'ls_id,ra,dec,type,flux_g,flux_r,flux_z,flux_ivar_g,flux_ivar_r,flux_ivar_z,psfdepth_r,nobs_r,maskbits'
    query = f"SELECT {optical} FROM ls_dr10.tractor WHERE q3c_radial_query(ra,dec,{ra},{dec},{30/3600:.8f})"
    ps1 = 'objID,raStack,decStack,primaryDetection,nDetections,gPSFMag,rPSFMag,rPSFMagErr,iPSFMag'
    return {
        'gaia': tap('https://gea.esac.esa.int/tap-server/tap/sync',f'SELECT {gaia} FROM gaiadr3.gaia_source WHERE {cone(90)}'),
        'allwise': tap('https://irsa.ipac.caltech.edu/TAP/sync',f'SELECT {wise} FROM allwise_p3as_psd WHERE {cone(30)}'),
        'legacy': Request('https://datalab.noirlab.edu/query/query?'+urlencode(dict(sql=query,ofmt='csv',**{'async':'false'})),headers={'X-DL-AuthToken':'anonymous.0.0.anon_access'}),
        'ps1': Request('https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/stack.csv?'+urlencode(dict(ra=ra,dec=dec,radius=30/3600,pagesize=5000,columns='['+ps1+']')))
    }

def get(task):
    sid,survey,request,folder=task
    dest=folder/f'{sid}-{survey}.json'
    if dest.exists() and json.loads(dest.read_text())['status']=='ok': return
    record={'srcid':sid,'survey':survey,'url':request.full_url}
    for attempt in range(2):
        try:
            with urlopen(request,timeout=40) as response:
                payload=response.read().decode('utf-8-sig')
                content_type=response.headers.get_content_type()
            if '<' in payload[:30] or payload.startswith('Error') or 'csv' not in content_type and survey!='legacy':
                raise ValueError(f'Unexpected response {content_type}: {payload[:130]}')
            reader=csv.DictReader(io.StringIO(payload))
            expected={'gaia':'source_id','allwise':'designation','legacy':'ls_id','ps1':'objID'}[survey]
            if expected not in (reader.fieldnames or []):
                # MAST sends an empty body for no rows. Keep this distinct from a verified empty list.
                raise ValueError(f'No catalogue header: {payload[:100]}')
            record.update(status='ok',rows=list(reader))
            break
        except Exception as err:
            record.update(status='error',error=str(err));time.sleep(0.3)
    if record['status']=='ok': record.pop('error',None)
    dest.write_text(json.dumps(record,separators=(',',':'))+'\n')

def load_targets(path):
    with open(path,newline='') as f:
        for row in csv.DictReader(f):
            r={k.lower():v for k,v in row.items()}
            yield r['srcid'],float(r['ra']),float(r['dec'])

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('targets');p.add_argument('cache',type=Path)
    p.add_argument('--limit',type=int);p.add_argument('--workers',type=int,default=8);args=p.parse_args()
    args.cache.mkdir(parents=True,exist_ok=True)
    tasks=[]
    for sid,ra,dec in list(load_targets(args.targets))[:args.limit]:
        for survey,req in requests_for(ra,dec).items():
            if survey=='ps1' and dec < -31: continue
            tasks.append((sid,survey,req,args.cache))
    with ThreadPoolExecutor(args.workers) as ex:
        for _ in ex.map(get,tasks): pass
