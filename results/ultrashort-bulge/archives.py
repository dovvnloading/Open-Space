"""Retrieve public image footprints and predicted moving objects at selected peaks."""
import argparse
import csv
import hashlib
import io
import json
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
from pathlib import Path
from run_grid import HERE


def coords(ra,dec):
    def convert(x):
        v=x.lstrip('+-').split(':')
        return (float(v[0])+float(v[1])/60+float(v[2])/3600)*(-1 if x.startswith('-') else 1)
    return 15*convert(ra),convert(dec)


def candidates():
    out={}
    for p in (HERE/'kmt').glob('catalog*.tsv'):
        for r in csv.DictReader(p.open(),delimiter='\t'):
            if r['t_0'] not in ('-',''):
                out[r['Event']]=(*coords(r['RA'],r['Dec']),float(r['t_0'])+2400000)
    for r in csv.DictReader((HERE/'ogle/catalog.csv').open()):
        out['OGLE-'+r['event']]=(*coords(r['ra_hms'],r['dec_dms']),float(r['t0_hjd']))
    for r in csv.DictReader((HERE/'other/index.csv').open()):
        try:
            ra,dec=coords(r['ra'],r['dec'])
            mt=Path('.survey-cache/records')/(r['event']+'.json')
            if mt.exists():out[r['event']]=(ra,dec,float(json.loads(mt.read_text()).get('page_t0_jd') or 0))
        except (KeyError,ValueError,IndexError):continue
    return out


def query(job):
    name,svc,url=job
    result={'event':name,'service':svc,'queried_utc':datetime.now(timezone.utc).isoformat(),'url':url}
    try:
        with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'Open-Space astronomy research'}),timeout=45) as f:
            b=f.read();result['status']=f.status
        result.update(bytes=len(b),sha256=hashlib.sha256(b).hexdigest())
        txt=b.decode('utf-8','replace')
        if svc=='SkyBoT':
            result['response']=txt[:4800]
        elif txt.lstrip().startswith('<'):
            result['service_response']=txt[:1700]
        else:
            rows=list(csv.DictReader(io.StringIO(txt)))
            result['rows']=rows
            if not rows:result['header']=txt[:180]
    except Exception as ex:
        result['error']=str(ex)
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('events',nargs='+');p.add_argument('--out',type=Path,required=True)
    a=p.parse_args(); targets=candidates()
    jobs=[]
    for name in a.events:
        ra,dec,t=targets[name]
        tap='https://mast.stsci.edu/vo-tap/api/v0.1/caom/sync?'
        sql=("SELECT obs_collection,obs_id,filters,t_min,t_max,calib_level FROM dbo.ObsPointing "
             f"WHERE obs_collection IN ('HST','HLA','JWST') AND 1=CONTAINS(POINT('ICRS',{ra},{dec}),s_region)")
        jobs.append((name,'MAST',tap+urllib.parse.urlencode({'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'csv','QUERY':sql})))
        fp=f"1=CONTAINS(POINT('ICRS',{ra:.9f},{dec:.9f}),s_region)"
        for service,base,extra in [('Euclid','https://eas.esac.esa.int/tap-server/tap/sync',''),('CADC','https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/argus/sync'," AND (facility_name IN ('HST','JWST') OR obs_collection IN ('HST','JWST','HSTHLA'))")]:
            sql=f"SELECT TOP 120 obs_collection,obs_id,instrument_name,t_min,t_max FROM ivoa.obscore WHERE {fp} AND dataproduct_type='image'"+extra
            jobs.append((name,service,base+'?'+urllib.parse.urlencode({'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'csv','QUERY':sql})))
        url='https://ssp.imcce.fr/webservices/skybot/api/conesearch.php?'
        pars={'-ep':f'{t:.5f}','-ra':f'{ra:.7f}','-dec':f'{dec:.7f}','-rd':0.01,'-mime':'text','-loc':'500','-output':'basic'}
        jobs.append((name,'SkyBoT',url+urllib.parse.urlencode(pars)))
    with ThreadPoolExecutor(3) as pool:
        a.out.write_text(json.dumps(list(pool.map(query,jobs)),indent=2)+'\n')


if __name__=='__main__':main()
