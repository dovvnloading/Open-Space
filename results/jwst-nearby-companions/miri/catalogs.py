"""Freeze catalog answers for the fields; compare every first-epoch trigger."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

HERE=Path(__file__).resolve().parent
SERVICES={
    "GaiaDR3":("https://gea.esac.esa.int/tap-server/tap/sync","source_id,ra,dec,pmra,pmdec,parallax,phot_g_mean_mag","gaiadr3.gaia_source"),
    "SIMBAD":("https://simbad.cfa.harvard.edu/simbad/sim-tap/sync","main_id,ra,dec,otype,pmra,pmdec","basic"),
    "NASA-PS":("https://exoplanetarchive.ipac.caltech.edu/TAP/sync","distinct pl_name,hostname,ra,dec","ps"),
    }

def get_service(task):
    target,name,meta,out,refresh=task
    dest=out/f"{target}_{name}.csv"
    stamp=out/f"{target}_{name}.json"
    ra,dec=meta['epochs'][0]['ra'],meta['epochs'][0]['dec']
    url,columns,table=SERVICES[name]
    query=f"SELECT {columns} FROM {table} WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra:.10f},{dec:.10f},0.0125))"
    url+="?"+urlencode(dict(REQUEST="doQuery",LANG="ADQL",FORMAT="csv",QUERY=query))
    if refresh:
        with urlopen(url,timeout=100) as resp:
            raw=resp.read()
            if not raw or b'<html' in raw.lower() or raw.lstrip().startswith(b'<'):
                raise ValueError((name,raw[:200]))
            dest.write_bytes(raw)
            stamp.write_text(json.dumps(dict(query=query,url=url,response_sha256=hashlib.sha256(raw).hexdigest(),
                http_status=resp.status,retrieved_utc=datetime.now(timezone.utc).isoformat()),indent=2)+"\n")
    raw=dest.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=json.loads(stamp.read_text())['response_sha256']:
        raise ValueError(dest)
    return (target,name),list(csv.DictReader(io.StringIO(raw.decode('utf-8-sig'))))

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("measurements",type=Path)
    parser.add_argument("out",type=Path);parser.add_argument("--refresh",action="store_true")
    a=parser.parse_args();a.out.mkdir(exist_ok=True,parents=True)
    meta=json.loads((a.measurements/'measurement.json').read_text())
    jobs=[(t,n,m,a.out,a.refresh) for t,m in meta.items() for n in SERVICES]
    with ThreadPoolExecutor(6) as pool:
        cats=dict(pool.map(get_service,jobs))
    with (a.measurements/'sources.csv').open() as f:rows=list(csv.DictReader(f))
    result=[]
    for row in rows:
        target=row['target']
        p=SkyCoord(float(row['ra_early_deg']),float(row['dec_early_deg']),unit='deg')
        r=dict(target=target,id=row['id'])
        for name in SERVICES:
            entries=cats[target,name]
            closest=None;dist=np.inf
            for x in entries:
                ra,dec=float(x['ra']),float(x['dec'])
                if name in ('GaiaDR3','SIMBAD') and x.get('pmra') and x.get('pmdec'):
                    # CDS basic coordinates are at J2000.0, Gaia DR3 at J2016.0.
                    from astropy.time import Time
                    dt=Time(meta[target]['epochs'][0]['date']).jyear-(2016. if name=='GaiaDR3' else 2000.)
                    ra+=float(x['pmra'])*dt/(3600000*np.cos(np.deg2rad(dec)))
                    dec+=float(x['pmdec'])*dt/3600000
                d=p.separation(SkyCoord(ra,dec,unit='deg')).arcsec
                if d<dist:closest=x;dist=d
            key={'GaiaDR3':'source_id','SIMBAD':'main_id','NASA-PS':'pl_name'}[name]
            r[name+'_nearest_in_cone']=closest[key] if closest else ''
            r[name+'_nearest_arcsec']=round(dist,4) if closest else ''
            r[name+'_within_1arcsec']=int(dist<=1.)
        result.append(r)
    with (a.out/'crossmatch.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,list(result[0]));w.writeheader();w.writerows(result)

if __name__=='__main__': main()
