"""Commanded Gaia DR3 field-of-view approaches; no image or detection claim."""
import csv,gzip,hashlib,io,json,sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
from datetime import datetime,timezone
import numpy as np
from astropy.time import Time
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'orbits'))
from propagate import family,sky

def query(k):
    sql=f'SELECT TOP 20000 jd_time,ra_fov{k},dec_fov{k},scan_angle_fov{k} FROM gaiadr3.commanded_scan_law WHERE ra_fov{k} BETWEEN 33.6 AND 36.8 AND dec_fov{k} BETWEEN -50.3 AND -48.5 ORDER BY jd_time'
    url='https://gea.esac.esa.int/tap-server/tap/sync?'+urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=sql));start=datetime.now(timezone.utc).isoformat()
    with urlopen(url,timeout=150) as f:b=f.read()
    rd=csv.DictReader(io.StringIO(b.decode()));assert 'jd_time' in rd.fieldnames,b[:400];rows=list(rd);assert len(rows)<20000
    (HERE/f'fov{k}.csv.gz').write_bytes(gzip.compress(b,mtime=0));return (dict(url=url,sha256=hashlib.sha256(b).hexdigest(),rows=len(rows),utc=start),rows)

def main():
    states=family(json.loads((HERE.parent/'orbits'/'orbit_input.json').read_text()));receipts=[];elements=[]
    with ThreadPoolExecutor(2) as ex:
        for k,(receipt,rows) in enumerate(ex.map(query,[1,2]),1):
            receipts.append(receipt);groups=[]
            for row in rows:
                t=float(row['jd_time']);
                if not groups or (t-float(groups[-1][-1]['jd_time']))*86400>90:groups.append([])
                groups[-1].append(row)
            for g in groups:
                at=np.array([float(x['jd_time']) for x in g]);time=Time(2455197.5,float(np.median(at)),format='jd',scale='tcb')
                a,d=sky(states,time.tt.jd-2451545.);aa=np.array([float(x[f'ra_fov{k}']) for x in g]);dd=np.array([float(x[f'dec_fov{k}']) for x in g]);sep=np.sqrt(((a[:,None]-aa[None,:])*np.cos(np.radians(d[:,None])))**2+(d[:,None]-dd[None,:])**2).min(axis=1)
                elements.append([k,time.utc.isot,len(g),round(float(sep.min()),5),int((sep<=.40).sum()),*np.round([a.min(),a.max(),d.min(),d.max()],6)])
    with (HERE/'nominal_passes.csv').open('w') as f:
        w=csv.writer(f);w.writerow(['fov','representative_utc','attitude_samples','closest_track_deg','tracks_within_0p40deg','ra_min','ra_max','dec_min','dec_max']);w.writerows(sorted(elements,key=lambda r:r[1]))
    (HERE/'scan_receipts.json').write_text(json.dumps(receipts,indent=2)+'\n');print(len(elements),sum(x[4]>0 for x in elements),flush=True)
if __name__=='__main__':main()
