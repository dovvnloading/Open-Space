from urllib.request import urlopen
from urllib.parse import urlencode
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import csv,sys,io,json
tap='https://jwst.esac.esa.int/server/tap/sync?'
coords={'a2744':(3.595,-30.40),'macs0416':(64.038,-24.068),'macs0647':(101.959,70.248),'elgordo':(15.735,-49.252),'a370':(39.97,-1.576)}
def run(item):
 k,(a,d)=item
 q=f"""SELECT TOP 10000 obs_id,t_min,t_max,calib_level,access_url,em_min,em_max,instrument_name FROM ivoa.ObsCore WHERE calib_level=3 AND instrument_name='NIRCAM/IMAGE' AND dataproduct_type='image' AND 1=CONTAINS(POINT('ICRS',{a},{d}),s_region) ORDER BY t_max DESC"""
 u=tap+urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=q))
 try:
  b=urlopen(u,timeout=100).read();Path(k+'.csv').write_bytes(b);t=list(csv.DictReader(io.StringIO(b.decode()))); print(k,len(t),[dict(i=x['obs_id'],m=x['t_max']) for x in t[:4]],flush=True)
 except Exception as e: print(k,str(e),e.read(800) if hasattr(e,'read') else '',flush=True)
list(ThreadPoolExecutor(5).map(run,coords.items()))
