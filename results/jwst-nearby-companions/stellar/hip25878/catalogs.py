from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
import json,hashlib
out=Path(__file__).resolve().parent
r,d=82.869603,-3.692178
cone=f"1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{r},{d},0.02222))"
queries=[
 ('simbad','https://simbad.cfa.harvard.edu/simbad/sim-tap/sync',"SELECT main_id,ra,dec,pmra,pmdec,plx_value,sp_type,otype FROM basic WHERE "+cone),
 ('allwise','https://irsa.ipac.caltech.edu/TAP/sync',"SELECT designation,ra,dec,w1mpro,w1sigmpro,w1sat,w2mpro,w2sigmpro,w2sat,w3mpro,w3sigmpro,w3sat,w4mpro,w4sigmpro,ph_qual,cc_flags FROM allwise_p3as_psd WHERE "+cone),
 ('gaia','https://gea.esac.esa.int/tap-server/tap/sync',"SELECT source_id,ra,dec,pmra,pmdec,parallax,phot_g_mean_mag FROM gaiadr3.gaia_source WHERE "+cone),
 ('twomass','https://irsa.ipac.caltech.edu/TAP/sync',"SELECT designation,ra,dec,k_m,k_cmsig,ph_qual FROM fp_psc WHERE "+cone),
 ('exoplanet_archive','https://exoplanetarchive.ipac.caltech.edu/TAP/sync',"SELECT DISTINCT pl_name,hostname FROM ps WHERE hip_name='HIP 25878' OR "+cone)]
def get(t):
 n,b,q=t;url=b+'?'+urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=q));data=urlopen(url,timeout=50).read()
 if data.lstrip().startswith(b'<'):raise ValueError(data[:200])
 (out/(n+'.csv')).write_bytes(data);return dict(name=n,file=n+'.csv',url=url,query=q,sha256=hashlib.sha256(data).hexdigest())
with ThreadPoolExecutor(6) as pool:rows=list(pool.map(get,queries))
(out/'catalog_inputs.json').write_text(json.dumps(rows,indent=2)+'\n')
