"""Query catalogues at the western XDF line. Save the exact response or the HTTP error."""
import csv,hashlib,json,sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib import request,parse
import xml.etree.ElementTree as ET
from urllib.error import HTTPError
pos='53.18577725,-27.78513879'
def tap(url,sql):
    return url+'?'+parse.urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=sql))
def fetch(task):
    name,url=task
    try:
        with request.urlopen(url,timeout=50) as r:
            d=r.read();status=r.status
        return dict(name=name,url=url,status=status,sha256=hashlib.sha256(d).hexdigest(),response=d.decode()[:12000],bytes=len(d))
    except HTTPError as e:
        return dict(name=name,url=url,status=e.code,error=e.read(1500).decode(errors='replace'))
    except Exception as e:return dict(name=name,url=url,error=str(e))
if __name__=='__main__':
    dl='https://datalab.noirlab.edu/tap/sync'
    sql=f"SELECT targetid,target_ra,target_dec,z,zwarn,spectype FROM desi_dr1.zpix WHERE 1=CONTAINS(POINT('ICRS',target_ra,target_dec),CIRCLE('ICRS',{pos},0.003))"
    casda='https://casda.csiro.au/casda_vo_tools/tap/sync'
    casql="SELECT TOP 50 table_name FROM TAP_SCHEMA.tables WHERE table_name LIKE '%racs%island%'"
    viz='https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?'+parse.urlencode({'-source':'J/ApJS/255/30','-c':'53.18577725 -27.78513879','-c.rs':15,'-out.max':30})
    tests=[('desi_dr1',tap(dl,sql)),('vlass_epoch1',viz),('RACS-table-discovery',tap(casda,casql))]
    if len(sys.argv)>2:
        racql=f"SELECT TOP 50 ra_deg_cont,dec_deg_cont FROM {sys.argv[2]} WHERE 1=CONTAINS(POINT('ICRS',ra_deg_cont,dec_deg_cont),CIRCLE('ICRS',{pos},0.004167))"
        tests.append(('RACS',tap(casda,racql)))
    with ThreadPoolExecutor(max_workers=4) as ex:data=list(ex.map(fetch,tests))
    Path(sys.argv[1]).write_text(json.dumps(data,indent=2)+'\n')
