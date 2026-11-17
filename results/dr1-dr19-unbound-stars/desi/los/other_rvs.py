"""Check official LAMOST DR12 LRS/MRS and SDSS Legacy spectral archives for new epochs."""
from concurrent.futures import ThreadPoolExecutor
from urllib.request import urlopen
from urllib.parse import urlencode, urljoin
from pathlib import Path
import csv, json, re, io, argparse, hashlib
OUT=Path(__file__).resolve().parent
def fetch(url):
    with urlopen(url,timeout=100) as r: return r.read()
def lam(pair):
    version,table=pair
    base='https://lamost.org/dr12/v2.0/'
    query="SELECT c.obsid,c.ra,c.dec FROM "+table+" c WHERE "+' OR '.join("spos(c.ra,c.dec) @ scircle '<(%sd,%sd),0.001388889d>'"%(ra,dec) for ra,dec in coords)+" LIMIT 100"
    if version=='MRS':
        ids=[r['source_id'] for r in csv.DictReader((OUT/'leads.csv').open())]
        query="SELECT obsid,ra,dec FROM med_catalogue WHERE gaia_source_id IN ("+','.join("'"+sid+"'" for sid in ids)+") LIMIT 100"
    # LAMOST returns an HTML page with the public result identifier, then paged JSON.
    opts=urlencode({'sql':query,'output.fmt':'html'})
    html=fetch(base+'sql/q?'+opts).decode()
    endpoint=re.search(r'url:"([^"]+output\.fmt=dbgrid)"',html)[1]
    raw=fetch(urljoin(base,endpoint)+'&rows=100&page=1')
    data=json.loads(raw)
    assert isinstance(data['rows'],list) and int(data['total'])==len(data['rows']),data
    return {'survey':'LAMOST DR12 v2.0 '+version,'endpoint':base+'sql/q','sql':query,'result':data}
def legacy(coords):
    box=' OR '.join(f"(ra between {float(ra)-.0073} and {float(ra)+.0073} and dec between {float(dec)-.0014} and {float(dec)+.0014})" for ra,dec in coords)
    base='https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch'
    q='select top 500 specObjID,ra,dec,z,zErr,class,mjd from SpecObjAll where '+box
    b=fetch(base+'?'+urlencode({'cmd':q,'format':'csv'}))
    assert b.startswith(b'#Table1\nspecObjID')
    return {'survey':'SDSS DR18 Legacy SpecObjAll','endpoint':base,'sql':q,'result':list(csv.DictReader(io.StringIO(b.decode().split('\n',1)[1])))}
if __name__=='__main__':
    with (OUT/'leads.csv').open() as f: coords=[(r['target_ra'],r['target_dec']) for r in csv.DictReader(f)]
    jobs={'LRS':lambda:lam(('LRS','catalogue')),'MRS':lambda:lam(('MRS','med_catalogue')),'legacy':lambda:legacy(coords)}
    def worker(name):
        outcome=jobs[name]()
        (OUT/('external_'+name+'.json')).write_text(json.dumps(outcome,indent=2)+'\n')
        print(name,str(outcome['result'])[:400],flush=True)
    with ThreadPoolExecutor(3) as pool:list(pool.map(worker,jobs))
