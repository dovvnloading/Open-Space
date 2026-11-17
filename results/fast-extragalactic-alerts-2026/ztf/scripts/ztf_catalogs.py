"""Independent Gaia DR3 and SIMBAD cone searches for the nine slope-selected IDs."""
import csv
import datetime as dt
import gzip
import hashlib
import json
import time
from concurrent.futures import ThreadPoolExecutor
from urllib.request import urlopen
from urllib.parse import urlencode
from ztf_tap import RAW,META,ROOT

def run(what):
    name,endpoint,q=what
    if (META/(name+'.json')).exists(): return
    u=endpoint+'?'+urlencode(dict(request='doQuery',lang='ADQL',format='csv',query=q))
    for attempt in range(3):
        try:
            with urlopen(u,timeout=75) as h:
                b=h.read(); status=h.status
                if b.lstrip().startswith(b'<'): raise ValueError(b[:500])
            break
        except (OSError,ValueError):
            if attempt==2: raise
            time.sleep(2)
    out=RAW/(name+'.csv.gz')
    out.write_bytes(gzip.compress(b,mtime=0))
    info=dict(url=u,query=q,method='GET',http_status=status,bytes=len(b),sha256=hashlib.sha256(b).hexdigest(),retrieved_utc=dt.datetime.now(dt.timezone.utc).isoformat(),file=str(out.relative_to(ROOT)))
    (META/(name+'.json')).write_text(json.dumps(info,indent=2)+'\n')
    print(name,len(b),flush=True)
def main():
    with (ROOT/'derived'/'ztf_pairs.csv').open() as f: rows=list(csv.DictReader(f))
    queries=[]
    for r in rows:
        ra,dec=map(float,(r['meanra'],r['meandec']))
        cone=f"1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra:.9f},{dec:.9f},{5/3600:.10f}))"
        queries.extend([
            (r['oid']+'_gaiadr3','https://gea.esac.esa.int/tap-server/tap/sync','SELECT source_id,ra,dec,parallax,parallax_error,pmra,pmra_error,pmdec,pmdec_error,phot_g_mean_mag,bp_rp,ruwe FROM gaiadr3.gaia_source WHERE '+cone),
            (r['oid']+'_simbad','https://simbad.harvard.edu/simbad/sim-tap/sync','SELECT main_id,otype,ra,dec FROM basic WHERE '+cone)])
    with ThreadPoolExecutor(max_workers=3) as p:list(p.map(run,queries))
if __name__=='__main__':main()
