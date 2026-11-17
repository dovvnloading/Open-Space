"""Retrieve object-linked references, 2MASS-normalization catalogue, and full field astrometry."""
import argparse, csv, hashlib, io, json, time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import urlopen
from urllib.parse import urlencode
HOSTS={'HD53143':('5479222240596469632',104.9963,-61.3344),'HD38206':('2967287596507432704',85.8404,-18.5576)}
def main(out):
    out.mkdir(exist_ok=True,parents=True)
    jobs=[]
    for n,(sid,ra,dec) in HOSTS.items():
        jobs.extend([
            (n+'-references','https://simbad.cfa.harvard.edu/simbad/sim-tap/sync',"SELECT DISTINCT bibcode,title,doi FROM ident JOIN has_ref ON ident.oidref=has_ref.oidref JOIN ref ON ref.oidbib=has_ref.oidbibref WHERE ident.id='Gaia DR3 "+sid+"' ORDER BY bibcode"),
            (n+'-gaia','https://gea.esac.esa.int/tap-server/tap/sync',f"SELECT source_id,ref_epoch,ra,dec,ra_error,dec_error,parallax,parallax_error,pmra,pmra_error,pmdec,pmdec_error,phot_g_mean_mag,ruwe FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec), CIRCLE('ICRS',{ra},{dec},0.012))"),
            (n+'-2mass','https://irsa.ipac.caltech.edu/TAP/sync',f"SELECT designation,ra,dec,k_m,k_cmsig,ph_qual FROM fp_psc WHERE 1=CONTAINS(POINT('ICRS',ra,dec), CIRCLE('ICRS',{ra},{dec},0.005))"),
            (n+'-allwise','https://irsa.ipac.caltech.edu/TAP/sync',f"SELECT designation,ra,dec,w1mpro,w1sigmpro,w2mpro,w2sigmpro,w2sat,cc_flags FROM allwise_p3as_psd WHERE 1=CONTAINS(POINT('ICRS',ra,dec), CIRCLE('ICRS',{ra},{dec},0.005))")
        ])
    def get(job):
        name,base,query=job
        params=dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',MAXREC=10000,QUERY=query)
        url=base+'?'+urlencode(params)
        for i in range(3):
            try:
                with urlopen(url,timeout=65) as r: raw=r.read()
                if not raw.strip() or raw.lstrip().startswith(b'<'):raise ValueError(raw[:200])
                break
            except OSError:
                if i==2:raise
                time.sleep(2**i)
        count=len(list(csv.DictReader(io.StringIO(raw.decode('utf-8-sig')))))
        assert count<10000
        (out/(name+'.csv')).write_bytes(raw)
        return dict(name=name,url=url,rows=count,sha256=hashlib.sha256(raw).hexdigest(),utc=datetime.now(timezone.utc).isoformat())
    with ThreadPoolExecutor(5) as pool:
        fut=[pool.submit(get,x) for x in jobs]
        receipt=[]
        for j,f in zip(jobs,fut):
            try:receipt.append(f.result())
            except Exception as e:receipt.append(dict(name=j[0],error=str(e)))
    (out/'catalogue-receipts.json').write_text(json.dumps(receipt,indent=2)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('out',type=Path);main(p.parse_args().out)
