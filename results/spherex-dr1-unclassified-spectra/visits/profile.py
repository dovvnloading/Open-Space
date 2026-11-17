"""Read all D4 footprints at one galaxy to trace its observed 3.4-um band."""
import argparse
import csv
import io
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from astropy.io import fits
from index_visits import approximate_xy, cut_url, read_url, sha
from measure import measure

SOURCE='J135607.93+300452.8'

def run(cache, science, out, ra=209.033074,dec=30.0813453):
    with (cache/'footprints-D4.csv').open() as f:
        rows=[]
        for r in csv.DictReader(f):
            x,y=approximate_xy(r['poly'],ra,dec)
            if 22<x<2016 and 22<y<2016 and any(f'/{rel}/' in r['uri'] for rel in ['dr1','qr3']):
                rows.append(r)
    with fits.open(science/'sci_D4.fits') as h:
        maps=[h[n].data for n in ('CWAVE','CBAND')]
        def fetch(r):
            url=cut_url(r['uri'],ra,dec)
            result=dict(url=url,release='qr3' if '/qr3/' in r['uri'] else 'dr1',mjd=r['time_bounds_lower'])
            path=cache/'cutouts'/(sha(url.encode())+'.fits')
            try:
                if not path.exists():
                    raw=read_url(url)
                    if not raw.startswith(b'SIMPLE  '):
                        raise ValueError('non-FITS response')
                    path.write_bytes(raw)
                raw=path.read_bytes()
                with fits.open(io.BytesIO(raw)) as hd:
                    result.update(sha256=sha(raw),**measure(hd,ra,dec,maps))
            except Exception as e:
                result['status']='read_error: '+str(e)
            return result
        with ThreadPoolExecutor(3) as pool:
            results=list(pool.map(fetch,rows))
    keys=list(dict.fromkeys(k for r in results for k in r))
    with out.open('w') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(results)
    print(SOURCE,len(rows),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('cache','science','out'):p.add_argument(n,type=Path)
    a=p.parse_args();run(a.cache,a.science,a.out)
