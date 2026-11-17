"""Retry unavailable IRSA cutouts with the corresponding original public AWS images."""
import sys,argparse,hashlib,json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from astropy.io import fits
from acquire import HERE,load_url

def repair(path,cache):
    ledger=json.loads(path.read_text());folder=cache/ledger['target']
    def fetch(item):
        p,kind=item; dest=folder/(p['id']+'_'+kind+'.fits');suffix=p[kind+'_url'].split('/p1bm_frm/')[1].split('?')[0]
        if not dest.exists():
            url='https://nasa-irsa-wise.s3.us-west-2.amazonaws.com/wise/neowiser/images/p1bm_frm/'+suffix
            try:
                raw=load_url(url)
                import io,gzip
                with fits.open(io.BytesIO(gzip.decompress(raw) if raw[:2]==bytes([31,139]) else raw)) as hdus: hdus.verify('exception'); hdus.writeto(dest)
                h=fits.getheader(dest);p[kind+'_original_url']=url;p[kind+'_original_sha256']=hashlib.sha256(raw).hexdigest();p[kind+'_sha256']=hashlib.sha256(dest.read_bytes()).hexdigest();p[kind+'_shape']=[h['NAXIS2'],h['NAXIS1']]
            except Exception as e:p[kind+'_fallback_error']=repr(e)
    with ThreadPoolExecutor(5) as ex:list(ex.map(fetch,[(p,k) for p in ledger['products'] for k in ('int','msk')]))
    path.write_text(json.dumps(ledger,indent=1)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--cache',type=Path,required=True);p.add_argument('targets',nargs='+');args=p.parse_args()
    for name in args.targets:repair(HERE/'receipts'/(name+'.json'),args.cache)
