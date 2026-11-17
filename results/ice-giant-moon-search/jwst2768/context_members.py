"""Read the 5x4 input layout and final drizzle parameters from real level-3 metadata."""
import argparse,csv,hashlib,re
from pathlib import Path
from astropy.io import fits
from mosaic_fetch import read_range

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('directory',type=Path);ap.add_argument('manifest',type=Path)
    a=ap.parse_args();rows=[]
    for name in sorted(a.directory.glob('jw02768-o*_i2d.fits')):
        url='https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:JWST/product/'+name.name
        with fits.open(name) as image:
            meta={}
            for kind in ('HDRTAB','ASDF'):
                h=image[kind];receipt=h.fileinfo();off=receipt['datLoc'];size=receipt['datSpan'];cache=name.with_suffix('.'+kind+'.bin')
                if not cache.exists():
                    raw,_=read_range(url,off,off+size-1);cache.write_bytes(raw)
                buf=cache.read_bytes()
                if len(buf)!=size:raise ValueError(kind)
                meta[kind.lower()+'_start']=off;meta[kind.lower()+'_bytes']=size;meta[kind.lower()+'_sha256']=hashlib.sha256(buf).hexdigest()
                if kind=='HDRTAB':
                    table=fits.BinTableHDU.fromstring(h.header.tostring().encode()+buf).data.copy()
                else:
                    yaml=buf.split(b'...\n')[0].decode('utf-8')
                    chunk=re.search(r'^  resample: \{(.*?)\}',yaml,re.M|re.S).group(1)
                    params={k.strip():v.strip() for k,v in (kv.split(':',1) for kv in chunk.split(','))}
            for idx,r in enumerate(table):
                dith=int(r['EXPOSURE']);assert dith==1+idx//4 and r['DETECTOR'].lower()==f'nrcb{1+idx%4}'
                rows.append(dict(mosaic=name.name,context_bit=idx,exposure=dith,detector=r['DETECTOR'],file=r['FILENAME'],exptime_s=r['EFFEXPTM'],
                                 pointings=params['pointings'],pixfrac=params['pixfrac'],kernel=params['kernel'],weight_type=params['weight_type'],**meta))
    if len(rows)!=60:raise ValueError('expected three products with five four-detector dithers')
    with a.manifest.open('w',newline='') as stream:
        wr=csv.DictWriter(stream,list(rows[0]),lineterminator='\n');wr.writeheader();wr.writerows(rows)
