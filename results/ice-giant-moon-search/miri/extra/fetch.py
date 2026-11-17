"""Read native CAL SCI, ERR and DQ arrays for outer Uranus and edge-of-chip Neptune."""
import argparse,csv,hashlib,io,sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
from astropy.io import fits
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from miri_audit import headers,source_range,write_csv


def read_product(row, cache, arrays):
    receipts={}
    _,primary,sh,start=headers(row,cache,receipts)
    data={}
    for name in ('SCI','ERR','DQ'):
        if name!='SCI':
            sh=fits.Header.fromfile(io.BytesIO(source_range(row,cache,start,2880,receipts)))
            if sh['EXTNAME']!=name:raise ValueError(sh['EXTNAME'])
            start+=2880
        if arrays:
            raw=source_range(row,cache,start,sh.data_size,receipts)
            im=np.frombuffer(raw,'>i4' if name=='DQ' else '>f4').reshape(sh['NAXIS2'],sh['NAXIS1'])
            if name=='DQ':im=im.astype('uint32')^np.uint32(sh.get('BZERO',0))
            data[name]=im
        start+=((sh.data_size+2879)//2880)*2880
    return data,receipts


def products():
    fs=list(csv.DictReader((ROOT/'miri/footprints.csv').open()))
    products={r['productFilename']:r for r in csv.DictReader((ROOT/'manifests/products.csv').open())}
    return [(f,products[f['file']]) for f in fs]


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('cache',type=Path);a=ap.parse_args();a.cache.mkdir(exist_ok=True,parents=True)
    sel=[r for f,r in products() if f['planet']=='Uranus' or f['visit']=='01604005001' or (f['proposal']=='01272' and f['planet_on_chip']=='1')]
    with ThreadPoolExecutor(10) as ex:res=list(ex.map(lambda r:read_product(r,a.cache,True)[1],sel))
    rec={k:v for d in res for k,v in d.items()}
    write_csv(Path(__file__).parent/'science_inputs.csv',['file','offset','length','sha256'],[dict(file=k[0],offset=k[1],length=k[2],sha256=v) for k,v in sorted(rec.items())])
    print(len(sel),len(rec),flush=True)
