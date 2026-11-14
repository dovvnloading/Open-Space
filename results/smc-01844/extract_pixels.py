"""Rebuild the released TESS photometry and block cross-products from TessCut FITS."""
import argparse
import hashlib
import io
import json
import zipfile
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from replay import HERE, EPOCH
from tess import SECTORS, FREQ, clock_matrix, detrend_blocks


def extract(path, sector, out):
    with fits.open(path, memmap=False) as hdus:
        data=hdus[1].data
        time=data['TIME'].astype(float)+7000-.0008
        cube=data['FLUX'].astype(float)
        good=(data['QUALITY']==0)&np.isfinite(time)&np.all(np.isfinite(cube),axis=(1,2))
        t,cube=time[good],cube[good]
        xy=WCS(hdus[2].header).world_to_pixel_values(12.57704804,-72.64897851)
    ix,iy=np.rint(xy).astype(int)
    flux=cube[:,iy,ix]-np.median(cube[:,[0,0,6,6],[0,6,0,6]],axis=1)
    np.savez_compressed(out/f's{sector}.npz',t=t,flux=flux)
    s=2*np.pi*(t-EPOCH)*FREQ
    x=np.column_stack((clock_matrix(t,6),np.cos(s),np.sin(s)))
    x,blocks=detrend_blocks(t,x)
    y,_=detrend_blocks(t,cube.reshape(len(t),-1))
    gram,cross=[],[]
    for k in np.unique(blocks):
        use=blocks==k
        if use.sum()>=25 and np.ptp(t[use])>=.35:
            gram.append(x[use].T@x[use])
            cross.append(x[use].T@y[use])
    return np.array(gram),np.array(cross)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('output',type=Path)
    ap.add_argument('--fits-dir',type=Path,help='directory with 1844-sN.fits, skips download')
    args=ap.parse_args()
    args.output.mkdir(exist_ok=True,parents=True)
    metadata=json.loads((HERE/'sources.json').read_text())['TESS']
    products={}
    for s in SECTORS:
        record=metadata[str(s)]
        path=(args.fits_dir or args.output)/f'1844-s{s}.fits'
        if not args.fits_dir and not path.exists():
            with urlopen(record['url'],timeout=240) as response:
                with zipfile.ZipFile(io.BytesIO(response.read())) as z:
                    names=[n for n in z.namelist() if n.endswith('.fits')]
                    if len(names)!=1:
                        raise ValueError(names)
                    path.write_bytes(z.read(names[0]))
        print(s,hashlib.sha256(path.read_bytes()).hexdigest(),flush=True)
        products[f'x{s}'],products[f'y{s}']=extract(path,s,args.output)
    np.savez_compressed(args.output/'pixel-blocks.npz',**products)


if __name__=='__main__':
    main()
