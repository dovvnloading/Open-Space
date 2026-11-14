"""Download public UVOT sky and exposure images and retain aligned sky cutouts."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from io import BytesIO
from urllib.request import urlopen
import traceback,gzip,hashlib,json,sys,os
import numpy as np
from astropy.table import Table
from astropy.wcs import WCS
from astropy.io import fits
from astropy.time import Time
BASE='https://nasa-heasarc.s3.amazonaws.com/swift/data/obs/'
D=Path(__file__).resolve().parent
POS=[21.06275,-73.32594444]
RECEIPTS = {r['obsid']: r for r in json.loads((D/'swift-receipts.json').read_text())}

def reduce(row):
    obs=str(row['DataLinkID']); date=Time(row['start_time'],format='mjd').isot[:7].replace('-','_')
    if obs not in RECEIPTS: return
    out=D/'swift_stamps'/obs
    if out.with_suffix('.json').exists():return
    prefix=BASE+date+'/'+obs+'/uvot/image/sw'+obs+'uw1_'
    record={'obsid':obs,'url':prefix+'sk.img.gz','planes':[]}
    try:
        with urlopen(record['url'],timeout=30) as h:b=h.read()
        record['sha256']=hashlib.sha256(b).hexdigest()
        if record['sha256'] != RECEIPTS[obs]['sha256_sky']:
            raise ValueError('Archive image checksum mismatch: ' + obs)
        arrays={}
        with fits.open(BytesIO(gzip.decompress(b)),memmap=False) as hd:
            for i,h in enumerate(hd):
                if h.data is None or h.data.ndim!=2:continue
                w=WCS(h.header).celestial
                if len(w.wcs.ctype)!=2:continue
                x,y=w.all_world2pix([POS],0)[0];ix,iy=np.rint([x,y]).astype(int)
                N=48;ny,nx=h.data.shape
                if ix<N or iy<N or ix+N>=nx or iy+N>=ny:continue
                stamp=h.data[iy-N:iy+N+1,ix-N:ix+N+1].copy()
                flat=stamp.ravel()
                rec=dict(idx=i,ci=[int(ix),int(iy)],dims=[nx,ny],target=[float(x),float(y)],wcs=w.to_header().tostring(),exposure=h.header['EXPOSURE'],date_obs=h.header['DATE-OBS'],date_end=h.header['DATE-END'],counts=float(flat.sum()))
                record['planes'].append(rec)
                arrays['image'+str(i)]=stamp
        if arrays:
            np.savez_compressed(out.with_suffix('.npz'),**arrays)
        out.with_suffix('.json').write_text(json.dumps(record))
    except Exception as e:
        raise RuntimeError(f'Could not recover Swift sky image {obs}') from e
if __name__=='__main__':
    (D/'swift_stamps').mkdir(exist_ok=True)
    t=Table.read(D/'inputs/swift-coverage.xml')
    rows=[r for r in t if r['uvot_expo_w1']>0]
    if len(sys.argv)>1:rows=rows[:int(sys.argv[1])]
    with ThreadPoolExecutor(10) as ex:list(ex.map(reduce,rows))
