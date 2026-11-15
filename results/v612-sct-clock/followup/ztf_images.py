"""Public ZTF image cutouts for a predeclared phase-dependent source-localization test.

Select at most four best-error exposures per (year, eighth-cycle phase bin),
using only archive metadata and the frozen PR #24 ephemeris. No image intensities
are inspected at selection. A manifest records errors and each source URL.
"""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen,Request
import hashlib,gzip,json
import numpy as np,pandas as pd

RA,DEC=277.941193,-14.315539
P=.15956200154633035
EPOCH=2459499.9473282346
ROOT=Path(__file__).resolve().parent.parent
OUT=Path('/tmp/v612-ztf-images')
OUT.mkdir(parents=True,exist_ok=True)
a=pd.read_csv(ROOT/'data/ztf-field.csv')
a=a[(a.oid==333114400007906)&(a.filtercode=='zg')&(a.catflags==0)&(a.magerr<.09)&(a.mag>12)&(a.hjd>2459000)].copy()
a['yr']=np.floor(2000+(a.hjd-2451544.5)/365.25).astype(int)
a['phase']=((a.hjd+.00080-EPOCH)/P)%1
a['phase_bin']=np.floor(a.phase*8).astype(int)
a=pd.concat([g.sort_values(['magerr','hjd']).head(4) for _,g in a.groupby(['yr','phase_bin'])]).sort_values('hjd').head(120)
print('ZTF_SELECTION_SIZE',len(a),flush=True)

def fetch(row):
    d=str(int(row.filefracday))
    base=(f"https://irsa.ipac.caltech.edu/ibe/data/ztf/products/sci/{d[:4]}/{d[4:8]}/{d[8:]}/"
          f"ztf_{d}_{int(row.field):06d}_{row.filtercode}_c{int(row.ccdid,16):02d}_o_q{int(row.qid,16)}_sciimg.fits")
    url=base+f'?center={RA},{DEC}&size=35arcsec&gzip=false'
    try:
        with urlopen(Request(url,headers={'User-Agent':'Open-Space public archival analysis'}),timeout=65) as response:
            blob=response.read()
        if not blob.startswith(b'SIMPLE'):
            raise ValueError('Unexpected non-FITS response')
        name=f'{d}_{row.expid}.fits.gz'
        (OUT/name).write_bytes(gzip.compress(blob,mtime=0))
        return dict(expid=str(row.expid),hjd=float(row.hjd),phase=float(row.phase),mag=float(row.mag),
                    magerr=float(row.magerr),yr=int(row.yr),file=name,sha256=hashlib.sha256(blob).hexdigest(),
                    url=url,bytes=len(blob),status='ok')
    except Exception as e:
        return dict(expid=str(row.expid),url=url,status='error',message=repr(e))

with ThreadPoolExecutor(max_workers=4) as pool:
    rows=list(pool.map(fetch,(row for row in a.itertuples())))
(OUT/'manifest.json').write_text(json.dumps(rows,indent=2)+'\n')
print('ZTF_IMAGE_SUMMARY',json.dumps(dict(expected=len(rows),fetched=sum(r['status']=='ok' for r in rows),
                                           errors=[r for r in rows if r['status']!='ok'][:5])),flush=True)
url=f'https://mast.stsci.edu/tesscut/api/v0.1/sector?ra={RA}&dec={DEC}'
try:
    with urlopen(url,timeout=45) as r:
        sectors=json.load(r)
    (OUT/'tess-sectors.json').write_text(json.dumps({'url':url,'result':sectors},indent=2)+'\n')
    print('TESS_SECTORS',json.dumps(sectors),flush=True)
except Exception as e:
    (OUT/'tess-sectors.json').write_text(json.dumps({'url':url,'error':repr(e)})+'\n')
    print('TESS_SECTORS_ERROR',repr(e),flush=True)
