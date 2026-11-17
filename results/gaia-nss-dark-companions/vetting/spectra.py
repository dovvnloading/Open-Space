"""Measure differential line motion directly on LAMOST DR11 extreme-RV nightly coadds."""
import argparse
import csv
import gzip
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

import numpy as np
from scipy.ndimage import percentile_filter, gaussian_filter1d
from scipy.optimize import minimize_scalar
from scipy.signal import find_peaks

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'epochs'))
from catalog_io import read_header, table_dtype

URL='https://www.lamost.org/dr11/v2.0/medspectrum/fits/'
# The two nights containing the most negative and most positive catalogued exposures.
IDS=('872107188','808707188')
WINDOWS=((4950,5000),(5000,5075),(5075,5150),(5150,5225),(5225,5300),(5300,5340))
C=299792.458


def spectrum(path, arm="B"):
    with gzip.open(path) as f:
        primary=read_header(f)
        if primary['HELIO']!='T' or primary['VACUUM']!='T' or primary['NAXIS']!='0':
            raise ValueError('wavelength frame differs from DR11')
        h=read_header(f)
        if arm=='R':
            assert h['EXTNAME']=='COADD_B' and int(h['NAXIS2'])==1
            f.read((int(h['NAXIS1'])+2879)//2880*2880)
            h=read_header(f)
        assert h['EXTNAME']=='COADD_'+arm and int(h['NAXIS2'])==1
        data=np.frombuffer(f.read(int(h['NAXIS1'])),dtype=table_dtype(h))[0]
    good=(data['IVAR']>0)&(data['ORMASK']==0)&np.isfinite(data['FLUX'])&(data['FLUX']>0)
    x=data['WAVELENGTH'].astype(float); flux=data['FLUX'].astype(float)
    # Estimate the pseudo-continuum sufficiently far above the ~0.7 A resolution width.
    flux=np.interp(x,x[good],flux[good])
    continuum=gaussian_filter1d(percentile_filter(flux,percentile=85,size=101),20)
    good &= np.isfinite(continuum)&(continuum>0)
    return x,1-flux/continuum,good


def run(cache):
    cache.mkdir(parents=True,exist_ok=True)
    data=[]; meta=[]
    for obsid in IDS:
        path=cache/(obsid+'.fits.gz')
        if not path.exists():
            path.write_bytes(urllib.request.urlopen(URL+obsid,timeout=90).read())
        meta.append(dict(obsid=obsid,url=URL+obsid,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        data.append(spectrum(path))
    x1,y1,ok1=data[0]; x2,y2,ok2=data[1]
    grid=np.arange(-140,141)
    traces={'dv_kms':grid}
    results=[]
    for lo,hi in (*WINDOWS,(4950,5340)):
        take=(x2>=lo)&(x2<=hi)&ok2
        # Use the same valid pixels for all lags, avoiding interpolation across any bad pixel.
        for v in grid[::2]:
            take &= np.interp(x2/(1+v/C),x1,ok1.astype(float))>=.999
        x,y=x2[take], y2[take]; y=y-y.mean()
        def ccf(v):
            t=np.interp(x/(1+v/C),x1,y1); t-=t.mean()
            return np.dot(t,y)/(np.linalg.norm(t)*np.linalg.norm(y))
        scores=np.array([ccf(v) for v in grid])
        j=np.argmax(scores); best=minimize_scalar(lambda v:-ccf(v),bounds=(grid[j]-1,grid[j]+1),method='bounded')
        peaks,_=find_peaks(scores,prominence=.05)
        traces[f'{lo}_{hi}']=scores
        results.append(dict(lo_A=lo,hi_A=hi,pixels=int(take.sum()),dv_kms=best.x,r_max=-best.fun,
                         local_peaks_kms=[int(grid[k]) for k in peaks]))
    with (HERE/'blue_ccf.csv').open('w',newline='') as f:
        w=csv.writer(f,lineterminator='\n'); w.writerow(traces)
        w.writerows(zip(*traces.values()))
    return dict(source_id='3378588057203660160',blue_coadds_low_then_high=meta,windows=results)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cache',type=Path,required=True)
    a=p.parse_args()
    (HERE/'spectral_motion.json').write_text(json.dumps(run(a.cache),indent=2)+'\n')
