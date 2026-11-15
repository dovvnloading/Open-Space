"""Rebuild the compact median bins directly from the MAST science ECSV."""
import csv
import hashlib
import io
import sys
from pathlib import Path
from urllib.request import urlopen
import numpy as np

HERE = Path(__file__).resolve().parent
NAME = 'jw12492-o001_t001_nirspec_clear-prism-s1600a1-sub512_whtlt.ecsv'
URL = 'https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:JWST/product/' + NAME
if __name__ == '__main__':
    buf = Path(sys.argv[1]).read_bytes() if len(sys.argv)>1 else urlopen(URL, timeout=60).read()
    print('MAST product SHA256:', hashlib.sha256(buf).hexdigest())
    col = np.loadtxt(line for line in io.BytesIO(buf) if line[:1].isdigit())
    utc,t,y = col.T
    ids = np.floor((t - t[0])*1440).astype(int)
    with (HERE/'one_minute.csv').open('w') as f:
        w = csv.writer(f); w.writerow(['mjd_utc','bjd_tdb_minus_2400000_5','relative_flux','count'])
        norm = np.nanmedian(y)
        for i in np.unique(ids):
            good = (ids==i) & np.isfinite(y)
            if good.sum()>=10:
                w.writerow([f'{np.mean(utc[good]):.10f}',f'{np.mean(t[good]):.10f}',f'{np.median(y[good])/norm:.9f}',int(good.sum())])
