"""Freeze complete development coadd arrays; never reads holdout flux."""
import hashlib,json,sys
from pathlib import Path
import numpy as np
from astropy.io import fits
src=Path(sys.argv[1]);records=json.loads(Path(sys.argv[2]).read_text());out=Path(sys.argv[3]);arrays={}
for r in records:
 if r['kind']=='hold':continue
 p=src/r['filename'];assert hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'];key=f"s{r['sdss_id']}_{r['kind']}"
 with fits.open(p) as h:
  d=h['COADD'].data
  for col in ('loglam','flux','ivar','and_mask'):arrays[key+'_'+col]=np.array(d[col])
np.savez_compressed(out,**arrays);print('spectra',len(arrays)//4,'bytes',out.stat().st_size,'sha256',hashlib.sha256(out.read_bytes()).hexdigest())
