"""Rebuild development-only FITS containers for offline numerical replay.

These are derived arrays and are not byte-identical to archive FITS. The new
manifest carries derived hashes; the original archive hashes remain unchanged.
"""
import hashlib,json,sys
from pathlib import Path
import numpy as np
from astropy.io import fits
root=Path(__file__).resolve().parent;out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True);source=json.loads((root/'manifest.json').read_text());arrays=np.load(root/'development-coadds.npz',allow_pickle=False);derived=[]
for row in source:
 if row['kind']=='hold':continue
 key=f"s{row['sdss_id']}_{row['kind']}";cols=[]
 for name,fmt in [('loglam','E'),('flux','E'),('ivar','E'),('and_mask','J')]:cols.append(fits.Column(name=name,format=fmt,array=arrays[key+'_'+name]))
 table=fits.BinTableHDU.from_columns(cols,name='COADD');dest=out/row['filename'];fits.HDUList([fits.PrimaryHDU(),table]).writeto(dest,overwrite=True);r=dict(row);r['source_sha256']=r['sha256'];r['sha256']=hashlib.sha256(dest.read_bytes()).hexdigest();derived.append(r)
(out/'derived-manifest.json').write_text(json.dumps(derived,indent=2)+'\n');print('rebuilt',len(derived),'development spectra; no holdout arrays included')
