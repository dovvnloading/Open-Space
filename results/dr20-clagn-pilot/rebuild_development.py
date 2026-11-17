"""Rebuild development-only FITS containers for offline numerical replay.

These are derived arrays and are not byte-identical to archive FITS. The new
manifest carries derived hashes; the original archive hashes remain unchanged.
"""
import hashlib,json,sys
from pathlib import Path
import numpy as np
from astropy.io import fits
root=Path(__file__).resolve().parent
archive=root/'development-coadds.npz'; data=archive.read_bytes()
if len(data)!=3133051 or hashlib.sha256(data).hexdigest()!='8fb4d1f8c9eef159aede7ce6901391c18f0b3373612d79ec64d450db16f0c78c':
 raise SystemExit('Frozen development archive size or SHA-256 mismatch; no output written')
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True);source=json.loads((root/'manifest.json').read_text());arrays=np.load(archive,allow_pickle=False);derived=[]
for row in source:
 if row['kind']=='hold':continue
 key=f"s{row['sdss_id']}_{row['kind']}";cols=[]
 for name,fmt in [('loglam','E'),('flux','E'),('ivar','E'),('and_mask','J')]:cols.append(fits.Column(name=name,format=fmt,array=arrays[key+'_'+name]))
 table=fits.BinTableHDU.from_columns(cols,name='COADD');dest=out/row['filename'];fits.HDUList([fits.PrimaryHDU(),table]).writeto(dest,overwrite=True);r=dict(row);r['source_sha256']=r['sha256'];r['sha256']=hashlib.sha256(dest.read_bytes()).hexdigest();derived.append(r)
(out/'derived-manifest.json').write_text(json.dumps(derived,indent=2)+'\n');print('rebuilt',len(derived),'development spectra; no holdout arrays included')
