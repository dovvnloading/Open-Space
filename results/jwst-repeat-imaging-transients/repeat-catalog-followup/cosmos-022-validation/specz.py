"""Search the official COSMOS spectroscopic compilation by coordinates and host ID."""
from hashlib import sha256
from urllib.request import urlopen
from io import BytesIO
from pathlib import Path
from astropy.io import fits
import numpy as np,json
D=Path(__file__).resolve().parent
root='https://media.githubusercontent.com/media/cosmosastro/speczcompilation/main/specz_compilation/'
url=root+'specz_compilation_COSMOS_DR1.1_unique.fits'
with urlopen(url,timeout=90) as f:b=f.read()
print(len(b),b[:8],flush=True)
with fits.open(BytesIO(b)) as h:
 a=h[1].data; cols={x.lower():x for x in a.names};print(a.names[:32],len(a),flush=True)
 ra,dec=a[cols['ra_corrected']],a[cols['dec_corrected']]
 sep=np.hypot((ra-150.122033)*np.cos(np.deg2rad(2.323386208)),dec-2.323386208)*3600
 idx=np.argsort(sep)[:5];c20=[k for k in a.names if 'cos20' in k.lower() and 'id' in k.lower()]
 for k in c20:idx=np.unique(np.r_[idx,np.where(a[k]==953994)[0]])
 def val(v):
  if isinstance(v,bytes):return v.decode().strip()
  if isinstance(v,np.generic):return v.item()
  return str(v)
 rows=[]
 for i in idx:
  row={k:val(a[k][i]) for k in a.names if any(q in k.lower() for q in ('ra','dec','specz','qf','cos20','survey','id_','flag','confiden','public_or_private','reference','catalog'))};row['separation_arcsec']=float(sep[i]);rows.append(row)
 out=dict(url=url,sha256=sha256(b).hexdigest(),bytes=len(b),rows_in_catalog=len(a),nearby=sorted(rows,key=lambda r:r['separation_arcsec']))
 (D/'specz.json').write_text(json.dumps(out,indent=2)+'\n')
 print(out['nearby'][0])
