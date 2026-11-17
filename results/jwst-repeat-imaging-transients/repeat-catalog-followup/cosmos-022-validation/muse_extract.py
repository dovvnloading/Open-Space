"""Extract a source-centered MUSE spectrum from the public 2015 COSMOS cube."""
import json,hashlib
from urllib.request import urlopen
from urllib.parse import urlencode
from pathlib import Path
from io import BytesIO
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
import pandas as pd

D=Path(__file__).resolve().parent
RA,DEC=150.122033,2.323386208
ID='ivo://eso.org/ID?ADP.2016-06-15T05:56:21.985'
url='https://dataportal.eso.org/dataPortal/soda/sync?'+urlencode({'ID':ID,'CIRCLE':f'{RA} {DEC} {2.3/3600:.8f}'})
with urlopen(url,timeout=100) as f:raw=f.read()
if not raw.startswith(b'SIMPLE'):raise ValueError(raw[:500])
print(len(raw),flush=True)
with fits.open(BytesIO(raw)) as h:
 print([(p.name,p.data.shape if p.data is not None else None) for p in h],flush=True)
 data=h[1].data.astype('f8'); hdr=h[1].header
 print(hdr.get('BUNIT'),hdr.get('CTYPE3'),hdr.get('CUNIT3'),flush=True)
 w=WCS(hdr,naxis=3)
 yy,xx=np.indices(data.shape[1:]);a,d=w.celestial.pixel_to_world_values(xx,yy)
 rad=np.hypot((a-RA)*np.cos(np.deg2rad(DEC)),d-DEC)*3600
 specw=w.spectral;lam=specw.pixel_to_world(np.arange(data.shape[0])).to_value(u.AA)
 ap=rad<=.7;sky=(rad>1.2)&(rad<2)
 spectrum=np.nansum(data[:,ap],axis=1)-ap.sum()*np.nanmedian(data[:,sky],axis=1)
 res=pd.DataFrame(dict(wavelength_A=lam,flux=spectrum))
 if len(h)>2 and h[2].data is not None and h[2].data.shape==data.shape:
  var=h[2].data
  res['formal_error']=np.sqrt(np.nansum(var[:,ap],axis=1)+ap.sum()**2*np.nanmedian(var[:,sky],axis=1)*np.pi/(2*sky.sum()))
 res.to_csv(D/'muse_spectrum.csv.gz',index=False,compression={'method':'gzip','mtime':0})
 (D/'muse_receipt.json').write_text(json.dumps(dict(url=url,raw_sha256=hashlib.sha256(raw).hexdigest(),bytes=len(raw),id=ID,unit=hdr.get('BUNIT'),wavelength=hdr.get('CTYPE3'),source_aperture_arcsec=.7,sky_annulus_arcsec=[1.2,2.0]),indent=2)+'\n')
