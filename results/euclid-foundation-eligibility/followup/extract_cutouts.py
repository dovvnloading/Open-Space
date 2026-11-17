from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import json,pathlib,numpy as np
p=pathlib.Path('results/euclid-foundation-eligibility/followup/inputs'); out=pathlib.Path(r'D:\Open-Space\campaigns\euclid-foundation-20261116-70b26a3-terra61\cutouts');out.mkdir(exist_ok=True)
for i in json.loads((p/'cutout-requests.json').read_text()):
 url=i['url'].split('?')[0]
 with fits.open(url,use_fsspec=True,fsspec_kwargs={'block_size':1000000}) as h:
  hdr=h[0].header; x,y=WCS(hdr).world_to_pixel(SkyCoord(273.167,64.3687,unit='deg')); x0=int(x)-150;y0=int(y)-150;arr=h[0].section[y0:y0+300,x0:x0+300];hdr['CRPIX1']-=x0;hdr['CRPIX2']-=y0;dest=out/(i['id']+'-'+i['kind']+'.fits');fits.PrimaryHDU(arr,header=hdr).writeto(dest,overwrite=True);print(dest.name,arr.shape, np.isfinite(arr).sum(),flush=True)
