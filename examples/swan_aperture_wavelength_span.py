"""Measure true BC variation over extraction aperture; not a full response integral."""
import json,argparse,hashlib
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy import log
from scipy.ndimage import map_coordinates
log.setLevel('ERROR')
p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--calibration-root',required=True);p.add_argument('--cache',required=True);p.add_argument('--output',required=True);a=p.parse_args()
bc=fits.getdata(Path(a.calibration_root)/'ParameterFiles/20260604_SSDC_BC_Band5.fits',1);out=[]
for row in json.load(open(a.input)):
 if not(row['start'].startswith('2025-12-25') or row['start'].startswith('2026-01-02 10:')):continue
 with fits.open(Path(a.cache)/(hashlib.sha256(row['url'].encode()).hexdigest()+'.fits')) as h:
  yy,xx=np.indices(h['IMAGE'].data.shape);rr=np.hypot(xx-row['x'],yy-row['y']);w=WCS(h['IMAGE'].header,key='A');w.sip=None;spans={}
  for rad in [2.5,3.5,5.,8.]:
   m=rr<=rad;ax,ay=w.wcs_pix2world(xx[m],yy[m],0);v=map_coordinates(bc,[ay,2039-ax],order=1);spans[str(rad)]={'min_um':float(v.min()),'max_um':float(v.max()),'span_nm':float(np.ptp(v)*1000)}
  out.append(dict(start=row['start'],url=row['url'],spans=spans))
Path(a.output).write_text(json.dumps(out,indent=2)+'\n')
for e in out:print(e['start'],e['spans'])
