"""Exploratory delta-line response check; not a CO2 fluorescence or production model."""
from pathlib import Path
import json,hashlib,sys
import numpy as np
from scipy.ndimage import map_coordinates
from scipy.integrate import trapezoid
from astropy.io import fits
from astropy.wcs import WCS
from astropy import log
log.setLevel('ERROR')
import argparse
ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--output',required=True);ap.add_argument('--calibration-root',required=True);ap.add_argument('--cache-root',required=True);args=ap.parse_args()
cal=Path(args.calibration_root)
bc=fits.getdata(cal/'ParameterFiles/20260604_SSDC_BC_Band5.fits',1);bw=fits.getdata(cal/'ParameterFiles/20260604_SSDC_BW_Band5.fits',1);ch=fits.getdata(cal/'ChannelMaps/spectral_channels_17ch_2040.fits',1)
summary=np.genfromtxt(next((cal/'FidSpecSummary/SubChannel17/Band5').glob('*.csv')),delimiter=',',names=True)
rows=json.load(open(args.input));out=[]
for r in rows:
 if r.get('detector')!='D5' or 'error' in r:continue
 file=Path(args.cache_root)/(hashlib.sha256(r['url'].encode()).hexdigest()+'.fits')
 if not file.exists():continue
 with fits.open(file) as h:
  wa=WCS(h['IMAGE'].header,key='A');wa.sip=None;x,y=wa.wcs_pix2world(r['x'],r['y'],0)
  x=float(x);y=float(y)
 cen=float(map_coordinates(bc,[[y],[2039-x]],order=1)[0]);width=float(map_coordinates(bw,[[y],[2039-x]],order=1)[0]);channel=int(ch[int(round(y)),int(round(x))]);
 if not 1<=channel<=17:continue
 spec=np.genfromtxt(cal/f'MiddlePixel/SubChannel17/Band5/FiducialSpec_Band5_SubCh{channel:02d}_17_Cent.csv',delimiter=',',names=True)
 fidcen=float(summary['BandCenter_um'][channel-1]);wave=spec['Wavelength_um']+(cen-fidcen);R=spec['NormTransmission'];normal=trapezoid(R/wave,wave)
 # Relative photon-count response for delta f_lambda line in F_nu-calibrated image: lambda_line R(lambda_line) / integral(R/lambda dlambda). Arbitrary common normalization.
 def line(l):return float(l*np.interp(l,wave,R,left=0,right=0)/normal)
 z={k:r[k] for k in ['observationid','start','url','wavelength_um','apertures','diagnostic_outlier_inclusive']};z.update(active_x=float(x),active_y=float(y),map_bc_um=cen,map_bw_um=width,channel=channel,fiducial_bc_um=fidcen,shift_um=cen-fidcen,response426=line(4.26),response4267=line(4.267),response_grid={str(l):line(l) for l in np.round(np.arange(4.24,4.291,.001),3)})
 out.append(z)
fn=Path(args.output);fn.write_text(json.dumps(out,indent=2));print('saved',fn)
for r in out:print(r['start'],round(r['map_bc_um'],6),round(r['map_bw_um'],5),r['channel'],round(r['response426'],3),round(r['response4267'],3),(r['diagnostic_outlier_inclusive'].get('2.5') or {}).get('sum_mjysr'))
