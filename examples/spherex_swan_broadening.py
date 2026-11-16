"""Descriptive same-image ePSF broadening; no discovery/significance inference."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates,gaussian_filter
from scipy.optimize import least_squares
p=argparse.ArgumentParser();p.add_argument('--manifest',required=True);p.add_argument('--cache',required=True);p.add_argument('--output',required=True);a=p.parse_args()
res=[]
for row in json.load(open(a.manifest)):
 if not row['start'].startswith('2025-12-25'):continue
 path=Path(a.cache)/(hashlib.sha256(row['url'].encode()).hexdigest()+'.fits')
 with fits.open(path) as h:
  im=h['IMAGE'].data;var=h['VARIANCE'].data;fl=h['FLAGS'].data;w=WCS(h['IMAGE'].header,h);x,y=w.all_world2pix([[row['ra'],row['dec']]],0)[0];ax,ay=WCS(h['IMAGE'].header,h,key='A').wcs_pix2world([[x,y]],0)[0]
  tb=h['EPSF'].data;k=np.argmin((tb['XCENTER']-ax)**2+(tb['YCENTER']-ay)**2);ps=np.pad(tb['EPSF'][k].astype(float),24);center=40
  yy,xx=np.indices(im.shape);r=np.hypot(xx-x,yy-y);good=np.isfinite(im)&np.isfinite(var)&(var>0)&((fl&~(2097152|524288))==0);m=r<=2.5
  out={'url':row['url'],'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'start':row['start'],'wave':row.get('wavelength_um'),'n':int(m.sum())}
  if (~good&m).any():out['rejected']='invalid core';res.append(out);continue
  X=xx[m];Y=yy[m];z=im[m];s=np.sqrt(var[m]);bg=float(np.median(im[good&(r>8)&(r<12)]))
  def model(v):
   amp,b,dx,dy=v[:4];sigma=v[4] if len(v)==5 else 0
   kernel=gaussian_filter(ps,sigma*5,mode='constant',truncate=4) if sigma>0 else ps
   return amp*map_coordinates(kernel,[(Y-y-dy)*5+center,(X-x-dx)*5+center],order=3,mode='constant')*25+b
  lo=[0,-np.inf,-.5,-.5];hi=[np.inf,np.inf,.5,.5];v0=[max(float((z-bg).sum()),.01),bg,0,0]
  def fit(init,low,high):return least_squares(lambda v:(model(v)-z)/s,init,bounds=(low,high),max_nfev=1000,ftol=1e-10,xtol=1e-10,gtol=1e-10)
  point=fit(v0,lo,hi);broad=min([fit(list(point.x)+[v],lo+[0],hi+[.75]) for v in [.05,.25,.6]],key=lambda f:np.sum(f.fun**2))
  out.update(point={'parameters':point.x.tolist(),'chi2':float((point.fun**2).sum())},broadened={'parameters':broad.x.tolist(),'chi2':float((broad.fun**2).sum())},delta_chi2=float((point.fun**2).sum()-(broad.fun**2).sum()))
  out['profile']=[]
  for sig in [0,.1,.2,.3,.4,.5,.6,.75]:
   def f(v):return (model(list(v)+[sig])-z)/s
   q=least_squares(f,broad.x[:4],bounds=(lo,hi));out['profile'].append({'sigma_pixels':sig,'chi2':float((q.fun**2).sum())})
  res.append(out)
Path(a.output).write_text(json.dumps({'parameter_order':['amplitude','background','dx','dy','gaussian_sigma_native_pixels'],'caveats':['OUTLIER-inclusive only; strict rejected','Same images, not independent observation','Correlated noise and model uncertainty excluded; delta chi2 not significance','EPSF zero-padded24 oversampled pixels; Gaussian truncated4sigma; unknown original wings not reconstructed','Broadening is descriptive, not resolved gas; fixedradius2.5; sigma capped.75px'], 'results':res},indent=2))
for o in res:print(o['wave'],o.get('point'),o.get('broadened'),o.get('delta_chi2'))
