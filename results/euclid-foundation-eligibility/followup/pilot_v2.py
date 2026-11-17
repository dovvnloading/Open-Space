"""Bounded, development-only two-visit aperture and control pilot."""
import argparse, csv, json
from pathlib import Path
import numpy as np
from scipy.ndimage import maximum_filter
from astropy.io import fits

ap=argparse.ArgumentParser();ap.add_argument('--cache',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=True)
ids=['102158276_013663548_VIS','102158276_028106832_VIS'];data=[]
allow=(1<<18)|(1<<24)
for ident in ids:
 images=[fits.getdata(args.cache/(ident+'-'+k+'.fits')) for k in ('science','rms','flag')]
 a,e,f=images;good=np.isfinite(a)&np.isfinite(e)&(e>0)&((f.astype(np.int64)&~allow)==0);data.append((a,e,good))
a,e,good=data[0];N=a.shape[0];peaks=np.argwhere(good&(a>8*e)&(a==maximum_filter(a,size=3)));selected=[]
yy,xx=np.mgrid[-8:9,-8:9];disk=xx*xx+yy*yy<=16; ring=(xx*xx+yy*yy>=36)&(xx*xx+yy*yy<=64)
for y,x in sorted(peaks,key=lambda q:float(a[tuple(q)]),reverse=True):
 if min(x,y,N-1-x,N-1-y)<10:continue
 if not good[y-8:y+9,x-8:x+9][disk].all():continue
 if any((x-u)**2+(y-v)**2<100 for u,v in selected):continue
 selected.append((int(x),int(y)))
rows=[]
for x,y in selected:
 measurements=[]
 for image,rms,mask in data:
  sl=np.s_[y-8:y+9,x-8:x+9];m=mask[sl];valid=m[disk].all() and m[ring].sum()>=50
  if not valid:measurements.append(None);continue
  sub=image[sl];sky=np.median(sub[ring&m]);pix=sub[disk]-sky;flux=float(pix.sum());err=float(np.sqrt((rms[sl][disk]**2).sum()));positive=np.maximum(pix,0);radius=float(np.sqrt(np.sum(positive*(xx[disk]**2+yy[disk]**2))/np.sum(positive))) if positive.sum()>0 else None
  measurements.append((flux,err,radius))
 if all(m is not None and m[0]>0 and m[2] is not None for m in measurements):
  m,n=measurements;rows.append({'x':x,'y':y,'flux1':m[0],'err1':m[1],'shape1':m[2],'flux2':n[0],'err2':n[1],'shape2':n[2]})
scale=float(np.median([r['flux2']/r['flux1'] for r in rows])) if rows else None
for r in rows:
 r['diff']=r['flux2']-scale*r['flux1']
for r in rows:
 matches=[z['diff'] for z in rows if z is not r and .5<=z['flux1']/r['flux1']<=2 and abs(z['shape1']-r['shape1'])<=.4]
 r['n_controls']=len(matches);r['control_median']=float(np.median(matches)) if matches else None;r['control_sigma']=float(1.4826*np.median(np.abs(np.array(matches)-np.median(matches)))) if matches else None
 r['passes']=bool(len(matches)>=10 and r['control_sigma']>0 and abs(r['diff']-r['control_median'])>5*r['control_sigma'] and abs(r['diff']-r['control_median'])>3*np.hypot(r['err2'],scale*r['err1'])) if matches else False
with (args.output/'objects.csv').open('w',newline='') as f:
 fields=['x','y','flux1','err1','shape1','flux2','err2','shape2','diff','n_controls','control_median','control_sigma','passes'];w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
summary={'initial_peaks':len(peaks),'selected_centers':len(selected),'paired_positive_flux':len(rows),'photometric_scale':scale,'eligible_ten_controls':sum(r['n_controls']>=10 for r in rows),'threshold_crossings':sum(r['passes'] for r in rows),'note':'Uncalibrated exploratory screen; formal error ignores correlated noise, PSF changes and source association.'}
(args.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
