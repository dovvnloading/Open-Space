"""Frozen V6 development test; centers come exclusively from first visit."""
import argparse,csv,json
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
p=argparse.ArgumentParser();p.add_argument('--stamps',required=True,type=Path);p.add_argument('--centers',required=True,type=Path);p.add_argument('--output',required=True,type=Path);args=p.parse_args();args.output.mkdir(parents=True,exist_ok=True)
z=np.load(args.stamps);data=[];allow=(1<<18)|(1<<24)
for i in (1,2):
 a,e,f=[z[f'{k}{i}'] for k in ('science','rms','flag')];data.append((a,e,np.isfinite(a)&np.isfinite(e)&(e>0)&((f.astype(np.int64)&~allow)==0)))
centers=[(int(r['x']),int(r['y'])) for r in csv.DictReader(args.centers.open())];yy,xx=np.mgrid[-8:9,-8:9];disk=xx**2+yy**2<=16;ring=(xx**2+yy**2>=36)&(xx**2+yy**2<=64);rows=[];counts={'centers':len(centers),'first_invalid':0,'first_nonpositive':0,'first_elongated':0,'compact_first':0,'second_invalid':0,'second_nonpositive':0}
def measure(x,y,entry):
 a,e,good=entry;sl=np.s_[y-8:y+9,x-8:x+9];m=good[sl]
 if not m[disk].all() or m[ring].sum()<50:return None
 im=a[sl];sky=np.median(im[ring&m]);v=im[disk]-sky;rr=e[sl];flux=float(np.sum(v));n=m[ring].sum();var=float(np.sum(rr[disk]**2)+len(v)**2*np.pi/2*np.sum(rr[ring&m]**2)/n**2);w=np.maximum(v,0)
 if w.sum()==0:return flux,np.sqrt(var),None,None
 x0=np.sum(w*xx[disk])/w.sum();y0=np.sum(w*yy[disk])/w.sum();vx=np.sum(w*(xx[disk]-x0)**2)/w.sum();vy=np.sum(w*(yy[disk]-y0)**2)/w.sum();xy=np.sum(w*(xx[disk]-x0)*(yy[disk]-y0))/w.sum();eig=np.linalg.eigvalsh([[vx,xy],[xy,vy]]);ratio=float(eig[1]/eig[0]) if eig[0]>0 else float('inf');return flux,np.sqrt(var),float(np.sqrt(vx+vy)),ratio
for idx,(x,y) in enumerate(centers):
 a=measure(8,8,tuple(a[idx] for a in data[0]))
 if a is None:counts['first_invalid']+=1;continue
 if a[0]<=0 or a[2] is None:counts['first_nonpositive']+=1;continue
 if a[3]>3:counts['first_elongated']+=1;continue
 counts['compact_first']+=1;b=measure(8,8,tuple(a[idx] for a in data[1]))
 if b is None:counts['second_invalid']+=1;continue
 if b[0]<=0:counts['second_nonpositive']+=1
 rows.append({'x':x,'y':y,'flux1':a[0],'err1':a[1],'shape1':a[2],'axis_ratio2':a[3],'flux2':b[0],'err2':b[1]})
for r in rows:
 controls=[z for z in rows if z is not r and .5<=z['flux1']/r['flux1']<=2 and abs(z['shape1']-r['shape1'])<=.4];r['n_controls']=len(controls);r['passes']=False
 if len(controls)>=10:
  scale=float(np.median([z['flux2']/z['flux1'] for z in controls]));dif=np.array([z['flux2']-scale*z['flux1'] for z in controls]);center=float(np.median(dif));sigma=float(1.4826*np.median(np.abs(dif-center)));res=r['flux2']-scale*r['flux1']-center;formal=float(np.hypot(r['err2'],scale*r['err1']));r.update(scale=scale,residual=res,control_sigma=sigma,formal_error=formal);r['passes']=bool(sigma>0 and abs(res)>5*sigma and abs(res)>3*formal and abs(res)>3*0.044814546038339294)
fields=['x','y','flux1','err1','shape1','axis_ratio2','flux2','err2','n_controls','scale','residual','control_sigma','formal_error','passes']
with (args.output/'objects.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
counts.update(paired=len(rows),eligible_ten_controls=sum(r['n_controls']>=10 for r in rows),threshold_crossings=sum(r['passes'] for r in rows));(args.output/'summary.json').write_text(json.dumps(counts,indent=2)+'\n');print(json.dumps(counts,indent=2))
