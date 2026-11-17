"""Export and inspect fixed exploratory diagnostics; needs numpy and astropy."""
import csv,gzip,io,json,math,sys
from pathlib import Path
import numpy as np
from astropy.io import fits
from measure import fisher,verify,EVENT_SHA,GTI_SHA
root=Path(sys.argv[1]); dest=Path(sys.argv[2]);dest.mkdir(exist_ok=True,parents=True)
ep=root/'xa202074010rsl_p0px1000_cl.evt.gz';gp=root/'xa202074010rsl_px1000_exp.gti.gz';verify(ep,EVENT_SHA);verify(gp,GTI_SHA)
with fits.open(ep) as h,fits.open(gp) as gh:
 e=h['EVENTS'].data; g=h['GTI'].data; pg=gh['GTIPIXEL'].data
 fields=['TIME','PIXEL','ITYPE','PI','DETX','DETY']
 with (dest/'events.csv.gz').open('wb') as f:
  with gzip.GzipFile(fileobj=f,mode='wb',filename='',mtime=0) as z:
   with io.TextIOWrapper(z,encoding='ascii',newline='') as text:
    w=csv.writer(text,lineterminator='\n');w.writerow(fields)
    for r in e:w.writerow([format(float(r[k]),'.9f') if k=='TIME' else int(r[k]) for k in fields])
 starts=g['START'];stops=g['STOP'];dt=stops-starts;total=float(dt.sum())
 (dest/'gti.json').write_text(json.dumps(dict(global_gti=[[float(x['START']),float(x['STOP'])] for x in g],pixel_gti=[[int(x['PIXEL']),float(x['START']),float(x['STOP'])] for x in pg]),separators=(',',':'))+'\n')
 def at(frac):
  target=total*frac;i=np.searchsorted(np.cumsum(dt),target);return float(starts[i]+target-dt[:i].sum()) if i<len(starts) else float(stops[-1])
 cut=at(.5);t=e['TIME'];pi=e['PI'];pix=e['PIXEL'];grade=e['ITYPE'];ok=(pix!=12)&(pix!=27)
 pixelok=np.zeros(len(e),bool)
 for p in range(36):
  inds=np.flatnonzero(pix==p); good=np.zeros(len(inds),bool)
  for row in pg[pg['PIXEL']==p]:good|=(t[inds]>=row['START'])&(t[inds]<row['STOP'])
  pixelok[inds]=good
 ok &= pixelok
 def n(m,lo,hi):return int(np.sum(m&(pi>=lo)&(pi<hi)))
 data={'half':[], 'shift':[], 'deciles':[]}
 for i,(lo,hi) in enumerate([(at(0),cut),(cut,at(1))]):
  both=ok&(t>=lo)&(t<hi);hp=both&(grade==0); inner=(e['DETX']>=3)&(e['DETX']<=4)&(e['DETY']>=3)&(e['DETY']<=4)
  data['half'].append(dict(half=i+1,left=n(hp,12400,12540),right=n(hp,12840,12980),low=n(hp,6000,10000),high=n(hp,14000,18000),inner_line=n(hp&inner,12620,12760),outer_line=n(hp&~inner,12620,12760),inner_side=n(hp&inner,12400,12540)+n(hp&inner,12840,12980),outer_side=n(hp&~inner,12400,12540)+n(hp&~inner,12840,12980),inner_hp=int(np.sum(both&inner&(grade==0))),inner_mp=int(np.sum(both&inner&(grade==1))),outer_hp=int(np.sum(both&~inner&(grade==0))),outer_mp=int(np.sum(both&~inner&(grade==1)))))
 for shift in [-10,0,10]:
  cells=[]
  for lo,hi in [(at(0),cut),(cut,at(1))]:
   m=ok&(grade==0)&(t>=lo)&(t<hi);cells.append([n(m,12620+shift,12760+shift),n(m,12400+shift,12540+shift)+n(m,12840+shift,12980+shift)])
  data['shift'].append(dict(shift_eV=shift*.5,counts=cells,p=fisher(*cells[0],*cells[1])))
 for i in range(10):
  lo,hi=at(i/10),at((i+1)/10);m=ok&(grade==0)&(t>=lo)&(t<hi)
  data['deciles'].append(dict(decile=i+1,start=lo,stop=hi,livetime_s=total/10,broad=n(m,6000,18000),line=n(m,12620,12760),side=n(m,12400,12540)+n(m,12840,12980)))
 (dest/'diagnostics.json').write_text(json.dumps(data,indent=2)+'\n');print(json.dumps(data,indent=2))
