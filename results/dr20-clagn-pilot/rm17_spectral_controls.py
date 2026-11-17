"""RM17 development-only [O III] doublet and observed-frame color controls."""
import csv,json,hashlib,sys
from pathlib import Path
import numpy as np
from astropy.io import fits
src=Path(sys.argv[1]);rec=json.loads(Path(sys.argv[2]).read_text());items=[r for r in rec if r['sdss_id']==62239511 and r['kind']!='hold'];z=next(r['z'] for r in items if r['kind']=='baseline');spectra={}
for r in items:
 p=src/r['filename'];assert hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256']
 with fits.open(p) as h:
  c=h['COADD'].data;obs=10**c['loglam'];flux=np.array(c['flux'],float);iv=np.array(c['ivar'],float);good=(iv>0)&np.isfinite(flux)&(c['and_mask']==0);spectra[r['kind']]=(obs,flux,iv,good)
def integrate(kind,core,left,right,offset):
 obs,flux,iv,good=spectra[kind];w=obs/(1+z)
 def window(a,b):
  all=(w>=a+offset)&(w<b+offset);sel=all&good;return sel, (sel.sum()/all.sum() if all.sum() else 0)
 # Core is held fixed; only continuum windows move.
 c,fc=window(*core);c=(w>=core[0])&(w<core[1])&good;full=(w>=core[0])&(w<core[1]);fc=c.sum()/full.sum() if full.sum() else 0
 a,fa=window(*left);b,fb=window(*right)
 if min(fc,fa,fb)<.8:return dict(coverage=False)
 xa=np.mean(w[a]);xb=np.mean(w[b]);ya=np.mean(flux[a]);yb=np.mean(flux[b]);t=(w[c]-xa)/(xb-xa);dw=np.gradient(w)[c];line=np.sum((flux[c]-(1-t)*ya-t*yb)*dw)
 variance=np.sum(dw**2/iv[c])+(np.sum(dw*(1-t))**2)*np.sum(1/iv[a])/a.sum()**2+(np.sum(dw*t)**2)*np.sum(1/iv[b])/b.sum()**2
 return dict(coverage=True,flux=float(line),formal_sigma=float(1.2*np.sqrt(variance)),ncore=int(c.sum()))
lines={}
for kind in spectra:
 lines[kind]={}
 for offset in (-5,0,5):
  lines[kind][str(offset)]={'4960':integrate(kind,(4950,4970),(4930,4943),(4975,4988),offset),'5008':integrate(kind,(4998,5018),(4975,4988),(5025,5038),offset)}
qa=[r for r in csv.DictReader(open(Path(__file__).with_name('development-qa.csv'),newline='')) if r['sdss_id']=='62239511'];scale=float(next(r['oiii_flux'] for r in qa if r['kind']=='baseline'))/float(next(r['oiii_flux'] for r in qa if r['kind']=='dev'))
ob,fb,ib,gb=spectra['baseline'];od,fd,idv,gd=spectra['dev'];colors=[]
for center in (5500,6500,7500,8500):
 def med(w,f,g):
  region=(w>=center-50)&(w<center+50);m=region&g;return (float(np.median(f[m])) if m.sum() and m.sum()/region.sum()>=.8 else None, float(m.sum()/region.sum()) if region.sum() else 0)
 x,cx=med(ob,fb,gb);y,cy=med(od,fd,gd);colors.append(dict(center=center,baseline=x,dev=y,baseline_coverage=cx,dev_coverage=cy,scaled_ratio=(scale*y/x if x and x>0 and y is not None else None)))
result=dict(sdss_id=62239511,oiii=lines,scale=scale,colors=colors);Path(sys.argv[3]).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
