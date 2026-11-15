import numpy as np
from scipy.optimize import minimize_scalar
from pathlib import Path
def rot(t,y,p,n=3):
 X=np.c_[np.ones(len(t)),*[f(2*np.pi*t*k/p) for k in range(1,n+1) for f in [np.cos,np.sin]]]
 return y-X@np.linalg.lstsq(X,y,rcond=None)[0]
def depth(t,r,p,w,c):
 ph=t/p-c;n=ph.round();d=abs(ph-n)
 rows=[]
 for k in np.unique(n):
  a=(n==k)&(d<w/2/p);b=(n==k)&(d>.27)
  if a.sum()>=2 and b.sum()>=4:
   rows.append([k,np.median(r[a])-np.median(r[b])])
 return np.array(rows)
def evaluate(z,id):
 ps,pe,w={'29489':(3.1679,12.35327485665821,2.0074),'3177':(3.547245,18.559880330703862,2.32)}[id]
 t=(z[:,0]-z[0,0])*24
 for j in (8,9):
  flux=z[:,j];m=(flux>0)&np.isfinite(flux)&(z[:,7]<np.quantile(z[:,7],.97));t2=t[m];y=-2.5*np.log10(flux[m])
  med=np.array([np.median(y[abs(t2-tt)<=15]) for tt in t2]);y-=med
  r=rot(t2,y,ps);sigma=np.median(abs(r-np.median(r)))*1.4826;m=abs(r-np.median(r))<4*sigma;t2=t2[m];y=y[m]
  s=minimize_scalar(lambda p:np.mean(rot(t2,y,p)**2),bounds=(.997*ps,1.003*ps),method='bounded').x;r=rot(t2,y,s)
  tcut=np.median(t2);first=t2<tcut
  print('extract',j,'p',round(s,4),'MAD',round(np.median(abs(r-np.median(r)))*1.4826,3),'rn',round(np.std(rot(t2,y,ps))/np.std(y),2), 'n',len(t2))
  for key,sel in [('all',np.ones(len(t2),bool)),('first',first),('second',~first)]:
   c=np.arange(100)/100
   def fit(cc):
    ph=abs((t2[sel]/pe-cc+.5)%1-.5);return np.mean(r[sel][ph<w/pe/2])-np.mean(r[sel][ph>.27])
   ds=[fit(cc) for cc in c];cb=c[np.argmax(ds)]
   if key=='first':cf=cb
   cc=cf if key=='second' else cb
   ev=depth(t2[sel],r[sel],pe,w,cc);print(key,'phase fit',cb,'mean fit',round(np.max(ds),3),'eval phase',cc,'events',len(ev),int(sum(ev[:,1]>0)),round(np.median(ev[:,1]),4))
if __name__=='__main__':
 for f in sorted(Path('tessex/independent').glob('stamps*.npz')):
  print(f.name);evaluate(np.load(f)['rows'],f.name[6:].split('_')[0])
