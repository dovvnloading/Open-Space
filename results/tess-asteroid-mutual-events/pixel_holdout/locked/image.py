from pathlib import Path
import numpy as np
from scipy.ndimage import gaussian_filter,median_filter
from scipy.stats import spearmanr
from numpy.lib.stride_tricks import sliding_window_view
def get(d):
 a=[]
 for f in sorted(Path(d).glob('*.npz')):
  try:
   z=np.load(f);a.append((float(z['t']),float(z['x']),float(z['y']),int(z['q']),z['flux'],f.name))
  except (KeyError,ValueError,OSError):continue
 return a
def align(a,b,R=4):
 # Absolute offset of second from first. A star moves oppositely on the stamp.
 best=(1e6,0,0)
 for dx in range(-R,R+1):
  for dy in range(-R,R+1):
   L=a[max(0,dy):min(13,13+dy),max(0,dx):min(13,13+dx)]; H=b[max(0,-dy):min(13,13-dy),max(0,-dx):min(13,13-dx)]
   d=L-H;d=d[np.isfinite(d)]
   if len(d)<45:continue
   cost=np.median(np.abs(d-np.median(d)))
   if cost<best[0]:best=(cost,dx,dy)
 return best
def reduce(frames):
 n=len(frames);offset=np.zeros((n,2),int);seg=np.zeros(n,int);scores=np.zeros(n)
 for i in range(1,n):
  if frames[i][0]-frames[i-1][0]>.065:seg[i]=seg[i-1]+1;continue
  sc,dx,dy=align(frames[i-1][4],frames[i][4]);scores[i]=sc
  offset[i]=offset[i-1]+[dx,dy];seg[i]=seg[i-1]
 out=[];y,x=np.mgrid[-6:7,-6:7];dif=[]
 for i,(t,fx,fy,q,arr,name) in enumerate(frames):
  if q:continue
  dist=np.linalg.norm(offset-offset[i],axis=1)
  cand=np.where((seg==seg[i])&(dist>=2.5)&(dist<12)&(np.arange(n)!=i))[0]
  refs=[]
  for j in cand:
   _,sx,sy,qj,im,_=frames[j]
   if qj:continue
   dx,dy=offset[j]-offset[i];rx=x-dx;ry=y-dy
   valid=(abs(rx)<=6)&(abs(ry)<=6)&(np.hypot(rx-sx,ry-sy)>2.6)
   m=np.full((13,13),np.nan);m[valid]=im[(ry[valid]+6).astype(int),(rx[valid]+6).astype(int)]
   comp=valid&(np.hypot(x-fx,y-fy)>2.6)
   if comp.sum()<10:continue
   m-=np.nanmedian(m[comp]-arr[comp])
   refs.append(m)
  if len(refs)<3:continue
  stack=np.array(refs)
  import warnings
  with warnings.catch_warnings():
   warnings.simplefilter('ignore'); ref=np.nanmedian(stack,axis=0)
  diff=arr-ref
  left=(np.isfinite(diff))&(np.hypot(x-fx,y-fy)>3)
  if left.sum()<8:continue
  back=np.median(diff[left]);diff-=back
  for cut in [0]:pass
  r=np.hypot(x-fx,y-fy); aper=r<1.65; prf=np.exp(-.5*r*r/.8**2);v=(r<2.8)&np.isfinite(diff)
  if sum(np.isfinite(diff[aper]))<aper.sum() or v.sum()<12:continue
  fit=diff[v]@prf[v]/(prf[v]@prf[v]);noise=1.4826*np.nanmedian(abs(diff[left]-np.median(diff[left])))
  out.append((t,fx,fy,offset[i,0],offset[i,1],seg[i],q,noise,fit,diff[aper].sum(),np.sum(refs[0]==0),len(refs),np.nanmedian(arr)))
  dif.append(diff)
 return np.array(out),np.array(dif),seg,offset,scores
if __name__=='__main__':
 import sys
 for d in sys.argv[1:]:
  a=get(d); out,dif,seg,offs,s=reduce(a); print(d,'frames',len(a),'stats',len(out),'segs',np.unique(seg,return_counts=True),'shifts',np.unique(np.diff(offs,axis=0),axis=0,return_counts=True)[1].max(),np.percentile(s,[0,50,80,99]),'fit,scatter,refs',np.nanpercentile(out[:,[7,8,11]],[5,50,95],axis=0),flush=True)
  np.savez_compressed('tessex/independent/'+Path(d).name+'.npz',rows=out,diffs=dif)
