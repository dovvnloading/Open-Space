"""Conditional planar-screen primitives; no implicit access to observations."""
import numpy as np
BOUND_NT=.5
MIN_SPEED=20.

def proposal(t,b,v,good,audit=None):
 t=np.asarray(t,float);b=np.asarray(b,float);v=np.asarray(v,float);good=np.asarray(good,bool)
 if b.shape!=(len(t),3) or v.shape!=b.shape or good.shape!=(len(t),) or not np.all(np.diff(t)>0):raise ValueError('shape/time')
 counts={'tested_1s_windows':0,'quality_cadence_windows':0,'field_endpoint_windows':0,'rotation_windows':0,'mva_failures':0,'speed_failures':0,'scanned_to_first_selection':False}
 if audit is not None:audit.update(counts);counts=audit
 for i in range(len(t)):
  j=np.searchsorted(t,t[i]+1.)
  if j>=len(t) or t[j]-t[i]>1.02:continue
  counts['tested_1s_windows']+=1
  k=(t>=t[i]-5)&(t<=t[i]+5)
  if t[i]-t[0]<5 or t[-1]-t[i]<5 or np.count_nonzero(k)<1000 or not good[k].all() or np.max(np.diff(t[k]))>.02:continue
  if not np.isfinite(b[k]).all() or not np.isfinite(v[k]).all():continue
  counts['quality_cadence_windows']+=1
  a=b[i];c=b[j];na=np.linalg.norm(a);nc=np.linalg.norm(c)
  if min(na,nc)<3:continue
  counts['field_endpoint_windows']+=1
  theta=np.degrees(np.arccos(np.clip(a@c/(na*nc),-1,1)))
  if theta<45:continue
  counts['rotation_windows']+=1
  e,u=np.linalg.eigh(np.cov(b[k].T))
  if e[0]<=0 or e[1]/e[0]<5:
   counts['mva_failures']+=1;continue
  n=u[:,0];vel=np.median(v[k],axis=0);vn=float(vel@n)
  if abs(vn)<MIN_SPEED:
   counts['speed_failures']+=1;continue
  counts['scanned_to_first_selection']=True
  return {'index':i,'time':float(t[i]),'normal':n.tolist(),'normal_speed_km_s':vn,'rotation_deg':float(theta),'mva_eigenvalues':e.tolist(),'B_before':a.tolist(),'B_after':c.tolist()}
 return None

def predict(proposal,positions):
 p=np.asarray(positions,float);n=np.asarray(proposal['normal']);speed=proposal['normal_speed_km_s']
 if p.shape!=(4,3) or not np.isfinite(p).all() or abs(speed)<MIN_SPEED:raise ValueError('invalid geometry')
 return ((p[1:]-p[0])@n/speed).tolist()

def sensitivity(t,b,v,verr,good,found,positions):
 """Finite declared perturbation scenarios, not probabilistic or certified bounds."""
 t=np.asarray(t);b=np.asarray(b);v=np.asarray(v);verr=np.asarray(verr);p=np.asarray(positions);i=found['index'];k=(t>=t[i]-5)&(t<=t[i]+5)
 if verr.shape!=v.shape or not np.isfinite(verr[k]).all() or np.min(verr[k])<0:return {'status':'inconclusive_velocity_errors'}
 base=np.median(v[k],axis=0);err=np.max(verr[k],axis=0);sign=np.where(t[k]<t[i],-1.,1.);lags=[]
 for axis in range(4):
  for sg in ((0,) if axis==3 else (-1,1)):
   bb=b[k].copy()
   if axis<3:bb[:,axis]+=sg*BOUND_NT*sign
   e,u=np.linalg.eigh(np.cov(bb.T))
   if e[0]<=0 or e[1]/e[0]<5:return {'status':'inconclusive_mva_sensitivity'}
   n=u[:,0]
   for mask in range(8):
    vv=base+err*np.array([1 if mask&(1<<j) else -1 for j in range(3)]);den=float(vv@n)
    if abs(den)<MIN_SPEED or den*float(base@n)<=0:return {'status':'inconclusive_normal_speed'}
    lags.append((p[1:]-p[0])@n/den)
 q=np.asarray(lags);lo=q.min(axis=0)-.3;hi=q.max(axis=0)+.3
 if not np.isfinite(q).all() or np.any(hi-lo>10):return {'status':'inconclusive_broad_lag'}
 return {'status':'qualified','lag_low_s':lo.tolist(),'lag_high_s':hi.tolist(),'scenarios':len(q),'meaning':'finite sensitivity scenarios only; not calibrated confidence bounds'}

def compare(t,b,good,source,lo,hi):
 """Search predefined 1s rotations within each prediction interval, preserving zero-lag alternative."""
 t=np.asarray(t,float);b=np.asarray(b,float);good=np.asarray(good,bool)
 if b.shape!=(len(t),3) or good.shape!=(len(t),) or not np.all(np.diff(t)>0):return {'status':'inconclusive_invalid_data'}
 center=source['time'];a=np.asarray(source['B_before']);z=np.asarray(source['B_after']);ref=z-a;nr=np.linalg.norm(ref)
 if not np.isfinite(ref).all() or nr<=0:return {'status':'inconclusive_source'}
 def scan(l,h):
  candidates=[];evaluated=0
  for i in np.flatnonzero((t>=center+l)&(t<=center+h)):
   j=np.searchsorted(t,t[i]+1)
   if j>=len(t) or t[j]-t[i]>1.02 or not good[i:j+1].all() or not np.isfinite(b[i:j+1]).all() or np.max(np.diff(t[i:j+1]))>.02:continue
   evaluated+=1
   delta=b[j]-b[i];nd=np.linalg.norm(delta)
   if min(np.linalg.norm(b[i]),np.linalg.norm(b[j]))<3 or nd==0:continue
   angle=float(np.degrees(np.arccos(np.clip(delta@ref/(nd*nr),-1,1))))
   if angle<=45 and .5<=nd/nr<=2:candidates.append(float(t[i]-center))
  return candidates,evaluated
 if center+lo<t[0] or center+hi+1>t[-1] or lo>hi:return {'status':'inconclusive_coverage'}
 match,n=scan(lo,hi);zero,nzero=scan(-.3,.3)
 if n==0 or nzero==0:return {'status':'inconclusive_quality'}
 return {'status':'matched' if match else 'no_match','matched_lags_s':match,'zero_lag_lags_s':zero,'zero_lag_alternative':bool(zero)}
