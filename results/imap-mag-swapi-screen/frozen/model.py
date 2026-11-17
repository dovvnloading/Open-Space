"""Pure-array, fixed-bin magnetic selection and proton-pressure comparison."""
import numpy as np
BIN=300

def select(t,b,good):
 t=np.asarray(t,float);b=np.asarray(b,float);good=np.asarray(good,bool)
 if b.shape!=(len(t),3) or good.shape!=(len(t),) or not np.all(np.diff(t)>0):raise ValueError('invalid arrays')
 summaries=[];bin_audit=[]
 for idx in range(144):
  take=(t>=idx*BIN)&(t<(idx+1)*BIN);tt=t[take];valid=good[take]&np.isfinite(b[take]).all(axis=1);tv=tt[valid]
  count_ok=len(tv)>=540;edge_ok=bool(len(tv) and tv[0]<=idx*BIN+2 and tv[-1]>=(idx+1)*BIN-2);gap_ok=bool(len(tv)>1 and np.max(np.diff(tv))<=2);usable=count_ok and edge_ok and gap_ok
  bin_audit.append({'bin':idx,'attempted':int(take.sum()),'valid':len(tv),'count_ok':count_ok,'edge_ok':edge_ok,'gap_ok':gap_ok,'usable':usable})
  summaries.append(float(np.median(np.linalg.norm(b[take][valid],axis=1))) if usable else None)
 gates={'triplets_tested':0,'quality_triplets':0,'both_controls_ge_3nT':0,'absolute_drop_ge_1nT':0,'both_controls_fraction_gate':0,'stopped_at_first_selection':False}
 for k in range(1,143):
  gates['triplets_tested']+=1;left,mid,right=summaries[k-1:k+2]
  if None in (left,mid,right):continue
  gates['quality_triplets']+=1;base=(left+right)/2
  if min(left,right)<3:continue
  gates['both_controls_ge_3nT']+=1
  if base-mid<1:continue
  gates['absolute_drop_ge_1nT']+=1
  if mid>.8*min(left,right):continue
  gates['both_controls_fraction_gate']+=1;gates['stopped_at_first_selection']=True
  return {'status':'selected','event_bin':k,'bins':[k-1,k,k+1],'mag_nT':[left,mid,right],'baseline_nT':base,'depletion_nT':base-mid,'valid_bins':sum(x is not None for x in summaries),'valid_triplets':gates['quality_triplets'],'bin_audit':bin_audit,'gates':gates}
 return {'status':'no_candidate' if gates['quality_triplets'] else 'inconclusive_coverage','valid_bins':sum(x is not None for x in summaries),'valid_triplets':gates['quality_triplets'],'bin_audit':bin_audit,'gates':gates}

def pressure_test(times,halfwidth,density,temp,densityerr,temperr,good,pred):
 t=np.asarray(times,float);h=np.asarray(halfwidth,float);n=np.asarray(density,float);T=np.asarray(temp,float);dn=np.asarray(densityerr,float);dT=np.asarray(temperr,float);good=np.asarray(good,bool)
 if any(len(x)!=len(t) for x in (h,n,T,dn,dT,good)):raise ValueError('shape')
 results=[];audit=[]
 for k in pred['bins']:
  contained=np.isfinite(t)&np.isfinite(h)&(h>=0)&(t-h>=k*BIN)&(t+h<=(k+1)*BIN)
  physical=np.isfinite(n)&np.isfinite(T)&np.isfinite(dn)&np.isfinite(dT)&(n>0)&(T>0)&(dn>=0)&(dT>=0)
  use=contained&good&physical
  reason='ok' if use.sum()>=3 else ('fewer_than_3_contained' if contained.sum()<3 else ('fewer_than_3_quality_valid' if (contained&good).sum()<3 else 'fewer_than_3_physical_valid'))
  audit.append({'bin':int(k),'whole_window_contained':int(contained.sum()),'quality_valid':int((contained&good).sum()),'quality_and_physical_valid':int(use.sum()),'reason':reason})
  if use.sum()>=3:
   fac=1.380649e-5;center=fac*n[use]*T[use];low=fac*np.maximum(n[use]-dn[use],0)*np.maximum(T[use]-dT[use],0);high=fac*(n[use]+dn[use])*(T[use]+dT[use]);results.append({'count':int(use.sum()),'p_pPa':float(np.median(center)),'low_scenario_pPa':float(np.median(low)),'high_scenario_pPa':float(np.median(high))})
 if len(results)<3:return {'status':'inconclusive_coverage','bin_audit':audit}
 l,c,r=results;baseline=(l['p_pPa']+r['p_pPa'])/2;controlhigh=(l['high_scenario_pPa']+r['high_scenario_pPa'])/2
 passed=c['p_pPa']>=1.1*baseline and c['low_scenario_pPa']>controlhigh
 return {'status':'conditional_increase' if passed else 'no_conditional_increase','bins':results,'bin_audit':audit,'baseline_pPa':baseline,'robust_to_error_scenario':bool(c['low_scenario_pPa']>controlhigh),'meaning':'finite reported-fit-error scenario; not calibrated significance'}
