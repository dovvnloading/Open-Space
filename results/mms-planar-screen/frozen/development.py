"""Review-authorized MMS1 only; MMS2-4 physical arrays remain closed."""
from pathlib import Path
import argparse,json
import cdflib,numpy as np
from plane import proposal,predict,sensitivity
from manifest import verify,sha,RAW,META

def authorize(path,expected,manifest_sha,phase,prediction_sha=None):
 if sha(path)!=expected:raise ValueError('authorization digest')
 a=json.loads(Path(path).read_text())
 if a.get('authorized') is not True or a.get('phase')!=phase or a.get('manifest_sha256')!=manifest_sha or (prediction_sha is not None and a.get('prediction_sha256')!=prediction_sha):raise ValueError('authorization content')

EXPECTED={"rotation_deg":45,"minimum_field_nt":3,"mva_half_window_s":5,"mva_eigen_ratio":5,"minimum_normal_speed_km_s":20,"fgm_sensitivity_nt":0.5,"velocity_error_multiplier":1,"time_margin_s":0.3,"max_lag_width_s":10,"heldout_vector_angle_deg":45,"heldout_delta_ratio_min":0.5,"heldout_delta_ratio_max":2,"fixed_utc_start":"2026-09-30T13:28:03.413441756","fixed_utc_end":"2026-09-30T13:30:12.157399823","quality":"FGM and DIS flags zero; nearest DIS within 0.08 s","source":"MMS1","heldout":["MMS2","MMS3","MMS4"]}
def checked(c,name,unit,dims):
 q=c.varinq(name);a=c.varattsget(name)
 if a.get('UNITS')!=unit or list(q.Dim_Sizes)!=dims or 'FILLVAL' not in a or 'VALIDMIN' not in a or 'VALIDMAX' not in a:raise ValueError('invalid CDF schema '+name)
 data=c.varget(name);low=np.asarray(a['VALIDMIN']);high=np.asarray(a['VALIDMAX']);good=np.isfinite(data)&(data!=a['FILLVAL'])&(data>=low)&(data<=high)
 return data,good

def read_sc1(f,d,cfg,coverage=None):
 ft=f.varget('Epoch');dt=d.varget('Epoch')
 if f.varinq('Epoch').Data_Type_Description!='CDF_TIME_TT2000' or d.varinq('Epoch').Data_Type_Description!='CDF_TIME_TT2000' or len(ft)==0 or len(dt)==0 or np.any(np.diff(ft)<=0) or np.any(np.diff(dt)<=0):raise ValueError('invalid TT2000 axis')
 def ns(x):
  import re
  z=[int(q) for q in re.findall(r'\d+',x)];z=z[:6]+[int(str(z[6]).ljust(9,'0')[:3]),int(str(z[6]).ljust(9,'0')[3:6]),int(str(z[6]).ljust(9,'0')[6:9])];return int(cdflib.cdfepoch.compute_tt2000(z))
 low,high=ns(cfg['fixed_utc_start']),ns(cfg['fixed_utc_end']);mask=(ft>=low)&(ft<=high)
 if coverage is not None:coverage.update(native_fgm_rows=len(ft),native_dis_rows=len(dt),cropped_fgm_rows=int(mask.sum()),cropped_first_fgm_native_index=int(np.flatnonzero(mask)[0]) if mask.any() else None,cropped_dis_rows=int(((dt>=low)&(dt<=high)).sum()))
 if mask.sum()<2:raise ValueError('source interval empty')
 allb,bvalid=checked(f,'mms1_fgm_b_gse_brst_l2','nT',[4]);allv,vvalid=checked(d,'mms1_dis_bulkv_gse_brst','km/s',[3]);alle,evalid=checked(d,'mms1_dis_bulkv_err_brst','km/s',[3]);fg=f.varget('mms1_fgm_flag_brst_l2');dg=d.varget('mms1_dis_errorflags_brst')
 dmask=(dt>=low)&(dt<=high)
 if not dmask.any():raise ValueError('DIS interval empty')
 dt=dt[dmask];allv=allv[dmask];vvalid=vvalid[dmask];alle=alle[dmask];evalid=evalid[dmask];dg=dg[dmask]
 ft=ft[mask];b=allb[mask,:3];bv=bvalid[mask,:3].all(axis=1);fg=fg[mask]==0
 ix=np.searchsorted(dt,ft);r=np.clip(ix,0,len(dt)-1);l=np.clip(ix-1,0,len(dt)-1);p=np.where(abs(dt[l]-ft)<=abs(dt[r]-ft),l,r)
 good=fg&bv&(dg[p]==0)&vvalid[p].all(axis=1)&evalid[p].all(axis=1)&(np.abs(dt[p]-ft)<=80000000)
 if coverage is not None:coverage.update(good_aligned_fgm_rows=int(good.sum()),nearest_dis_within_tolerance=int((np.abs(dt[p]-ft)<=80000000).sum()))
 origin=int(ft[0]);return origin,(ft-origin)/1e9,b,allv[p],alle[p],good

def empty_classification(gates):
 return 'inconclusive_no_valid_coverage' if gates['quality_cadence_windows']==0 else ('inconclusive_model_gates' if gates['rotation_windows'] else 'no_qualified_candidate')

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--manifest',required=True);ap.add_argument('--manifest-sha256',required=True);ap.add_argument('--authorization',required=True);ap.add_argument('--authorization-sha256',required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
 m=verify(a.manifest,a.manifest_sha256);cfg=json.loads((Path(__file__).parent/'config.json').read_text());
 if cfg!=EXPECTED:raise ValueError('configuration differs from implemented gates')
 authorize(a.authorization,a.authorization_sha256,a.manifest_sha256,'sc1_development')
 if a.output.resolve().drive.upper()!='D:' or a.output.exists():raise ValueError('new D output required')
 fs=sorted(n for n in m['inputs'] if n.startswith('mms1_fgm_'));ds=sorted(n for n in m['inputs'] if n.startswith('mms1_fpi_'))
 if len(fs)!=1 or len(ds)!=1:raise ValueError('source identity')
 coverage={};f=cdflib.CDF(str(RAW/fs[0]));d=cdflib.CDF(str(RAW/ds[0]));origin,t,field,velocity,err,good=read_sc1(f,d,cfg,coverage)
 gates={};found=proposal(t,field,velocity,good,gates)
 empty_status=empty_classification(gates)
 out={'manifest_sha256':a.manifest_sha256,'git_head':m['git_head'],'fgm_filename':fs[0],'dis_filename':ds[0],'origin_tt2000':origin,'selected':found,'coverage':coverage,'gates':gates,'sensitivity_failures':0,'status':empty_status if found is None else 'inconclusive','scope':'SC1 only; engineering sensitivity not confidence'}
 if found:
  found['native_fgm_index']=coverage['cropped_first_fgm_native_index']+found['index']
  when=origin+round(found['time']*1e9);pos=[]
  for craft in range(1,5):
   name=next(n for n in m['inputs'] if n.startswith(f'mms{craft}_fgm_'));c=cdflib.CDF(str(RAW/name));pt=c.varget('Epoch_state');xyz=c.varget(f'mms{craft}_fgm_r_gse_brst_l2');
   if not pt[0]<=when<=pt[-1]:raise ValueError('position outside coverage')
   pos.append([float(np.interp(0,pt.astype(float)-when,xyz[:,j])) for j in range(3)])
  s=sensitivity(t,field,velocity,err,good,found,pos);out.update(positions_gse_km=pos,nominal_lags_s=predict(found,pos),sensitivity=s,status=s['status'])
  if s['status']!='qualified':out['sensitivity_failures']=1
  if s['status']=='qualified':
   meta=json.loads(META.read_text());segments=meta['joint_segments'];start=when+round(min(s['lag_low_s'])*1e9)-5_000_000_000;end=when+round(max(s['lag_high_s'])*1e9)+6_000_000_000
   if not any(q['start_ns']<=start and q['end_ns']>=end for q in segments):out['status']='inconclusive_joint_quality_coverage'
 with a.output.open('x',encoding='utf-8') as h:json.dump(out,h,indent=2,allow_nan=False);h.write('\n')
if __name__=='__main__':main()
