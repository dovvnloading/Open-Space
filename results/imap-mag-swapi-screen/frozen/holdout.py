"""Separately authorized SWAPI proton test; reads only frozen three bins."""
import argparse,json
from pathlib import Path
import cdflib,numpy as np
from manifest import RAW,FILES,verify,sha
from development import authorize,EXPECTED
from model import pressure_test

def read_selected(c,pred):
 ep=c.varinq('epoch');hs=c.varinq('epoch_delta');ha=c.varattsget('epoch_delta')
 if ep.Data_Type_Description!='CDF_TIME_TT2000' or ep.Dim_Sizes!=[] or hs.Dim_Sizes!=[] or hs.Data_Type_Description!='CDF_INT8' or ha.get('UNITS')!='ns' or any(k not in ha for k in ('FILLVAL','VALIDMIN','VALIDMAX')):raise ValueError('epoch_delta schema')
 epoch=c.varget('epoch');h=c.varget('epoch_delta');flags=c.varget('swp_flags');start=int(cdflib.cdfepoch.compute_tt2000([2026,7,31,0,0,0,0,0,0]));t=(epoch-start)/1e9;half=h/1e9
 if c.varinq('epoch').Data_Type_Description!='CDF_TIME_TT2000' or np.any(np.diff(epoch)<=0):raise ValueError('epoch schema')
 mask=np.zeros(len(t),bool)
 hvalid=(h!=ha['FILLVAL'])&(h>=ha['VALIDMIN'])&(h<=ha['VALIDMAX'])
 for k in pred['bins']:mask|=(t-half>=k*300)&(t+half<=(k+1)*300)&hvalid
 idx=np.flatnonzero(mask)
 if not len(idx):result={'status':'inconclusive_coverage','bin_audit':[{'bin':int(k),'whole_window_contained':0,'quality_valid':0,'quality_and_physical_valid':0,'reason':'fewer_than_3_contained'} for k in pred['bins']]}
 else:
  names=['proton_sw_density','proton_sw_temperature','proton_sw_density_uncert','proton_sw_temperature_uncert'];units=['#/cm^3','K','#/cm^3','K'];vals=[];valid=np.ones(len(idx),bool)
  for name,unit in zip(names,units):
   spec=c.varinq(name);attr=c.varattsget(name)
   if spec.Dim_Sizes!=[] or attr.get('UNITS')!=unit or any(k not in attr for k in ('FILLVAL','VALIDMIN','VALIDMAX')):raise ValueError('SWAPI schema '+name)
   v=np.asarray([c.varget(name,startrec=int(i),endrec=int(i)).reshape(-1)[0] for i in idx]);vals.append(v);valid &= np.isfinite(v)&(v!=attr['FILLVAL'])&(v>=attr['VALIDMIN'])&(v<=attr['VALIDMAX'])
  for v in vals:v[~valid]=np.nan
  result=pressure_test(t[idx],half[idx],*vals,flags[idx]==0,pred)
 result['time_audit']={'native_records':len(epoch),'invalid_epoch_delta_records':int((~hvalid).sum()),'selected_contained_records':len(idx)}
 return result,idx

def main():
 p=argparse.ArgumentParser();p.add_argument('--manifest',required=True);p.add_argument('--manifest-sha256',required=True);p.add_argument('--authorization',required=True);p.add_argument('--authorization-sha256',required=True);p.add_argument('--prediction',required=True);p.add_argument('--prediction-sha256',required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();verify(a.manifest,a.manifest_sha256)
 if json.loads((Path(__file__).parent/'config.json').read_text())!=EXPECTED:raise ValueError('config mismatch')
 if sha(a.prediction)!=a.prediction_sha256:raise ValueError('prediction digest')
 pred=json.loads(Path(a.prediction).read_text());authorize(a.authorization,a.authorization_sha256,'swapi_holdout',a.manifest_sha256,a.prediction_sha256)
 if pred['manifest_sha256']!=a.manifest_sha256 or pred['status']!='selected':raise ValueError('no selected prediction')
 if a.output.exists() or a.output.resolve().drive.upper()!='D:':raise ValueError('new D output required')
 result,idx=read_selected(cdflib.CDF(str(RAW/FILES[1])),pred)
 result.update(prediction_sha256=a.prediction_sha256,selected_record_indices=idx.tolist(),scope='conditional proton pressure; not total pressure or calibrated significance')
 with a.output.open('x',encoding='utf-8') as h:json.dump(result,h,indent=2,allow_nan=False);h.write('\n')
if __name__=='__main__':main()
