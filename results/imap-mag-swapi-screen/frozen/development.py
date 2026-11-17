"""Authorized MAG first-12h reader. SWAPI physical arrays stay sealed."""
import argparse,json
from pathlib import Path
import cdflib,numpy as np
from manifest import RAW,FILES,verify,sha
from model import select
EXPECTED={'date':'2026-07-31','mag_development_utc_start':'2026-07-31T00:00:00','mag_development_utc_end_exclusive':'2026-07-31T12:00:00','bin_seconds':300,'mag_min_good_fraction':0.9,'mag_max_gap_s':2,'mag_min_control_nt':3,'mag_depletion_fraction':0.2,'fraction_reference':'each_control','minimum_field_applies_to':'each_control','mag_depletion_absolute_nt':1,'proton_min_valid_records_per_bin':3,'proton_increase_fraction':0.1,'holdout':'SWAPI physical arrays sealed until prediction and separate authorization'}
def authorize(path,digest,phase,manifest_sha,pred_sha=None):
 if sha(path)!=digest:raise ValueError('authorization digest')
 x=json.loads(Path(path).read_text())
 if x.get('authorized') is not True or x.get('phase')!=phase or x.get('manifest_sha256')!=manifest_sha or (pred_sha is not None and x.get('prediction_sha256')!=pred_sha):raise ValueError('authorization content')
def mag_arrays(c):
 epoch=c.varget('epoch');q=c.varinq('epoch');a=c.varattsget('b_gse');bq=c.varinq('b_gse')
 if q.Data_Type_Description!='CDF_TIME_TT2000' or bq.Dim_Sizes!=[3] or a.get('UNITS')!='nT':raise ValueError('MAG schema')
 start=int(cdflib.cdfepoch.compute_tt2000([2026,7,31,0,0,0,0,0,0]));end=int(cdflib.cdfepoch.compute_tt2000([2026,7,31,12,0,0,0,0,0]));i=int(np.searchsorted(epoch,start));j=int(np.searchsorted(epoch,end));
 if i>=j or np.any(np.diff(epoch)<=0):raise ValueError('missing or unsorted interval')
 b=c.varget('b_gse',startrec=i,endrec=j-1);flags=c.varget('quality_flags',startrec=i,endrec=j-1);bits=c.varget('quality_bitmask',startrec=i,endrec=j-1);t=(epoch[i:j]-start)/1e9;good=(flags==0)&(bits==0)&np.isfinite(b).all(axis=1)&(b!=a['FILLVAL']).all(axis=1)&(b>=a['VALIDMIN']).all(axis=1)&(b<=a['VALIDMAX']).all(axis=1)
 return t,b,good,{'first_record':i,'last_record':j-1,'records':j-i,'good':int(good.sum())}
def main():
 p=argparse.ArgumentParser();p.add_argument('--manifest',required=True);p.add_argument('--manifest-sha256',required=True);p.add_argument('--authorization',required=True);p.add_argument('--authorization-sha256',required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();verify(a.manifest,a.manifest_sha256)
 if json.loads((Path(__file__).parent/'config.json').read_text())!=EXPECTED:raise ValueError('config mismatch')
 authorize(a.authorization,a.authorization_sha256,'mag_development',a.manifest_sha256)
 if a.output.exists() or a.output.resolve().drive.upper()!='D:':raise ValueError('new D output required')
 t,b,good,coverage=mag_arrays(cdflib.CDF(str(RAW/FILES[0])));prediction=select(t,b,good);prediction.update(manifest_sha256=a.manifest_sha256,coverage=coverage,scope='MAG-only first 12 hours; SWAPI sealed')
 with a.output.open('x',encoding='utf-8') as h:json.dump(prediction,h,indent=2,allow_nan=False);h.write('\n')
if __name__=='__main__':main()
