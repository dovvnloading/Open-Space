import unittest,tempfile,sys
from pathlib import Path
from unittest.mock import patch
import numpy as np,cdflib
from model import select,pressure_test
import development,holdout
class Tests(unittest.TestCase):
 def test_depletion_and_pressure(self):
  t=np.arange(0,43200,.5);b=np.zeros((len(t),3));b[:,0]=6;b[(t>=300)&(t<600),0]=4;g=np.ones(len(t),bool);p=select(t,b,g);self.assertEqual(p['bins'],[0,1,2]);tt=np.array([60,120,180,360,420,480,660,720,780]);n=np.array([5,5,5,8,8,8,5,5,5]);T=np.ones(9)*100000;out=pressure_test(tt,np.ones(9)*10,n,T,np.ones(9)*.1,np.ones(9)*1000,np.ones(9,bool),p);self.assertEqual(out['status'],'conditional_increase');n[3:6]=5;self.assertEqual(pressure_test(tt,np.ones(9)*10,n,T,np.ones(9)*.1,np.ones(9)*1000,np.ones(9,bool),p)['status'],'no_conditional_increase')
 def test_no_coverage_and_gap(self):
  t=np.arange(0,900,.5);b=np.ones((len(t),3));g=np.zeros(len(t),bool);self.assertEqual(select(t,b,g)['status'],'inconclusive_coverage')
 def test_native_mag_cdf_scope_flags(self):
  with tempfile.TemporaryDirectory(dir='D:/Open-Space/campaigns') as tmp:
   p=Path(tmp)/'synthetic.cdf';w=cdflib.cdfwrite.CDF(p,delete=True);start=cdflib.cdfepoch.compute_tt2000([2026,7,31,0,0,0,0,0,0]);times=start+np.arange(86410,dtype=np.int64)*500000000;b=np.tile([6.,0,0],(len(times),1));b[700:1200,0]=4;b[100,0]=-1e31;flag=np.zeros(len(times),np.uint8);flag[101]=1;bits=np.zeros(len(times),np.uint16);bits[102]=4
   for name,data,kind,dims,attrs in [('epoch',times,33,[],{'UNITS':'ns'}),('b_gse',b,45,[3],{'UNITS':'nT','FILLVAL':-1e31,'VALIDMIN':-100000.,'VALIDMAX':100000.}),('quality_flags',flag,11,[],{}),('quality_bitmask',bits,12,[],{})]:w.write_var({'Variable':name,'Data_Type':kind,'Num_Elements':1,'Rec_Vary':True,'Dim_Sizes':dims},var_attrs=attrs,var_data=data)
   w.close();c=cdflib.CDF(p);t,b,good,counts=development.mag_arrays(c);self.assertEqual(counts['records'],86400);self.assertEqual(counts['good'],86397);self.assertLess(t[-1],43200);self.assertEqual(select(t,b,good)['bins'],[0,1,2]);del c
   import gc;gc.collect()
 def test_authorization_before_cdf(self):
  args=['x','--manifest','x','--manifest-sha256','x','--authorization','x','--authorization-sha256','x','--output','D:/x']
  with patch.object(sys,'argv',args),patch.object(development,'verify'),patch.object(development,'authorize',side_effect=ValueError('sealed')),patch.object(development.cdflib,'CDF',side_effect=AssertionError('opened')):
   with self.assertRaisesRegex(ValueError,'sealed'):development.main()
 def test_native_swapi_selected_bins_only(self):
  with tempfile.TemporaryDirectory(dir='D:/Open-Space/campaigns') as tmp:
   p=Path(tmp)/'synthetic-swapi.cdf';w=cdflib.cdfwrite.CDF(p,delete=True);start=cdflib.cdfepoch.compute_tt2000([2026,7,31,0,0,0,0,0,0]);secs=np.arange(30,1530,60,dtype=np.int64);epoch=start+secs*1000000000
   for name,data,kind,unit,lo,hi in [('epoch',epoch,33,'ns',0,2**63-1),('epoch_delta',np.ones(len(secs),np.int64)*20000000000,8,'ns',0,86000000000000),('swp_flags',np.zeros(len(secs),np.uint16),12,' ',0,65534),('proton_sw_density',np.where((secs>=300)&(secs<600),8.,5.),45,'#/cm^3',0.,150.),('proton_sw_temperature',np.ones(len(secs))*1e5,45,'K',0.,7e6),('proton_sw_density_uncert',np.ones(len(secs))*.1,45,'#/cm^3',0.,150.),('proton_sw_temperature_uncert',np.ones(len(secs))*1000,45,'K',0.,7e6)]:
    attr={'UNITS':unit,'FILLVAL':-1e31,'VALIDMIN':lo,'VALIDMAX':hi} if kind==45 else ({'UNITS':unit,'FILLVAL':-9223372036854775808,'VALIDMIN':0,'VALIDMAX':86000000000000} if name=='epoch_delta' else {'UNITS':unit})
    w.write_var({'Variable':name,'Data_Type':kind,'Num_Elements':1,'Rec_Vary':True,'Dim_Sizes':[]},var_attrs=attr,var_data=data)
   w.close();c=cdflib.CDF(p);result,idx=holdout.read_selected(c,{'bins':[0,1,2]});
   with patch.object(c,'varattsget',side_effect=lambda n: {'UNITS':'s'} if n=='epoch_delta' else {}):
    with self.assertRaisesRegex(ValueError,'epoch_delta schema'):holdout.read_selected(c,{'bins':[0,1,2]})
   self.assertEqual(result['status'],'conditional_increase');self.assertEqual(idx.tolist(),list(range(15)));del c
   import gc;gc.collect()
 def test_holdout_authorization_prevents_cdf(self):
  args=['x','--manifest','x','--manifest-sha256','x','--authorization','x','--authorization-sha256','x','--prediction','x','--prediction-sha256','x','--output','D:/x']
  with patch.object(sys,'argv',args),patch.object(holdout,'verify'),patch.object(holdout,'sha',return_value='x'),patch.object(holdout.Path,'read_text',side_effect=[__import__('json').dumps(development.EXPECTED),__import__('json').dumps({'status':'selected','manifest_sha256':'x'})]),patch.object(holdout,'authorize',side_effect=ValueError('sealed')),patch.object(holdout.cdflib,'CDF',side_effect=AssertionError('opened')):
   with self.assertRaisesRegex(ValueError,'sealed'):holdout.main()
 def test_asymmetric_boundary_rejected(self):
  t=np.arange(0,43200,.5);b=np.zeros((len(t),3));b[:,0]=4;b[(t>=0)&(t<300),0]=10;b[(t>=300)&(t<600),0]=3.5
  self.assertLessEqual(3.5,.8*((10+4)/2));out=select(t,b,np.ones(len(t),bool));self.assertNotEqual(out.get('event_bin'),1);self.assertGreater(out['gates']['absolute_drop_ge_1nT'],0);self.assertEqual(len(out['bin_audit']),144)
 def test_swapi_coverage_audit(self):
  pred={'bins':[0,1,2]};t=np.array([60,120,180,360,420,660,720,780]);v=np.ones(len(t));out=pressure_test(t,np.ones(len(t))*10,v*5,v*1e5,v*.1,v*1000,np.ones(len(t),bool),pred)
  self.assertEqual(out['status'],'inconclusive_coverage');self.assertEqual(out['bin_audit'][1]['whole_window_contained'],2);self.assertEqual(out['bin_audit'][1]['reason'],'fewer_than_3_contained')
if __name__=='__main__':unittest.main()
