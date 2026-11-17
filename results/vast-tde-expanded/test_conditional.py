import unittest,tempfile,pathlib,numpy as np
from conditional_pipeline import *

def row(source,image,value=10.,error=.1,forced=False):
 return dict(source=source,image_id=image,flux_peak=value,flux_peak_err=error,local_rms=.1,beam_bmaj=.004,beam_bmin=.003,forced=forced,flag_c4=False,has_siblings=False)

class TestConditional(unittest.TestCase):
 def setUp(self):
  self.groups=[[1,2],[3,4],[5,6],[7,8],[9,10],[11,12]];self.controls=set(range(100,125));self.target=1
  self.rows=[row(s,i,10+s*.01) for s in self.controls for i in range(1,13)]+[row(1,i,-.5 if i==1 else 1.,.2,True) for i in range(1,13)]
 def test_full_dev_and_two_session_holdout(self):
  allowed=set(range(1,9));dev,drops=phase_rows([r for r in self.rows if r['image_id'] in allowed],allowed,self.controls|{1});self.assertEqual(len(drops),0)
  selected=stable_controls(dev,self.controls,self.groups[:4]);self.assertEqual(len(selected),25)
  g=image_gains(dev,selected,allowed);v=session_measurements(dev,1,self.groups[:4],g);self.assertEqual(len(v),4);self.assertLess(v[0][0],1)
  pred=baseline(v);future=set(range(9,13));held,_=phase_rows([r for r in self.rows if r['image_id'] in future],future,self.controls|{1});hg=image_gains(held,selected,future);h=session_measurements(held,1,self.groups[-2:],hg);self.assertEqual(len(h),2);self.assertFalse(compare(pred,h)['both_at_least_3'])
 def test_forced_signed_and_invalid(self):
  self.assertTrue(valid(row(1,1,-3.,.2,True))[0]);bad=row(1,1);bad['local_rms']=float('nan');self.assertFalse(valid(bad)[0]);del bad['beam_bmin'];self.assertFalse(valid(bad)[0]);bad=row(1,1);bad['flag_c4']=True;self.assertFalse(valid(bad)[0])
 def test_phase_duplicate_and_missing_control(self):
  with self.assertRaises(ValueError):phase_rows([row(1,9)],{1},{1})
  with self.assertRaises(ValueError):phase_rows([row(1,1),row(1,1)],{1},{1})
  with self.assertRaises(ValueError):image_gains([row(100,1)],{100:11},[1])
 def test_correlated_pairs_not_independent(self):
  g={1:dict(gain=1.,fractional_sigma=.03),2:dict(gain=1.,fractional_sigma=.03)};v=session_measurements([row(1,1,1,.2),row(1,2,1,.2)],1,[[1,2]],g)[0];self.assertGreaterEqual(v[1],.2**2)
 def test_frozen_development_identity(self):
  with tempfile.TemporaryDirectory() as d:
   p=pathlib.Path(d)/'baseline.json';digest=freeze_development(p,{'source':1},'packet','code')
   self.assertEqual(load_reviewed_development(p,digest,'packet','code'),{'source':1})
   with self.assertRaises(FileExistsError):freeze_development(p,{},'packet','code')
   with self.assertRaises(ValueError):load_reviewed_development(p,digest,'changed','code')
   p.write_text('{}')
   with self.assertRaises(ValueError):load_reviewed_development(p,digest,'packet','code')
 def test_toy_null_and_injection(self):
  rng=np.random.default_rng(708);n=50000;e=rng.normal(size=(n,6));h=rng.normal(size=(n,2));z=(h-e.mean(axis=1)[:,None])/np.sqrt(1+1/6);self.assertLess(np.mean((z>=3).all(axis=1)),.002);self.assertGreater(np.mean((z+4/np.sqrt(1+1/6)>=3).all(axis=1)),.4)
if __name__=='__main__':unittest.main()
