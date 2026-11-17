import unittest,numpy as np
from plane import predict,proposal,sensitivity,compare
from development import empty_classification
class PlaneTest(unittest.TestCase):
 def setUp(self):
  self.t=np.arange(0,20,.01);self.v=np.tile([0.,0.,100.],(len(self.t),1));self.err=np.ones_like(self.v);self.good=np.ones(len(self.t),bool);self.pos=np.array([[0,0,0],[0,0,20],[0,0,40],[0,0,60.]])
  self.b=self.field(0.)
 def field(self,delay):
  a=np.clip(self.t-9.5-delay,0,1)*np.pi/2;return np.column_stack((8*np.cos(a),8*np.sin(a),.02*np.sin(self.t*5)))
 def test_full_plane_and_zero_control(self):
  x=proposal(self.t,self.b,self.v,self.good);self.assertIsNotNone(x);s=sensitivity(self.t,self.b,self.v,self.err,self.good,x,self.pos);self.assertEqual(s['status'],'qualified');lags=predict(x,self.pos)
  for i,lag in enumerate(lags):
   q=compare(self.t,self.field(lag),self.good,x,s['lag_low_s'][i],s['lag_high_s'][i]);self.assertEqual(q['status'],'matched')
  q=compare(self.t,self.b,self.good,x,s['lag_low_s'][0],s['lag_high_s'][0]);self.assertTrue(q['zero_lag_alternative'])
 def test_missing_and_low_speed(self):
  self.assertIsNone(proposal(self.t,self.b,self.v*.01,self.good));self.assertIsNone(proposal(self.t,self.b,self.v,np.zeros(len(self.t),bool)))
  x=proposal(self.t,self.b,self.v,self.good);s=sensitivity(self.t,self.b,self.v,np.ones_like(self.v)*100,self.good,x,self.pos);self.assertEqual(s['status'],'inconclusive_normal_speed')
 def test_constant_curved_wave(self):
  self.assertIsNone(proposal(self.t,np.ones_like(self.b)*5,self.v,self.good));x=proposal(self.t,self.b,self.v,self.good);s=sensitivity(self.t,self.b,self.v,self.err,self.good,x,self.pos)
  q=compare(self.t,self.field(5),self.good,x,s['lag_low_s'][0],s['lag_high_s'][0]);self.assertEqual(q['status'],'no_match')
  q=compare(self.t,self.field(0),self.good,x,-.3,.3);self.assertTrue(q['zero_lag_alternative'])
 def test_coverage_counts_missing_vs_constant(self):
  audit={};self.assertIsNone(proposal(self.t,self.b,self.v,np.zeros(len(self.t),bool),audit));self.assertEqual(audit['quality_cadence_windows'],0);self.assertGreater(audit['tested_1s_windows'],0);self.assertEqual(empty_classification(audit),'inconclusive_no_valid_coverage')
  audit={};self.assertIsNone(proposal(self.t,np.ones_like(self.b)*5,self.v,self.good,audit));self.assertGreater(audit['quality_cadence_windows'],0);self.assertGreater(audit['field_endpoint_windows'],0);self.assertEqual(audit['rotation_windows'],0);self.assertEqual(empty_classification(audit),'no_qualified_candidate')
 def test_spatial_curvature_four_craft(self):
  source=proposal(self.t,self.b,self.v,self.good);s=sensitivity(self.t,self.b,self.v,self.err,self.good,source,self.pos);lags=predict(source,self.pos)
  bent=[self.field(lags[0]),self.field(lags[1]),self.field(lags[2])]
  bent[1][:,1]*=-1
  verdicts=[compare(self.t,b,self.good,source,s['lag_low_s'][i],s['lag_high_s'][i])['status'] for i,b in enumerate(bent)]
  self.assertEqual(verdicts,['matched','no_match','matched'])
 def test_normal_sign_and_coverage(self):
  p={'normal':[1.,0.,0.],'normal_speed_km_s':50.};r=np.array([[0,0,0],[50,0,0],[-100,0,0],[25,0,0]]);np.testing.assert_allclose(predict(p,r),[1,-2,.5]);p['normal']=[-1.,0.,0.];p['normal_speed_km_s']=-50.;np.testing.assert_allclose(predict(p,r),[1,-2,.5]);x=proposal(self.t,self.b,self.v,self.good);self.assertEqual(compare(self.t,self.b,self.good,x,-20,0)['status'],'inconclusive_coverage')
if __name__=='__main__':unittest.main()
