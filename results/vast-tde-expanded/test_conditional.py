import unittest,tempfile,pathlib,numpy as np
from unittest.mock import patch
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
 def test_full_develop_freeze_evaluate(self):
  packet={'sources':{'1':{'development':self.groups[:4],'holdout':self.groups[4:]}}}
  dev=[r for r in self.rows if r['image_id']<=8];hold=[r for r in self.rows if r['image_id']>8]
  result=develop(1,packet,dev,self.controls)
  with tempfile.TemporaryDirectory() as d:
   f=pathlib.Path(d)/'dev.json';digest=freeze_development(f,result,'packet','code');frozen=load_reviewed_development(f,digest,'packet','code')
   self.assertFalse(evaluate_holdout(frozen,hold)['result']['both_at_least_3'])
   injected=[dict(r,flux_peak=3.) if r['source']==1 else r for r in hold]
   self.assertTrue(evaluate_holdout(frozen,injected)['result']['both_at_least_3'])
 def test_rms_and_median_standard_error(self):
  self.assertFalse(valid(dict(row(1,1),local_rms=1.821))[0])
  values=[row(s,1,10*(1+(s-112)*.04)) for s in range(100,125)]
  gains=image_gains(values,{s:10 for s in range(100,125)},[1]);x=np.array([v['flux_peak']/10 for v in values]);expected=max(.03,1.2533*1.4826*np.median(abs(x-np.median(x)))/(np.median(x)*np.sqrt(len(x))))
  self.assertAlmostEqual(gains[1]['fractional_sigma'],expected)
 def test_wrong_phase_prevents_reader(self):
  import development_cli as cli
  with patch.object(cli,'snapshot',side_effect=AssertionError('must not inspect environment')):
   with self.assertRaises(ValueError):cli.verify_manifest('none','none','holdout')
 def test_wrong_image_hash_precedes_parquet_decode(self):
  import read_development as reader
  with patch.object(reader,'guard',return_value={}),patch.object(reader,'digest',return_value='wrong'),patch.object(reader.pq,'read_table',side_effect=AssertionError('decoded')):
   with self.assertRaisesRegex(ValueError,'image metadata'):reader.read_development(1,'fake.parquet','wrong.csv')
 def test_worktree_blob_mismatch_rejected(self):
  import development_cli as cli
  with tempfile.TemporaryDirectory() as d:
   repo=pathlib.Path(d);(repo/'file').write_bytes(b'wrong\r\n')
   with patch.object(cli,'REPO',repo),patch.object(cli.subprocess,'check_output',return_value=b'wrong\n'):
    with self.assertRaisesRegex(ValueError,'worktree/Git blob bytes differ'):cli.verify_blob_bytes(['file'],'abc')
   (repo/'file').write_bytes(b'wrong\n')
   with patch.object(cli,'REPO',repo),patch.object(cli.subprocess,'check_output',return_value=b'wrong\n'):
    cli.verify_blob_bytes(['file'],'abc')
 def test_manifest_code_mismatch_prevents_execution(self):
  import development_cli as cli,hashlib,json
  with tempfile.TemporaryDirectory() as d:
   p=pathlib.Path(d)/'m.json';p.write_text(json.dumps({'git_head':'old'}));digest=hashlib.sha256(p.read_bytes()).hexdigest()
   with patch.object(cli,'snapshot',return_value={'git_head':'different'}):
    with self.assertRaisesRegex(ValueError,'code, packet'):cli.verify_manifest(p,digest,'development')
   with self.assertRaisesRegex(ValueError,'manifest SHA'):cli.verify_manifest(p,'0'*64,'development')
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
