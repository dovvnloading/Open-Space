import json,os,tempfile,unittest
from pathlib import Path
from replay_guard import fresh,verified,sha,identity,check_checkpoint
class GuardTests(unittest.TestCase):
 def test_raw_and_checkpoint_binding(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp)/'raw';p.write_bytes(b'one');expected=sha(p);verified(p,expected)
   ident=identity(p,'event',True,{'model':'abc'},{'python':'3'}, {'t0':1.0});self.assertEqual(check_checkpoint({'_identity':ident},ident)['_identity'],ident)
   for key,value in [('cap',False),('code_sha256',{'model':'changed'}),('parameters',{'t0':2.0})]:
    altered=dict(ident);altered[key]=value
    with self.assertRaises(ValueError):check_checkpoint({'_identity':ident},altered)
   p.write_bytes(b'two')
   with self.assertRaises(ValueError):verified(p,expected)
 def test_real_consolidation_checkpoint_identity(self):
  from unittest.mock import patch
  import consolidate
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);code=root/'code';code.write_bytes(b'original');stage=root/'stage';stage.mkdir()
   spec={'KMT-2024-BLG-0001':{'input_sha256':'rawhash'}};code_hashes={str(code):sha(code)};runtime={'python':'test'}
   meta={'cap':True,'code_sha256':code_hashes,'runtime':runtime};(stage/'campaign.json').write_text(json.dumps(meta))
   ident={'event':'KMT-2024-BLG-0001','input_sha256':'rawhash','cap':True,'code_sha256':code_hashes,'runtime':runtime,'parameters':{'t0':1}}
   item=stage/'KMT-2024-BLG-0001.json';item.write_text(json.dumps({'status':'ok','_identity':ident}))
   with patch.object(consolidate,'external',side_effect=lambda p:Path(p)):
    self.assertEqual(len(consolidate.read_results(stage,spec,True)),1)
    with self.assertRaises(ValueError):consolidate.read_results(stage,spec,False)
    ident['input_sha256']='changed';item.write_text(json.dumps({'status':'ok','_identity':ident}))
    with self.assertRaises(ValueError):consolidate.read_results(stage,spec,True)
    ident['input_sha256']='rawhash';item.write_text(json.dumps({'status':'ok','_identity':ident}));code.write_bytes(b'changed')
    with self.assertRaises(ValueError):consolidate.read_results(stage,spec,True)
 def test_real_entrypoints_reject_cached_or_corrupt_data(self):
  import importlib.util
  from unittest.mock import patch
  root=Path(__file__).resolve().parent
  def module(name,path):
   spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
  kmt=module('guarded_kmt',root/'kmt'/'freeze.py');other=module('guarded_other',root/'other'/'fetch.py');ogle=module('guarded_ogle',root/'ogle'/'freeze.py')
  with tempfile.TemporaryDirectory() as tmp:
   d=Path(tmp);(d/'KMT-2024-BLG-0001.tar.gz').write_bytes(b'stale')
   with self.assertRaises(FileExistsError):kmt.fetch({'event':'KMT-2024-BLG-0001','url':'unused'},d)
   (d/'records').mkdir();(d/'curves').mkdir();(d/'records'/'MOA-2024-BLG-001.json').write_text('{}')
   with self.assertRaises(FileExistsError):other.process({'event':'MOA-2024-BLG-001','page':'unused'},d)
   (d/'2025').mkdir();(d/'2025'/'blg-0001.dat').write_bytes(b'bad')
   ogle.FROZEN_CURVES={'2025-BLG-0001':{'sha256':'0'*64,'bytes':'3'}}
   with patch.object(ogle,'get',side_effect=AssertionError('network requested')):
    result=ogle.fetch_curve(d,{'event':'2025-BLG-0001','t0_hjd':'2450000'})
   self.assertIn('frozen input mismatch',result['status'])
 @unittest.skipUnless(os.name=="nt", "requires external Windows D: temporary root")
 def test_exclusive_normalized_paths(self):
  with tempfile.TemporaryDirectory() as tmp:
   dest=fresh(Path(tmp)/'one'/'..'/'out');self.assertTrue(dest.is_dir())
   with self.assertRaises(FileExistsError):fresh(Path(tmp)/'out')
   with self.assertRaises(ValueError):fresh(Path('relative'))
if __name__=='__main__':unittest.main()
