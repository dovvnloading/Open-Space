import importlib.util,json,os,sys,tempfile,types,unittest
import numpy
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
  import consolidate,replay_guard
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);stage=root/'stage';stage.mkdir();pkg=root/'VBMicrolensing';(pkg/'data').mkdir(parents=True)
   package=pkg/'__init__.py';native=pkg/'VBMicrolensing.pyd';table=pkg/'data'/'ESPL.tbl';lib=root/'grid.so';code=root/'code.py'
   for path in (package,native,table,lib,code):path.write_bytes(path.name.encode())
   sources={'run_grid':code};roles={'run_grid':str(code),'vbm_package':str(package),'vbm_native':str(native),'espl_table':str(table),'grid_library':str(lib)}
   spec={'KMT-2024-BLG-0001':{'input_sha256':'rawhash'}};code_hashes={str(p):sha(p) for p in (code,package,native,table,lib)};runtime={'python':'3','platform':'test','numpy':'1','vbm':'5.6.1'}
   timing={'KMT-2024-BLG-0001':{'t0':1}}
   meta={'schema':replay_guard.CAMPAIGN_SCHEMA,'model':'frozen grid','cap':True,'asset_roles':roles,'code_sha256':code_hashes,'runtime':runtime,'parameters_sha256':consolidate.parameters_digest(timing)}
   campaign=stage/'campaign.json';campaign.write_text(json.dumps(meta));ident={'event':'KMT-2024-BLG-0001','input_sha256':'rawhash','cap':True,'code_sha256':code_hashes,'runtime':runtime,'parameters':{'t0':1}}
   item=stage/'KMT-2024-BLG-0001.json';item.write_text(json.dumps({'status':'ok','_identity':ident}))
   with patch.object(consolidate,'external',side_effect=lambda p:Path(p)),patch.object(consolidate,'frozen_parameters',return_value=timing),patch.object(replay_guard,'source_assets',return_value=sources):
    self.assertEqual(len(consolidate.read_results(stage,spec,True)),1)
    for key,value in [('code_sha256',{}),('runtime',{}),('asset_roles',{}),('schema','old')]:
     bad=dict(meta);bad[key]=value;campaign.write_text(json.dumps(bad))
     with self.assertRaises(ValueError):consolidate.read_results(stage,spec,True)
    for role in ('vbm_package','vbm_native','espl_table','grid_library'):
     bad=dict(meta);bad['asset_roles']=dict(roles);del bad['asset_roles'][role];campaign.write_text(json.dumps(bad))
     with self.assertRaises(ValueError):consolidate.read_results(stage,spec,True)
    campaign.write_text(json.dumps(meta))
    with self.assertRaises(ValueError):consolidate.read_results(stage,spec,False)
    ident['input_sha256']='changed';item.write_text(json.dumps({'status':'ok','_identity':ident}))
    with self.assertRaises(ValueError):consolidate.read_results(stage,spec,True)
    ident['input_sha256']='rawhash';ident['parameters']={'t0':2};item.write_text(json.dumps({'status':'ok','_identity':ident}))
    with self.assertRaises(ValueError):consolidate.read_results(stage,spec,True)
    ident['parameters']={'t0':1};item.write_text(json.dumps({'status':'ok','_identity':ident}));table.write_bytes(b'changed table')
    with self.assertRaises(ValueError):consolidate.read_results(stage,spec,True)
    table.write_bytes(b'ESPL.tbl');code.write_bytes(b'changed')
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
 def test_duplicate_jobs_rejected_before_output(self):
  from unittest.mock import patch
  root=Path(__file__).resolve().parent;spec=importlib.util.spec_from_file_location('guarded_run_grid',root/'run_grid.py');module=importlib.util.module_from_spec(spec)
  with patch.dict(sys.modules,{'fit':types.ModuleType('fit')}):spec.loader.exec_module(module)
  duplicate=("KMT-2024-BLG-0001",Path('unused'),1.,1.)
  with tempfile.TemporaryDirectory() as tmp:
   out=Path(tmp)/'out'
   with patch.object(module,'external',side_effect=lambda p:Path(p)),patch.object(module,'verify_receipt'),patch.object(module,'inputs',return_value=[duplicate,duplicate]),patch.object(module,'fresh',side_effect=AssertionError('output created')),patch.object(sys,'argv',['run_grid.py','--cache',tmp,'--checkpoint',str(out),'--grid-library','unused','--kmt-bound']):
    with self.assertRaisesRegex(ValueError,'duplicate'):module.main()
   self.assertFalse(out.exists())
 def test_engine_assets_include_native_and_numerical_table(self):
  from unittest.mock import patch
  from replay_guard import engine_assets
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);(root/'data').mkdir();pkg=root/'__init__.py';native=root/'engine.pyd';table=root/'data'/'ESPL.tbl'
   for p in (pkg,native,table):p.write_bytes(p.name.encode())
   module=types.ModuleType('mock_vbm_native');module.__file__=str(native);klass=type('Engine',(),{'__module__':'mock_vbm_native'});vbm=types.SimpleNamespace(__file__=str(pkg),VBMicrolensing=klass)
   with patch.dict(sys.modules,{'mock_vbm_native':module}):
    assets=engine_assets(vbm);self.assertEqual(set(assets),{pkg,native,table});old={str(p):sha(p) for p in assets};table.write_bytes(b'changed');self.assertNotEqual(old[str(table)],sha(table));native.unlink()
    with self.assertRaises(ValueError):engine_assets(vbm)
 @unittest.skipUnless(os.name=="nt", "requires external Windows D: temporary root")
 def test_exclusive_normalized_paths(self):
  with tempfile.TemporaryDirectory() as tmp:
   dest=fresh(Path(tmp)/'one'/'..'/'out');self.assertTrue(dest.is_dir())
   with self.assertRaises(FileExistsError):fresh(Path(tmp)/'out')
   with self.assertRaises(ValueError):fresh(Path('relative'))
if __name__=='__main__':unittest.main()
