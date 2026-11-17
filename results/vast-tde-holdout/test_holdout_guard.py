import sys,unittest,tempfile,pathlib,json,hashlib
from unittest.mock import patch
import holdout_cli as h

class GuardTests(unittest.TestCase):
 def test_wrong_phase_before_snapshot(self):
  with patch.object(h,'snapshot',side_effect=AssertionError('snapshot reached')):
   with self.assertRaisesRegex(ValueError,'wrong phase'):h.verify_review('none','none','none','none','development')
 def test_unauthorized_before_any_read(self):
  class A:authorize_holdout=False
  with patch.object(h,'verify_review',side_effect=AssertionError('verification reached')):
   with self.assertRaisesRegex(ValueError,'authorization'):h.run(A())
 def test_wrong_input_digest_before_decoder(self):
  with patch.object(h,'sha',return_value='wrong'),patch.dict(sys.modules,{'pyarrow.parquet':None}):
   with self.assertRaisesRegex(ValueError,'input digest'):h.holdout_rows(1,{}, {},'fake','bad')
 def test_wrong_development_hash(self):
  with patch.object(h,'sha',return_value='wrong'):
   with self.assertRaisesRegex(ValueError,'original manifest mismatch'):h.verify_old('bad','bad')
 def test_wrong_baseline_file_digest(self):
  old={'git_head':h.OLD_HEAD,'worktree_sha256':{'results/vast-tde-expanded/conditional_pipeline.py':'code'}}
  summary={'manifest_sha256':h.OLD_MANIFEST_SHA,'git_head':h.OLD_HEAD,'sources':{k:{'status':'qualified-for-review','sha256':v,'file':k+'.json'} for k,v in h.DEV_SHA.items()}}
  with tempfile.TemporaryDirectory() as d:
   p=pathlib.Path(d);(p/'old.json').write_text(json.dumps(old));(p/'summary.json').write_text(json.dumps(summary))
   for k in h.DEV_SHA:(p/(k+'.json')).write_text('{}')
   def fake(path):
    s=str(path)
    if s.endswith('old.json'):return h.OLD_MANIFEST_SHA
    if s.endswith('summary.json'):return h.SUMMARY_SHA
    if 'conditional_pipeline.py' in s:return 'code'
    return 'wrong'
   with patch.object(h,'sha',side_effect=fake):
    with self.assertRaisesRegex(ValueError,'development file mismatch'):h.verify_old(p/'old.json',p)
if __name__=='__main__':unittest.main()
