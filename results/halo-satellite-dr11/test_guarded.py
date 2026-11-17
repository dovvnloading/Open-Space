import hashlib,json,pathlib,tempfile,unittest
from unittest.mock import patch
import run_guarded as g
class GuardTests(unittest.TestCase):
 def setUp(self):
  self.t=tempfile.TemporaryDirectory(dir='D:/Open-Space/campaigns');self.r=pathlib.Path(self.t.name);(self.r/'data').write_bytes(b'fixed');self.m={'files':{'data':g.sha(self.r/'data')},'runtime':{}}
 def tearDown(self): self.t.cleanup()
 def test_modified_input_rejected(self):
  (self.r/'data').write_bytes(b'changed')
  with self.assertRaisesRegex(ValueError,'input/code mismatch'):g.verify(self.r,self.m)
 def test_changed_runtime_rejected(self):
  with self.assertRaisesRegex(ValueError,'runtime mismatch'):g.verify(self.r,self.m)
 def test_collision_does_not_overwrite(self):
  source=self.r/'source';source.mkdir();out=self.r/'exists';out.mkdir();(out/'keep').write_text('original')
  with patch.object(g,'verify'):
   with self.assertRaises(FileExistsError):g.create(source,out,self.m,'search')
  self.assertEqual((out/'keep').read_text(),'original')
 def test_relative_and_source_output_denied(self):
  with self.assertRaises(ValueError):g.create(self.r,'relative',self.m,'search')
  with self.assertRaises(ValueError):g.create(self.r,self.r/'child',self.m,'search')
if __name__=='__main__':unittest.main()
